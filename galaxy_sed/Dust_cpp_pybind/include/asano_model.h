#ifndef ASANO_MODEL_H
#define ASANO_MODEL_H

#include <cfloat>
#include <iostream>
#include <chrono>

#include <constants.h>
#include <file_util.h>
#include <SED_setting_util.h>
#include <calculation_util.h>

#include "function.h"
#include "read_file.h"

namespace asano_model {
using val = std::valarray<double>; // doubleの配列を表す
using val2 = std::valarray<val>; // valの配列なので、2次元配列である
using val3 = std::valarray<val2>;

/**
 int main() {
    val array1 = {1.0, 2.0, 3.0};  // double の配列
    val2 array2 = {array1, {4.0, 5.0, 6.0}};  // val の配列（2次元配列）

    // ここで array1 や array2 を使用するなどの処理が続く...

    return 0;
}
　のようにして型名を簡潔にし、配列を表すのに使う。
*/

using SEDSetting = my_util::SED_util::SEDSettingUtil;

/**
 * Calculate dust evolution with Asano model (Asano et al., 2013a, b, 2014; Nozawa et al., 2014)
 */
class AsanoModel {
  private:
    const val  a_val_         = SEDSetting::DustRadiusDustModel();
    const val3 M_SN_val3_     = ReadSNDustFile(a_val_); // Dust mass produced by supernova
    const val4 M_AGB_val4_    = ReadAGBDustFile(); // Dust mass produced by AGB
    const val3 Y_Z_val3_      = ReadMetalYieldFile(); // Metal yield from stars
    const val2 M_rem_val2_    = ReadRemnantFile(); // Stellar remnant mass
    const val4 dest_eta_val4_ = ReadSNDestructionFile(); // Dust destroyed efficiency by SN shocks
    const val  tau_star_val_  = StellarLifeTime();

  public:
    AsanoModel() = default;

    ~AsanoModel() = default;

    /**
     * Calculate dust evolution
     * @param free_params
     */
    inline void Calculate(const FreeParameter& free_params,
                          const std::string& ofn_total_dust_mass,
                          const std::string& ofn_dust_mass_distribution,
                          const std::string& ofn_dust_number_distribution) const {
        const auto n_age = static_cast<size_t>(free_params.age_max_ / TIME_BIN_DUST_MODEL) + 1;

        const auto calc_start = std::chrono::system_clock::now();

        const auto IMF_val = SFH::NormalizedIMF();

        auto M_gas_val = val(n_age);
        if (free_params.is_infall_) M_gas_val[0] = DBL_MIN; // If the initial mass is 0, it will be troublesome. Set it to 1.
        else M_gas_val[0] = 1.0;
//        else M_gas_val[0] = free_params.M_total_gal_;

        auto M_star_val    = val(n_age);
        auto M_Z_val       = val(n_age);
        auto M_C_ISM_val   = val(n_age); // Carbon mass in ISM (including non dust)
        auto M_sil_ISM_val = val(n_age); // Silicate mass in ISM (including non dust)
        auto SNR_val       = val(n_age); // Supernova rate
        auto Y_d_val3      = val3(val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES),
                                  n_age); // Dust yield by stars

        for (auto i_age = std::size_t(0); i_age < n_age - 1; ++i_age) {
            const auto i_m_min      = MinimumStellarMassIndex(i_age, tau_star_val_);
            const auto i_birth_valt = StellarBirthIndex(i_age, i_m_min, tau_star_val_); // (t - tau)

            const auto i_Z_valt = StellarMetallicityIndex(i_m_min, i_birth_valt,
                                                          M_Z_val / M_gas_val);
            const auto SFR_val  = SFH::SFRValarray(i_age, free_params.schmidt_index_,
                                                   free_params.tau_SF_, M_gas_val);

            // Equ. (4.7) of Nishida D-thesis
            const auto R = ReturnMass(i_m_min, i_birth_valt, i_Z_valt, IMF_val, SFR_val,
                                      M_rem_val2_);

            // Equ. (4.9) of Nishida D-thesis
            const auto[Y_Z, Y_C, Y_sil] = MetalYield(i_m_min, i_birth_valt, i_Z_valt, IMF_val,
                                                     SFR_val, Y_Z_val3_);

            // Equ. (4.10) of Nishida D-thesis
            Y_d_val3[i_age + 1] = DustYield(i_m_min, i_birth_valt, M_gas_val, M_Z_val, IMF_val,
                                            SFR_val, M_SN_val3_, M_AGB_val4_);
            SNR_val[i_age + 1]  = SNR(i_m_min, i_birth_valt, IMF_val, SFR_val);

            // Equ. (4.1) of Nishida D-thesis
            M_star_val[i_age + 1] = StellarMass(M_star_val[i_age], M_gas_val[i_age], R,
                                                free_params.schmidt_index_, free_params.tau_SF_);

            // Equ. (4.2) of Nishida D-thesis
            M_gas_val[i_age + 1] = GasMass(free_params.is_infall_, i_age, M_gas_val[i_age], R,
                                           free_params.schmidt_index_, free_params.tau_SF_,
                                           free_params.tau_infall_);

            // Equ. (4.3) of Nishida D-thesis
            M_Z_val[i_age + 1]       = MetalMass(M_gas_val[i_age], SFR_val[i_age],
                                                 M_Z_val[i_age],
                                                 Y_Z);
            M_C_ISM_val[i_age + 1]   = MetalMass(M_gas_val[i_age], SFR_val[i_age],
                                                 M_C_ISM_val[i_age],
                                                 Y_C);
            M_sil_ISM_val[i_age + 1] = MetalMass(M_gas_val[i_age],
                                                 SFR_val[i_age], M_sil_ISM_val[i_age], Y_sil);

        }

        my_util::MyTime::PrintElapsedTime(calc_start, "sec",
                                          "Asano model (calc before dust evolution): ");

        const auto volume_val = Map(a_val_, Volume);

        // Make output file stream
        auto ofp_total = MakeTotalMassFile(ofn_total_dust_mass);
        auto ofp_m     = MakeMassDistributionFile(ofn_dust_mass_distribution);
        auto ofp_n     = MakeNumberDistributionFile(ofn_dust_number_distribution);

        auto M_dust_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);

        for (auto i_age = std::size_t(0); i_age < n_age; ++i_age) {
            if (i_age % static_cast<int>(TIME_BIN / TIME_BIN_DUST_MODEL) == 0 && i_age != 0) {
                WriteFiles(i_age, M_Z_val[i_age], M_gas_val[i_age], M_C_ISM_val[i_age],
                           M_sil_ISM_val[i_age], M_star_val[i_age], a_val_, volume_val,
                           M_dust_val2, ofp_total, ofp_n, ofp_m);
            }
            if (i_age == n_age - 1) break;

            const auto M_evolution_val2 = DustMassEvolution(M_gas_val[i_age],
                                                            {M_C_ISM_val[i_age],
                                                             M_sil_ISM_val[i_age]},
                                                            a_val_, volume_val, M_dust_val2);

            const auto M_SND_val2 = SNDestruction(M_gas_val[i_age], M_Z_val[i_age], SNR_val[i_age],
                                                  a_val_, M_dust_val2, dest_eta_val4_);

            M_dust_val2 = DustMass(M_gas_val[i_age], free_params.schmidt_index_,
                                   free_params.tau_SF_, M_dust_val2, Y_d_val3[i_age], M_SND_val2,
                                   M_evolution_val2);
        }

        fclose(ofp_total);
        fclose(ofp_n);
        fclose(ofp_m);

        my_util::MyTime::PrintElapsedTime(calc_start, "sec", "Asano model (calc): ");
    }
};


}
#endif // ASANO_MODEL_H
