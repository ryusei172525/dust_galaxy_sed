#ifndef SED_MODEL_MAIN_H
#define SED_MODEL_MAIN_H

#include <cstdio>
#include <cstring>

#include <file_util.h>
#include <SED_setting_util.h>
#include <utility>
#include <valarray>
#include <calculation_util.h>
#include <future>

#include "read.h"
#include "spectra.h"

namespace stellar_spectrum {
using vec = std::vector<double>;
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;
using valt = std::valarray<std::size_t>;

using SEDSetting = my_util::SED_util::SEDSettingUtil;

class StellarSpectrum {
  private:
    static constexpr auto N_MAX_TIME_    = std::size_t(20002);
    static constexpr auto N_MAX_TIME_SSP = std::size_t(600);
    static constexpr auto n_lambda_      = std::size_t(1221);

    const std::vector<std::string> fn_track_vec_ = ReadSSPZFileName(); // track の名前が書いてあるファイルからその名前を読み込む 例. Chabrier_SSPs.dat, Chabrier_x=0.5_SSPs.datなど
    // file 名から金属量を読み取る
    const std::vector<double>      Z_SSP_vec_    = ReadMetallicityVector(fn_track_vec_);
    const std::size_t              n_Z_          = fn_track_vec_.size(); // The number of metallicity

    val lambda_angstrom_val_;
    val lambda_cm_output_val_;

    valt t_inv_val_     = valt(N_MAX_TIME_);
    val2 ejecta_val2_   = val2(val(N_MAX_TIME_), n_Z_);
    val2 ejecta_Z_val2_ = val2(val(N_MAX_TIME_), n_Z_);
    val2 m_BHNS_val2_   = val2(val(N_MAX_TIME_), n_Z_);
    val2 m_WD_val2_     = val2(val(N_MAX_TIME_), n_Z_);
    val2 m_alive_val2_  = val2(val(N_MAX_TIME_), n_Z_);
    val  F_neb_val      = val(n_lambda_); // Flux of nebular continuum
    val  beta_val_      = val(N_MAX_TIME_); // t_SSP bin での右側の値
    val2 n_Lym_val2_    = val2(val(N_MAX_TIME_), n_Z_);
    val3 F_SSP_val3     = val3(val2(val(n_lambda_), N_MAX_TIME_SSP), n_Z_);

    /**
     * Read stellar spectrum, nebular spectrum, and track files.
     */
    void ReadFiles() {
        const auto start_time = my_util::MyTime::GetTime();

        auto F_stellar_LCB_val2 = val2(); // Flux of stellar (< 5000 K)
        auto F_stellar_CM_val2  = val2(); // Flux of stellar (> 5000 K)
        ReadStellarSpectrumFile(lambda_angstrom_val_, F_stellar_LCB_val2, F_stellar_CM_val2);

        auto lambda_line_val = val(); // wavelength of line unit in angstrom
        auto F_line_val      = val(); // Flux of nebular lines
        ReadNebularSpectrumFile(lambda_angstrom_val_, F_neb_val, lambda_line_val, F_line_val);

        auto n_spitzer            = std::size_t(0);
        auto tau_dust_Spitzer_val = val(20);
        auto y_Spitzer_val        = val(20);
        ReadSpitzer(n_spitzer, tau_dust_Spitzer_val, y_Spitzer_val);

        // Total luminosity with a metallicity Z, time t, and wavelength lambda.
        auto t_SSP_valt = valt(N_MAX_TIME_SSP); // Almost same as index

        auto F_bol_SSP_val2 = val2(val(N_MAX_TIME_SSP), n_Z_);

        auto n_SNII_val2 = val2(val(N_MAX_TIME_), n_Z_);
        auto n_SNIa_val2 = val2(val(N_MAX_TIME_), n_Z_);

        // PEGASEで作った星進化トラックdatファイルの読み込み
        ReadStellarFluxFile(fn_track_vec_, F_stellar_LCB_val2, F_stellar_CM_val2, t_SSP_valt,
                            F_bol_SSP_val2, F_SSP_val3, m_alive_val2_, t_inv_val_, beta_val_,
                            n_Lym_val2_, n_SNII_val2, n_SNIa_val2, ejecta_val2_, ejecta_Z_val2_,
                            m_BHNS_val2_, m_WD_val2_);

        std::cout << "Stellar spectrum (read files): ";
        my_util::MyTime::PrintElapsedTime(start_time, "sec");
    }

  public:
    /**
     * Set wavelength member variable and read files
     * @param lambda_cm_output_val
     */
    inline explicit StellarSpectrum(val lambda_cm_output_val) :
    lambda_cm_output_val_(std::move(lambda_cm_output_val)) {
        ReadFiles();
    }

    /**
     * Calculate stellar luminosity
     * @param free_params
     */
    inline void Calculate(const FreeParameter& free_params) const {
        const auto start_time = my_util::MyTime::GetTime();

        auto Z_gas_val = val(N_MAX_TIME_); // Metallicity of gas
        Z_gas_val[0] = INITIAL_Z;
        auto Z_SFR_val = val(N_MAX_TIME_); // Metallicity of formed star
        Z_SFR_val[0] = INITIAL_Z;

        // Initialization
        auto M_gal_val = val(N_MAX_TIME_); // Unit in [/Msys]
        if (free_params.is_infall_) M_gal_val[0] = DBL_MIN;
//        else M_gal_val[0] = free_params.M_total_gal_;
        else M_gal_val[0] = 1.0;

        auto sigma_gas_val = val(N_MAX_TIME_); // Mass of gas
        sigma_gas_val[0] = M_gal_val[0];

        double sigma_Z = sigma_gas_val[0] * Z_gas_val[0]; // Mass of metal in gas phase

        auto sigma_BHNS_val = val(N_MAX_TIME_);
        auto sigma_WD_val   = val(N_MAX_TIME_);
        auto sigma_sub_val  = val(N_MAX_TIME_);

        auto ofp = my_util::FileUtil::OfpOpen(StellarContinuumFileName(free_params));
        fprintf(ofp, "wavelength[cm] L_gal[erg/s/cm]\n");
        auto ofs = my_util::FileUtil::OfsOpen(FYoungFileName(free_params));
        ofs << "wavelength[cm] f_y[t_y=10Myr] f_y[t_y=100Myr] f_y[t_y=1Gyr]" << std::endl;

        auto SFR_lum_val    = val(N_MAX_TIME_);
        auto age_star_val   = val(N_MAX_TIME_); // Averaged stellar age
        auto sigma_star_val = val(N_MAX_TIME_);
        auto Z_star_val     = val(N_MAX_TIME_);
        auto SFR_val        = val(N_MAX_TIME_);

        // 各 age での Z_SFR が track file の metallicity を超える最初のindexが i_Z + 1
        // example: Z_SFR = 2.5, Z_SSP = {1, 2, 3, 4} -> i_Z = 1
        auto i_Z_val   = valt(N_MAX_TIME_);
        auto alpha_val = val(N_MAX_TIME_); // alpha は i_Z にどれだけ Z_SFR が近いかを表す

        MassEvolution(n_Z_, free_params.age_max_, free_params.schmidt_index_,
                      free_params.is_infall_, free_params.tau_infall_, free_params.tau_SF_,
                      sigma_Z, Z_SSP_vec_, ejecta_val2_, ejecta_Z_val2_, m_BHNS_val2_, m_WD_val2_,
                      m_alive_val2_, i_Z_val, SFR_val, sigma_gas_val, Z_SFR_val, alpha_val,
                      Z_gas_val, SFR_lum_val, age_star_val, sigma_star_val, Z_star_val,
                      sigma_sub_val, M_gal_val, sigma_BHNS_val, sigma_WD_val);
        //auto i0 = std::size_t(0);

        /////// Calculate stellar luminosity ///////////////////////////////////////////
        // Unit in [Myr]
        const auto age = (free_params.ai_ + 1) * static_cast<std::size_t>(TIME_BIN / 1e6);

        const auto f_atte_nebular = 1.0;
        // 宇宙年齢がtのときにできた星による宇宙年齢ageでの放射
        const auto F_gal_val2     = ContinuumFluxFromBornAtTStar(n_lambda_, age, f_atte_nebular,
                                                                 i_Z_val, t_inv_val_, alpha_val,
                                                                 beta_val_, F_neb_val, SFR_lum_val,
                                                                 n_Lym_val2_, F_SSP_val3);

        const auto F_gal_val = ContinuumFlux(n_lambda_, F_gal_val2);

        // BolometricFlux(age, i_Z_val, alpha_val, beta_val_, t_inv_val_, SFR_lum_val, F_bol_SSP_val2, Z_SFR_val);

        // const auto F_line_tot_val = LineFlux(age, F_line_val, SFR_lum_val, alpha_val, beta_val_,
        //                                 n_Lym_val2_, f_atte_nebular, i_Z_val, t_inv_val_);

        // [[maybe_unused]] const auto n_Lym_tot = NumberOfLymanContinuumPhoton(age, i_Z_val,            //                                                                 SFR_lum_val, alpha_val,                                                                 n_Lym_val2_);

        // [[maybe_unused]] const auto n_SNII_tot = NumberOfSNII(age, i_Z_val, SFR_lum_val, alpha_val, n_SNII_val2);

        // [[maybe_unused]] const auto n_SNIa_tot = NumberOfSNIa(age, i_Z_val, SFR_lum_val, alpha_val, n_SNIa_val2);

        WriteContinuumFile(lambda_angstrom_val_ * 1e-8, F_gal_val, lambda_cm_output_val_, ofp);

        const auto f_young_star_val2 = YoungStellarFraction(n_lambda_, lambda_angstrom_val_,
                                                            age, F_gal_val2, F_gal_val);

        WriteYoungStellarFractionFile(lambda_cm_output_val_, f_young_star_val2,
                                      lambda_angstrom_val_ * 1e-8, ofs);
        fclose(ofp);

        std::cout << "Stellar spectrum (calc): ";
        my_util::MyTime::PrintElapsedTime(start_time, "sec");
    }
};
} // namespace
#endif //SED_MODEL_MAIN_H
