#ifndef DUST_PARAMETER_H
#define DUST_PARAMETER_H

#include <valarray>

#include <calculation_util.h>
#include <SED_file_util.h>
#include <SED_setting_util.h>

namespace dust_property {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using SEDFile = my_util::SED_util::SEDFileUtil;

/**
 * Read dust property file each dust species
 * @param param_name e.g., G, Albedo, Qabs, Qext
 * @param dust_species_vector e.g., Sil, Gra, PAHneu, and PAHion.
 * @return val3 [dust species][wavelength][dust radius]
 */
val3 ReadDustPropertyFiles(const std::vector<std::string>& dust_species_vector,
                           const std::string& param_name) {
    auto param_val3 = val3(dust_species_vector.size());

    for (auto i = std::size_t(0); i < dust_species_vector.size(); ++i) {
        param_val3[i] = SEDFile::ReadDustPropertyFile(dust_species_vector[i], param_name);
    }
    return param_val3;
}

/**
 * Write averaged dust properties file
 * @param ofn Output file name
 * @param n_age Number of age
 * @param lambda_cm_val Wavelength [cm]
 * @param k_val Extinction coefficient per dust mass [g^-1]
 * @param omega_val Scattering albedo
 * @param g_val asymmetry parameter
 */
void WriteAveragedPropertyFile(const std::string& ofn, std::size_t n_age, const val& lambda_cm_val,
                               const val2& k_val, const val2& omega_val, const val2& g_val,
                               const val2& Cabs_sum_val2) {
    auto ofs = my_util::FileUtil::OfsOpen(ofn);
    ofs << "wavelength[cm] k[g^-1] omega g" << std::endl;

    for (auto i_age = std::size_t(0); i_age < n_age; ++i_age) {
        for (auto i_lambda = std::size_t(0); i_lambda < N_LAMBDA; ++i_lambda) {
            ofs << lambda_cm_val[i_lambda] << " " << k_val[i_age][i_lambda]
                << " " << omega_val[i_age][i_lambda] << " " << g_val[i_age][i_lambda]
                << " " << Cabs_sum_val2[i_age][i_lambda] << std::endl;
        }
    }
}

/**
 * This averaged dust parameters dose not depend on total galaxy mass
 */
class AveragedDustProperty {
    using SEDSetting = my_util::SED_util::SEDSettingUtil;
  public:
    /**
     * Calculate averaged dust properties.
     * @param a_cm_val Dust radius
     * @param free_params Free parameter class object
     * @param n_val3 Dust number distribution. [dust species][galaxy age][dust radius]
     * @return val3 Extinction coefficient in unit mass [Msun/g], scattering albedo, and asymmetry parameter
     * [param species][galaxy age][wavelength]
     * Return value dose not depend on total galaxy mass
     */
    inline static val3
    Calculate(const val& a_cm_val, const FreeParameter& free_params, const val3& n_val3,
              const std::string& ofn) {
        // Read extinction coefficient per unit dust mass, albedo, and asymmetry parameter from file each dust species
        const auto dust_species_vector = std::vector<std::string>{"Sil", "Gra", "PAHneu", "PAHion"};

        const auto g_val3    = ReadDustPropertyFiles(dust_species_vector, "G");
        const auto Qabs_val3 = ReadDustPropertyFiles(dust_species_vector, "Qabs");
        const auto Qext_val3 = ReadDustPropertyFiles(dust_species_vector, "Qext");

        auto      Qsca_val3 = val3(dust_species_vector.size());
        for (auto i_species = std::size_t(0); i_species < dust_species_vector.size(); ++i_species) {
            Qsca_val3[i_species] = Qext_val3[i_species] - Qabs_val3[i_species];
        }

        const auto n_age         = free_params.n_age_;
        const auto lambda_cm_val = SEDSetting::LambdaCm();

        auto k_val2     = val2(val(N_LAMBDA), n_age);
        auto g_val2     = val2(val(N_LAMBDA), n_age);
        auto omega_val2 = val2(val(N_LAMBDA), n_age);

        auto Cabs_sum_val2 = val2(val(N_LAMBDA), n_age);

        for (auto i_age = std::size_t(0); i_age < n_age; ++i_age) {
            for (auto i_lambda = std::size_t(0); i_lambda < N_LAMBDA; ++i_lambda) {
                const auto m_total = (4. / 3 * M_PI * a_cm_val * a_cm_val * a_cm_val
                                      * (n_val3[0][i_age] * RHO_GRAIN[10] +
                                         (n_val3[1][i_age] + n_val3[2][i_age] +
                                          n_val3[3][i_age]) * RHO_GRAIN[0])).sum();

                Cabs_sum_val2[i_age][i_lambda] = (M_PI * a_cm_val * a_cm_val
                                                  * (Qabs_val3[0][i_lambda] * n_val3[0][i_age]
                                                     + Qabs_val3[1][i_lambda] * n_val3[1][i_age]
                                                     + Qabs_val3[2][i_lambda] * n_val3[2][i_age]
                                                     + Qabs_val3[3][i_lambda] *
                                                       n_val3[3][i_age])).sum();

                const auto Csca_total = (M_PI * a_cm_val * a_cm_val
                                         * (Qsca_val3[0][i_lambda] * n_val3[0][i_age]
                                            + Qsca_val3[1][i_lambda] * n_val3[1][i_age]
                                            + Qsca_val3[2][i_lambda] * n_val3[2][i_age]
                                            + Qsca_val3[3][i_lambda] * n_val3[3][i_age])).sum();

                const auto g_Csca_total = (M_PI * a_cm_val * a_cm_val
                                           * (Qsca_val3[0][i_lambda] * g_val3[0][i_lambda] *
                                              n_val3[0][i_age]
                                              + Qsca_val3[1][i_lambda] * g_val3[1][i_lambda] *
                                                n_val3[1][i_age]
                                              + Qsca_val3[2][i_lambda] * g_val3[2][i_lambda] *
                                                n_val3[2][i_age]
                                              + Qsca_val3[3][i_lambda] * g_val3[3][i_lambda] *
                                                n_val3[3][i_age])).sum();

                const auto k_abs = Cabs_sum_val2[i_age][i_lambda] / m_total;
                const auto k_sca = Csca_total / m_total;

                k_val2[i_age][i_lambda]     = k_abs + k_sca;
                omega_val2[i_age][i_lambda] = k_sca / k_val2[i_age][i_lambda];
                g_val2[i_age][i_lambda]     = g_Csca_total / Csca_total;

                if (std::isnan(k_val2[i_age][i_lambda])) k_val2[i_age][i_lambda]         = 0;
                if (std::isnan(omega_val2[i_age][i_lambda])) omega_val2[i_age][i_lambda] = 0;
                if (std::isnan(g_val2[i_age][i_lambda])) g_val2[i_age][i_lambda]         = 0;
            }
        }

        WriteAveragedPropertyFile(ofn, n_age, lambda_cm_val, k_val2, omega_val2, g_val2,
                                  Cabs_sum_val2);

        return {k_val2, omega_val2, g_val2};
    }

};

}
#endif //SED_MODEL_MAIN_H
