//
// Created by 西田和樹 on 2020/08/05.
//

#ifndef SED_MODEL_SED_FILE_UTIL_H
#define SED_MODEL_SED_FILE_UTIL_H

#include <valarray>

#include <constants.h>
#include <file_name.h>
#include <file_util.h>

namespace my_util::SED_util {
class SEDFileUtil : public my_util::FileUtil {
    using val = std::valarray<double>;
    using val2 = std::valarray<val>;
    using val3 = std::valarray<val2>;

    using valt = std::valarray<std::size_t>;

  public:
    /**
     * Read only dust to gas mass ratio from total dust mass file calculated by Asano model
     * @param free_params Free parameter class object
     * @param ifn Input file name
     * @return val Dust to gas mass ratio each galaxy age.
     */
    static val ReadDustToGasMassRatio(const std::string& ifn, std::size_t n_age) {
        auto ifs = IfsOpen(ifn);
        SkipLine(ifs); // skip header
        auto dust_to_gas_mass_ratio_val = val(n_age);
        for (auto&& dust_to_gas_mass_ratio : dust_to_gas_mass_ratio_val) {
            for (auto i = std::size_t(0); i < 12; ++i) ReadDouble(ifs); // Skip unused values
            dust_to_gas_mass_ratio = ReadDouble(ifs);
            SkipLine(ifs);
        }
        return dust_to_gas_mass_ratio_val;
    }

    // TODO: いちいちファイルから読み込むのは変なので、Infall モデルから計算する形にした方が良いと思う。
    /**
    * Read only total galaxy mass from total dust mass file calculated by Asano model
    * @param free_params Free parameter class object
    * @param ifn Input file name
    * @return val Total galaxy mass each galaxy age unit in Msun.
    */
    static val ReadTotalGalaxyMassRatio(const std::string& ifn, std::size_t n_age) {
        auto ifs = IfsOpen(ifn);
        SkipLine(ifs); // skip header
        auto M_total_gal_val = val(n_age);
        for (auto&& M_total_gal : M_total_gal_val) {
            ReadDouble(ifs); // Skip age
            M_total_gal = ReadDouble(ifs);
            SkipLine(ifs);
        }
        return M_total_gal_val;
    }

    /**
     * Read dust number/mass distribution file calculated by Asano model
     * This number/mass distribution is normalized by total galaxy mass
     * @param ifn Input file name
     * @param n_age Number of age
     * @return val3 [carbon or silicate][galaxy age][dust radius]
     * 0: silicate, 1: carbon
     */
    inline static val3
    ReadSilicateAndCarbonDustDistributionFile(const std::string& ifn, std::size_t n_age) {
        auto ifs = IfsOpen(ifn);
        SkipLine(ifs);

        auto      C_val2   = val2(val(N_DUST_RADIUS), n_age);
        auto      sil_val2 = val2(val(N_DUST_RADIUS), n_age);
        for (auto i_age    = std::size_t(0); i_age < n_age; ++i_age) {
            for (auto i_a = std::size_t(0); i_a < N_DUST_RADIUS; ++i_a) {
                ReadDouble(ifs); // skip dust radius
                ifs >> C_val2[i_age][i_a] >> sil_val2[i_age][i_a];
                SkipLine(ifs);
            }
        }
        return {sil_val2, C_val2};
    }

    static inline val2
    ReadDustPropertyFile(const std::string& dust_species, const std::string& parameter) {
        auto ifs = IfsOpen(DustPropertyFileName(dust_species, parameter));

        SkipLine(ifs); // skip header

        auto n_a = N_DUST_RADIUS;
        if (dust_species == "PAHneu" || dust_species == "PAHion") n_a = N_PAH_RADIUS;

        // In PAH case, large radius regions are filled by 0.
        auto param_val2 = val2(val(N_DUST_RADIUS), N_LAMBDA);

        for (auto i = std::size_t(0); i < N_LAMBDA; ++i) {
            ReadDouble(ifs); // Skip wavelength
            for (auto j = std::size_t(0); j < n_a; ++j) {
                ifs >> param_val2[i][j];
            }
        }
        return param_val2;
    }

    /**
     * Read averaged dust property file
     * @param ifn Input file name
     * @param k_val2 Extinction coefficient pur unit mass
     * @param omega_val2 Scattering albedo
     * @param g_val2 Asymmetry parameter
     */
    inline static void
    ReadDustAveragedPropertyFile(const std::string& ifn, val2& k_val2, val2& omega_val2,
                                 val2& g_val2) {
        auto ifs = IfsOpen(ifn);
        SkipLine(ifs); // Skip header

        for (auto ai = std::size_t(0); ai < k_val2.size(); ++ai) {
            for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
                auto lambda = 0.0;
                ifs >> lambda >> k_val2[ai][li] // extinction coefficient per unit mass [g^-1]
                    >> omega_val2[ai][li] // omega: scattering albedo
                    >> g_val2[ai][li]; // g: asymmetry parameter
            }
        }
    }

    /**
     * Read sum of absorption coefficient of dust from file
     * @param ifn Input file name
     * @param n_age Number of age
     * @return val2 [galaxy age][wavelength]
     */
    inline static val2 ReadSumAbsorptionCoefficient(const std::string& ifn, std::size_t n_age) {
        auto ifs = IfsOpen(ifn);
        SkipLine(ifs); // Skip header

        auto Cabs_sum_val2 = val2(val(N_LAMBDA), n_age);

        for (auto i_age = std::size_t(0); i_age < n_age; ++i_age) {
            for (auto i_lambda = std::size_t(0); i_lambda < N_LAMBDA; ++i_lambda) {
                auto buff = 0.0;
                ifs >> buff >> buff >> buff >> buff >> Cabs_sum_val2[i_age][i_lambda];
            }
        }
        return Cabs_sum_val2;
    }

    [[nodiscard]] inline static val ReadStellarContinuumFile(const std::string& ifn) {
        auto ifs = IfsOpen(ifn);
        SkipLine(ifs); // skip header
        auto L_pegase_val = val(N_LAMBDA);
        for (auto&& L_pegase : L_pegase_val) {
            ReadDouble(ifs); // Skip wavelength
            ifs >> L_pegase;
        }
        return L_pegase_val;
    }

    [[nodiscard]] inline static val ReadYoungStellarFraction(const std::string& ifn) {
        auto       ifs         = IfsOpen(ifn);
        const auto head        = SkipLine(ifs); // skip header
        auto       f_young_val = val(N_LAMBDA);

        for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
            auto d = 0.0;
            ifs >> d >> f_young_val[li] >> d >> d;
        }
        return f_young_val;
    }

    [[nodiscard]] inline static val
    CalcEnergyDensity(const val& L_pegase_val, const val& Cabs_sum_val, double M_gal,
                      const val& transmission_rate_val) noexcept {

        auto      u_cgs = val(N_LAMBDA);
        for (auto i     = std::size_t(0); i < N_LAMBDA; ++i) {
            u_cgs[i] =
            (1 - transmission_rate_val[i]) * L_pegase_val[i] / (cl * M_gal * Cabs_sum_val[i]);
            if (u_cgs[i] < DBL_MIN) u_cgs[i] = 0;
        }

        return u_cgs;
    }

    [[nodiscard]] static inline val2 ReadQabsFile(std::size_t n_a, const std::string& fn) {
        auto ifs       = IfsOpen(fn);
        auto Qabs_val2 = val2(val(n_a), N_LAMBDA);

        SkipLine(ifs); // skip header

        for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
            auto wavelength = 0.0;
            ifs >> wavelength;
            for (auto di = std::size_t(0); di < n_a; ++di) {
                ifs >> Qabs_val2[li][di];
            }
        }
        return Qabs_val2;
    }

    [[nodiscard]] static inline val ReadTransmissionRateFile(const FreeParameter& free_params) {
        auto ifs = my_util::FileUtil::IfsOpen(TransmissionRateFileName(free_params));
        SkipLine(ifs);

        auto transmission_rate_val = val(N_LAMBDA);
        for (auto&& transmission_rate : transmission_rate_val) {
            auto buff = 0.0;
            ifs >> buff >> transmission_rate >> buff;
        }
        return transmission_rate_val;
    }

};
}
#endif //SED_MODEL_SED_FILE_UTIL_H
