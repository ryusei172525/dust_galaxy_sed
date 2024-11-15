#ifndef SED_MODEL_SED_SETTING_UTIL_H
#define SED_MODEL_SED_SETTING_UTIL_H

#include <iomanip>
#include <sstream>
#include <valarray>

#include <constants.h>

namespace my_util {

using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using valt = std::valarray<std::size_t>;

namespace SED_util {

class SEDSettingUtil {
  public:
    // unit in [cm]
    /*
     * Calculate dust radius list
     *
     * @return std::valarray<double> Dust radius [cm]
     */
    static inline val DustRadiusCm() noexcept {
        auto radius_val = val(N_DUST_RADIUS);

        for (auto i = std::size_t(0); i < N_DUST_RADIUS; ++i) {
            auto ss = std::stringstream();
            ss << std::setprecision(5) << std::pow(10.0, -7.9 + static_cast<double>(i + 4) * 0.1);
            ss >> radius_val[i];
        }
        return radius_val;
    }

    /**
     * Make dust radius vector for Asano model
     * The difference between this and normal a_val is minimum radius.
     * @return val
     */
    static inline val DustRadiusDustModel() noexcept {
        auto radius_val = val(N_MAX_DUST_RADIUS);

        for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i) {
            auto ss = std::stringstream();
            ss << std::setprecision(5) << std::pow(10.0, -7.9 + static_cast<double>(i) * 0.1);
            ss >> radius_val[i];
        }
        return radius_val;
    }

    /**
     * Calculate wavelength [cm]
     * @return val
     */
    static inline val LambdaCm() noexcept {
        auto      lambda_cm_val = val(N_LAMBDA);
        for (auto i             = std::size_t(0); i < N_LAMBDA; ++i)
            lambda_cm_val[i] = 1e-7 * std::pow(10, static_cast<double>(i) / 300.0); //[um]
        return lambda_cm_val;
    }

};
}
}

#endif //SED_MODEL_SED_SETTING_UTIL_H
