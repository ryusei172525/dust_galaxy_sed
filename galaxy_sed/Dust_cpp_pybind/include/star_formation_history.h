//
// Created by 西田和樹 on 2020/12/03.
//

#ifndef SED_MODEL_STAR_FORMATION_HISTORY_H
#define SED_MODEL_STAR_FORMATION_HISTORY_H

#include <constants.h>

namespace SFH {
using val = std::valarray<double>;
using val2 = std::valarray<val>;

inline double IMF(double m, int IMF_TYPE) noexcept {
    switch (IMF_TYPE) {
        case 0://Salpeter IMF (25/9/2012)
            return std::pow(m, -2.35);
        case 1://Larson IMF
            return std::pow(m, -2.35) * std::exp(-M_CH / m);
        case 2://Chabrier IMF
            if (m<= 1.0)
            return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
            return A_pow*std::pow(m, -2.3);
        case 3://Chabrier y=1
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, -2);
        case 4://Chabrier y=0.5
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, -1.5);
        case 5://Chabrier y=0.1
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, -1.1);
        case 6://Chabrier y=0
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, 0);
        case 7://Chabrier y=-0.1
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, 0.1);
        case 8://Chabrier y=-0.5
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, 0.5);
        case 9://Chabrier y=-1
            if (m<= 1.0)
                return A_exp*std::exp(-std::pow(log10(m)- log10(mc), 2.0)/(2*std::pow(sigma,2.0)) );
            else
                return A_pow*std::pow(m, 1);
        default:
            return 0;
    }
}

inline double IMFNormalization(double m_min, double m_max, int IMF_TYPE) noexcept {
    const auto dm = 0.01;
    const auto n  = std::size_t((m_max - m_min) / dm);

    auto      s1 = 0.0;
    for (auto i  = std::size_t(1); i <= n - 1; i += 2) {
        const auto m = m_min + dm * static_cast<double>(i);
        s1 += m * IMF(m, IMF_TYPE);
    }

    auto      s2 = 0.0;
    for (auto i  = std::size_t(2); i <= n - 2; i += 2) {
        const auto m = m_min + dm * static_cast<double>(i);
        s2 += m * IMF(m, IMF_TYPE);
    }
    return dm / 3.0 * (0.1 * IMF(m_min, IMF_TYPE) + 4.0 * s1 + 2.0 * s2 + 100 * IMF(m_max, IMF_TYPE));
}


inline val NormalizedIMF(int IMF_TYPE) noexcept {
    auto normed_IMF = val(N_M_STAR);

    const auto IMF_0 = IMFNormalization(0.1, 100, IMF_TYPE);

    for (auto i = std::size_t(0); i < N_M_STAR; ++i) {
        const auto m = 1.0 + 0.01 * static_cast<double>(i);
        normed_IMF[i] = IMF(m, IMF_TYPE) / IMF_0;
    }
    return normed_IMF;
}

inline double SFR(double M_gas, double schmidt_index, double tau_SF) noexcept {
    if (SFR_TYPE == 3) { // Schmitd law
        return std::pow(M_gas, schmidt_index) / tau_SF;
    } else {
        return 0;
    }
}

inline val
SFRValarray(std::size_t i_age, double schmidt_index, double tau_SF, const val& M_gas_val) noexcept {
    auto SFR_val = val(i_age + 1);

    for (auto i = std::size_t(0); i <= i_age; ++i) {
        SFR_val[i] = SFR(M_gas_val[i], schmidt_index, tau_SF);
    }
    return SFR_val;
}

// age is units in yr.
//Although Gamma function in Yoshii et al.(1996) is the second incomplete gamma function, this is wrong.
//The correct one is the first incomplete gamma function. (28/9/2012)
//Inoue (2011)
/**
 * Calculate infall mass dM/dt
 * @param age
 * @param M_total_infall
 * @param tau_infall
 * @return
 */
inline double InfallMass(double age, double M_total_infall, double tau_infall) noexcept {
    //return M_total_infall * 0.5 * (std::exp(-age / tau_infall) / tau_infall + std::exp(-age / 1e9) / 1e9); // Two component model case
    return M_total_infall * (std::exp(-age / tau_infall) / tau_infall);
}
}
#endif //SED_MODEL_STAR_FORMATION_HISTORY_H
