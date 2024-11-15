#ifndef SPECTRA_H
#define SPECTRA_H

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cfloat>

#include <constants.h>
#include <star_formation_history.h>

namespace stellar_spectrum {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;
using valt = std::valarray<std::size_t>;

double linearhokan(double x1, double x2, double y1, double y2, double x) {
    if (x1 - x2 < DBL_EPSILON) return (y1 + y2) * 0.5;

    const double a = (y2 - y1) / (x2 - x1);

    const double y = a * (x - x1) + y1;

    if (y < std::pow(10.0, -100)) return 0;

    return y;
}

//(x1,y1) (x2,y2) から両対数グラフを線形近似して点xでのyを求める.
double InterPolate(double x1, double x2, double y1, double y2, double x) {
    if (x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0) return linearhokan(x1, x2, y1, y2, x);
    const double a = (std::log(y2) - std::log(y1)) / (std::log(x2) - std::log(x1));
    const double y = std::exp(a * (std::log(x) - std::log(x1)) + std::log(y1));
    if (y < 1e-100) return 0;
    return y;
}

/*
// もともとの SSPhokan とかいう気持ち悪いやつ
val2 SSPhokan(const val& wl_PEGASE, const val& input) {
    auto      SSPhokan      = val2(val(2701), 4);
    auto      lambda_varray = val(2701); //[um]
    auto      j             = std::size_t(0);
    for (auto i             = std::size_t(0); i < 2701; ++i) {
        lambda_varray[i] = 1e-3 * std::pow(10, static_cast<double>(i) / 300.0); //[A]
        while (wl_PEGASE[j] <= lambda_varray[i] || j == 1221 - 1) {
            ++j;
        }

        if (j >= 1221 - 1) {
            j = 1221 - 1;
            for (auto o = std::size_t(0); o < 4; ++o) {
                SSPhokan[o][i] = hokan(wl_PEGASE[j - 1], wl_PEGASE[j],
                                       input[o][j - 1], input[o][j], lambda_varray[i]);
            }
        } else if (lambda_varray[i] < 0.0912) {
            for (auto o = std::size_t(0); o < 4; ++o)SSPhokan[o][i] = 0;
        } else {
            for (auto o = std::size_t(0); o < 4; ++o) {
                SSPhokan[o][i] = hokan(wl_PEGASE[j - 1], wl_PEGASE[j],
                                       input[o][j - 1], input[o][j], lambda_varray[i]);
            }
        }
    }
    return SSPhokan;
}
*/

val InterPolateVarray(const val& x_pre_val, const val& x_new_val, const val& y_pre_val) noexcept {
    const auto n_x_new = std::size(x_new_val);
    const auto n_x_pre = std::size(x_pre_val);

    auto y_new_val = val(n_x_new);

    for (auto i = std::size_t(0); i < n_x_new; ++i) {
        auto j = static_cast<std::size_t>(std::upper_bound(std::begin(x_pre_val),
                                                           std::end(x_pre_val), x_new_val[i]) -
                                          std::begin(x_pre_val));
        if (j >= n_x_pre) {
//        if (j >= 1000) {
            j = n_x_pre - 1;
//            j = 1000;
            y_new_val[i] = InterPolate(x_pre_val[j - 1], x_pre_val[j],
                                       y_pre_val[j - 1], y_pre_val[j], x_new_val[i]);
        } else if (x_new_val[i] < 912e-8) { // Lyman break
            y_new_val[i] = 0;
        } else {
            y_new_val[i] = InterPolate(x_pre_val[j - 1], x_pre_val[j],
                                       y_pre_val[j - 1], y_pre_val[j], x_new_val[i]);
        }
    }
    return y_new_val;
}

// Search index number cooresponding to Z_SSP
std::tuple<std::size_t, double> WeightZ(std::size_t n_Z, double Z, const vec& Z_SSP_vec) noexcept {
    if (n_Z <= 1) return std::make_tuple(0, 1.0);
    if (Z <= Z_SSP_vec[0]) return std::make_tuple(0, 1.0);
    if (Z >= Z_SSP_vec.back()) return std::make_tuple(n_Z - 2, 0.0);

    // i は Z と Z_SSP_varray の大小関係が入れ替わるところのインデックス
    const auto i_Z = static_cast<size_t>(std::lower_bound(std::begin(Z_SSP_vec),
                                                          std::end(Z_SSP_vec), Z)
                                         - std::begin(Z_SSP_vec)) - 1;
    // Debug
    if (i_Z + 1 == n_Z) std::cout << "Erorr: i_Z_sup = n_Z_" << std::endl;
    // alpha は小さい方から何割大きい方に近づいているかを表す
    if (Z_SSP_vec[i_Z + 1] > 0 && Z > 0) {
        const auto alpha = std::log(Z_SSP_vec[i_Z + 1] / Z)
                           / std::log(Z_SSP_vec[i_Z + 1] / Z_SSP_vec[i_Z]);
        return std::make_tuple(i_Z, alpha);
    } else { // 負の値が出てくるときは線形で補間だがでないはず
        // Debug
        std::cout << "Z_SSP is minus" << std::endl;
        const auto alpha = (Z_SSP_vec[i_Z + 1] - Z) / (Z_SSP_vec[i_Z + 1] - Z_SSP_vec[i_Z]);
        return std::make_tuple(i_Z, alpha);
    }
}

double Steffen(std::size_t n, const val& x, const val& y, bool condition, double t, std::size_t i) {
    //return z?? i =0~n-2?? condition=FALSEorTRUE
    double a, b, d;
    double h_i, h_im1, h_ip1;
    double s_i, s_im1, s_ip1;
    double p_i, p_ip1, y1_i, y1_ip1;
    double tmx_i;
    double z;

    if (i == 0) {
        h_i    = x[i + 1] - x[i];
        s_i    = (y[i + 1] - y[i]) / h_i;
        h_ip1  = x[i + 2] - x[i + 1];
        s_ip1  = (y[i + 2] - y[i + 1]) / h_ip1;
        p_ip1  = (s_i * h_ip1 + s_ip1 * h_i) / (h_i + h_ip1);
        y1_ip1 = (copysign(1.0, s_i) + copysign(1.0, s_ip1))
                 * std::fmin(std::fmin(std::fabs(s_i), std::fabs(s_ip1)), 0.5 * std::fabs(p_ip1));

        if (condition) {
            y1_i = (6.0 * s_i * x[0] * (x[0] + h_i)
                    - y1_ip1 * x[0] * (3.0 * x[0] + 2. * h_i))
                   / (3.0 * x[0] * x[0] + 4.0 * x[0] * h_i + h_i * h_i);
        } else {
            y1_i = 1.5 * s_i - 0.5 * y1_ip1;
        }
    } else if (i == n - 2) {
        h_i    = x[i + 1] - x[i];
        s_i    = (y[i + 1] - y[i]) / h_i;
        h_im1  = x[i] - x[i - 1];
        s_im1  = (y[i] - y[i - 1]) / h_im1;
        p_i    = (s_im1 * h_i + s_i * h_im1) / (h_im1 + h_i);
        y1_i   = (copysign(1.0, s_im1) + copysign(1.0, s_i))
                 * std::fmin(std::fmin(std::fabs(s_im1), std::fabs(s_i)), 0.5 * std::fabs(p_i));
        y1_ip1 = 1.5 * s_i - 0.5 * y1_i;
    } else {
        h_i    = x[i + 1] - x[i];
        s_i    = (y[i + 1] - y[i]) / h_i;
        h_im1  = x[i] - x[i - 1];
        s_im1  = (y[i] - y[i - 1]) / h_im1;
        h_ip1  = x[i + 2] - x[i + 1];
        s_ip1  = (y[i + 2] - y[i + 1]) / h_ip1;
        p_i    = (s_im1 * h_i + s_i * h_im1) / (h_im1 + h_i);
        p_ip1  = (s_i * h_ip1 + s_ip1 * h_i) / (h_i + h_ip1);
        y1_i   = (copysign(1.0, s_im1) + copysign(1.0, s_i))
                 * std::fmin(std::fmin(std::fabs(s_im1), std::fabs(s_i)), 0.5 * std::fabs(p_i));
        y1_ip1 = (copysign(1.0, s_i) + copysign(1.0, s_ip1))
                 * std::fmin(std::fmin(std::fabs(s_i), std::fabs(s_ip1)), 0.5 * std::fabs(p_ip1));
    }

    if (t > x[n - 1]) {
        z = y[n - 1] + y1_ip1 * (t - x[n - 1]);
    } else {
        tmx_i = t - x[i];
        a     = (y1_i + y1_ip1 - 2.0 * s_i) / (h_i * h_i);
        b     = (3.0 * s_i - 2.0 * y1_i - y1_ip1) / h_i;
        d     = y[i];
        z     = ((a * tmx_i + b) * tmx_i + cl * 1e8) * tmx_i + d;
    }

    return z;
}

std::size_t Bracket(std::size_t n, const val& x, double t, std::size_t i) {
    //return i??

    auto   di = std::size_t(0);
    int    niter;
    double xinf, xsup;

    auto iinf = std::size_t(0);

    if (i == 0) {
        niter = 0;
        iinf  = std::size_t(0);
        auto isup = n - 1;
        if (t > x[iinf]) {
            if (t >= x[isup]) {
                iinf = isup - 1;
            } else {
                while (iinf + 1 < isup) {
                    niter           = niter + 1;
                    const auto imed = (iinf + isup) / 2;
                    if (t <= x[imed]) {

                        isup = imed;
                    } else {
                        iinf = imed;
                    }
                }
            }
        }
    } else {
        niter = 0;
        iinf  = std::size_t(0);
        auto isup = n - 1;
        if (t > x[iinf]) {
            if (t >= x[isup]) {
                iinf = isup - 1;
            } else {
                di = 1;
                if (t >= x[i]) {
                    iinf = i;
                    isup = static_cast<std::size_t>(std::fmin(n, i + di));
                    xsup = x[isup];
                    while (t > xsup) {
                        niter = niter + 1;
                        di    = 2 * di;
                        isup  = static_cast<std::size_t>(std::fmin(n, i + di));
                        xsup  = x[isup];
                    }
                } else {
                    isup = i;
                    iinf = static_cast<std::size_t>(std::fmax(1, i - di));
                    xinf = x[iinf];
                    while (t < xinf) {
                        niter = niter + 1;
                        di    = 2 * di;
                        iinf  = static_cast<std::size_t>(std::fmax(1, i - di));
                        xinf  = x[iinf];
                    }
                }
                while (iinf + 1 < isup) {
                    niter           = niter + 1;
                    const auto imed = (iinf + isup) / 2;
                    if (t <= x[imed]) {
                        isup = imed;
                    } else {
                        iinf = imed;
                    }
                }

            }
        }
    }
    i = iinf;
    return i;
}


// t is Myr
double
StarFormationRate(double t, double SFR_param1, double SFR_param2, double sigma_gas) noexcept {
    auto SFR = 0.0;

    //instantaneous burst
    if (SFR_TYPE == 0) {
        return 0;
    }

    //constant SFR
    if (SFR_TYPE == 1) {
        if (t <= SFR_param2) return SFR_param1;
        else return 0;
    }

    // Exponentially decreasing or increasing SFR
    if (SFR_TYPE == 2) {
        SFR = SFR_param2 * std::exp(-t / SFR_param1) / SFR_param1;
    }

    // SFR proportional to a power of the mass of gas
    if (SFR_TYPE == 3) {
        SFR = std::pow(sigma_gas, SFR_param1) / SFR_param2;
    }

    if (t >= T_WIND) {
        SFR = 0;
    }

    //maximal SFR

    if (SFR >= sigma_gas) {
        SFR = sigma_gas;
        std::cout << "Error: SFR is Excess!" << std::endl;
    }
    return SFR;
}

void DustComposition(const val& Zext, const val& frac, double Z, val& coeffextZ) {
    auto i = std::size_t(0);
    i = Bracket(5, Zext, Z, i);
    coeffextZ[0] = Steffen(5, Zext, frac, false, Z, i);
    coeffextZ[1] = 1.0 - coeffextZ[0];
}

double NebularFraction(std::size_t i0, std::size_t n_spitzer, double Z_gas,
                       const val& tau_dust_Spitzer_val, const val& y_Spitzer_val) noexcept {
    if (isNEBULAR_EMISSION) {
        const auto code_ext = 0;
        if (code_ext == false) return 1.0;
        const double tau_dust = 0.5 * Z_gas / Zsun;
        i0 = Bracket(n_spitzer, tau_dust_Spitzer_val, tau_dust, i0);
        const double y = Steffen(n_spitzer, tau_dust_Spitzer_val, y_Spitzer_val,
                                 false, tau_dust, i0);
        return y * y * y;
    }
    return 0.0;
}

void
MassEvolution(std::size_t n_Z, double max_age, double schmidt_index, bool is_infall,
              double tau_infall, double tau_SF, double sigma_Z, const vec& Z_SSP_vec,
              const val2& ejecta_val2, const val2& ejecta_Z_val2, const val2& m_BHNS_val2,
              const val2& m_WD_val2, const val2& m_alive_val2, valt& i_Z_val, val& SFR_val,
              val& sigma_gas_val, val& Z_SFR_val, val& alpha_val, val& Z_gas_val,
              val& SFR_lum_val, val& age_star_val, val& sigma_star_val, val& Z_star_val,
              val& sigma_sub_val, val& m_gal_val, val& sigma_BHNS_val, val& sigma_WD_val) noexcept {
    for (auto age = std::size_t(0); age < static_cast<std::size_t>(max_age / 1e6) + 1; ++age) {
        const auto tua_SF_Myr = tau_SF / 1e6; // convert from yr to Myr

        SFR_val[age] = SFH::SFR(sigma_gas_val[age], schmidt_index, tua_SF_Myr); // Msun/Myr

        const auto[i_Z_age, alpha_age] = WeightZ(n_Z, Z_SFR_val[age], Z_SSP_vec);
        i_Z_val[age]   = i_Z_age;
        alpha_val[age] = alpha_age;

        SFR_lum_val[age] = (1.0 - F_SUB_STELLAR) * SFR_val[age];

        auto d_ejecta     = 0.0;
        auto d_ejectaZ    = 0.0;
        auto d_sigma_BHNS = 0.0;
        auto d_sigma_WD   = 0.0;

        auto m_Z_star = 0.0; // 星の中の金属の質量

        // 銀河の始まりから今までの積分
        for (auto i = std::size_t(0); i <= age; ++i) {
            const auto i_Z      = i_Z_val[i];
            const auto age_star = age - i; // age_star represents the age of star
            const auto alpha    = alpha_val[i]; // corresponding to Z
            const auto SFR      = SFR_lum_val[i];

            const auto ejecta_i = alpha * SFR * ejecta_val2[i_Z][age_star] +
                                  (1 - alpha) * SFR * ejecta_val2[i_Z + 1][age_star];
            d_ejecta += ejecta_i;
            d_ejectaZ += alpha * SFR * ejecta_Z_val2[i_Z][age_star]
                         + (1 - alpha) * SFR * ejecta_Z_val2[i_Z + 1][age_star]
                         + Z_gas_val[i] * ejecta_i;
            d_sigma_BHNS += alpha * SFR * m_BHNS_val2[i_Z][age_star]
                            + (1 - alpha) * SFR * m_BHNS_val2[i_Z + 1][age_star];
            d_sigma_WD += alpha * SFR * m_WD_val2[i_Z][age_star]
                          + (1 - alpha) * SFR * m_WD_val2[i_Z + 1][age_star];

            const auto sigma_star_i = alpha * SFR * m_alive_val2[i_Z][age_star] +
                                      (1 - alpha) * SFR * m_alive_val2[i_Z + 1][age_star];
            sigma_star_val[age] += sigma_star_i;
            m_Z_star += Z_SFR_val[i] * sigma_star_i;
            age_star_val[age] += static_cast<double>(age_star) * sigma_star_i;
        }

        auto d_m_infall = 0.0;
        if (is_infall)
            d_m_infall =
            //SFH::InfallMass(static_cast<double>(age) * 1e6, M_total_infall, tau_infall) * 1e6;
            SFH::InfallMass(static_cast<double>(age) * 1e6, 1.0, tau_infall) * 1e6;

        const auto d_sigma_gas = d_ejecta - SFR_val[age] + d_m_infall;
        const auto d_sigma_Z   = d_ejectaZ - SFR_val[age] * Z_gas_val[age] + d_m_infall * Z_INFALL;

        if (sigma_star_val[age] > 0) {
            age_star_val[age] /= sigma_star_val[age];
            Z_star_val[age] = m_Z_star / sigma_star_val[age];
        } else {
            age_star_val[age]   = 0;
            sigma_star_val[age] = 0;
            Z_star_val[age]     = m_Z_star;
        }

        sigma_sub_val[age + 1] = sigma_sub_val[age] + F_SUB_STELLAR * SFR_val[age];

        if (static_cast<double>(age) >= T_WIND) {
            m_gal_val[age + 1]      = m_gal_val[age] - sigma_gas_val[age] - d_ejecta;
            sigma_BHNS_val[age + 1] = sigma_BHNS_val[age] + d_sigma_BHNS;
            sigma_WD_val[age + 1]   = sigma_WD_val[age] + d_sigma_WD;
            sigma_gas_val[age + 1]  = 0;
            sigma_Z = 0;
            Z_gas_val[age + 1] = 0;
        } else {
            m_gal_val[age + 1]      = m_gal_val[age] + d_m_infall;
            //std::cout << "m_gal[" << age << "] = " << m_gal_val[age] <<
            //        ", m_star = " << sigma_star_val[age] << ", sigma_gas = "
            //      << sigma_gas_val[age] << std::endl;
            sigma_BHNS_val[age + 1] = sigma_BHNS_val[age] + d_sigma_BHNS;
            sigma_WD_val[age + 1]   = sigma_WD_val[age] + d_sigma_WD;
            sigma_gas_val[age + 1]  = std::max(sigma_gas_val[age] + d_sigma_gas, 0.0);
            sigma_Z = std::max(sigma_Z + d_sigma_Z, 0.0);
            if (sigma_gas_val[age + 1] <= 0) {
                Z_gas_val[age + 1] = 0;
            } else
                Z_gas_val[age + 1] = sigma_Z / sigma_gas_val[age + 1];
        }
        const auto is_code_Z = true;
        if (is_code_Z == 1) {
            Z_SFR_val[age + 1] = Z_gas_val[age + 1];
        } else {
            Z_SFR_val[age + 1] = Z_SFR_val[0];
        }
    }
}

val2 ContinuumFluxFromBornAtTStar(std::size_t n_lambda, std::size_t age, double f,
                                  const valt& i_Z_val, const valt& t_inv_val, const val& alpha_val,
                                  const val& beta_val, const val& F_neb_val, const val& SFR_lum_val,
                                  const val2& n_Lym_val2, const val3& F_SSP_val3) noexcept {
    auto F_gal_val2 = val2(val(age + 1), n_lambda);

    //for (auto li = std::size_t(0); li < n_lambda; ++li) {
    for (auto li = std::size_t(121); li < n_lambda; ++li) {
        auto F_f_neb = f * F_neb_val[li];

        for (auto t = std::size_t(0); t <= age; ++t) {
            const auto i_Z   = i_Z_val[t];
            const auto q     = t_inv_val[age - t]; // 星の年齢
            const auto beta  = beta_val[age - t]; // age
            const auto alpha = alpha_val[t]; // metallicity
            F_gal_val2[li][t] = SFR_lum_val[t] * (
            alpha * beta * (F_SSP_val3[i_Z][q][li] + n_Lym_val2[i_Z][q] * F_f_neb)
            + (beta - alpha * beta) *
              (F_SSP_val3[i_Z + 1][q][li] + n_Lym_val2[i_Z + 1][q] * F_f_neb)
            + (alpha - beta * alpha) *
              (F_SSP_val3[i_Z][q + 1][li] + n_Lym_val2[i_Z][q + 1] * F_f_neb)
            + (1 - alpha - beta + alpha * beta) *
              (F_SSP_val3[i_Z + 1][q + 1][li] + n_Lym_val2[i_Z + 1][q + 1] * F_f_neb));
        }
    }
    return F_gal_val2;
}

val ContinuumFlux(std::size_t n_lambda, const val2& F_gal_val2) {
    auto      F_gal_val = val(n_lambda);
    for (auto i         = std::size_t(0); i < n_lambda; ++i) {
        F_gal_val[i] = F_gal_val2[i].sum();
    }
    return F_gal_val;
}

void
BolometricFlux(std::size_t age, const valt& i_Z_val, const val& alpha_val, const val& beta_val,
               const valt& t_inv_val, const val& SFR_lum_val, const val2& F_bol_SSP_val2,
               const val& Z_SFR_val) {
    double F_bol     = 0.0;
    double F_bol_Z   = 0.0;
    double F_bol_age = 0.0;

    for (auto i = std::size_t(0); i <= age; ++i) {
        const auto i_Z   = i_Z_val[i];
        const auto alpha = alpha_val[i];
        const auto q     = age - i;
        const auto beta  = beta_val[q];
        const auto t_inv = t_inv_val[q];
        F_bol += SFR_lum_val[i] * (beta * (alpha * F_bol_SSP_val2[i_Z][t_inv] +
                                           (1.0 - alpha) * F_bol_SSP_val2[i_Z + 1][t_inv])
                                   + (1.0 - beta) * (alpha * F_bol_SSP_val2[i_Z][t_inv + 1] +
                                                     (1.0 - alpha) *
                                                     F_bol_SSP_val2[i_Z + 1][t_inv + 1]));
        F_bol_Z += SFR_lum_val[i] * Z_SFR_val[i]
                   * (beta * (alpha * F_bol_SSP_val2[i_Z][t_inv] +
                              (1.0 - alpha) * F_bol_SSP_val2[i_Z + 1][t_inv])
                      + (1.0 - beta) * (alpha * F_bol_SSP_val2[i_Z][t_inv + 1] +
                                        (1.0 - alpha) * F_bol_SSP_val2[i_Z + 1][t_inv + 1]));
        F_bol_age += SFR_lum_val[i] * static_cast<double>(q)
                     * (beta * (alpha * F_bol_SSP_val2[i_Z][t_inv] +
                                (1.0 - alpha) * F_bol_SSP_val2[i_Z + 1][t_inv])
                        + (1.0 - beta) * (alpha * F_bol_SSP_val2[i_Z][t_inv + 1] +
                                          (1.0 - alpha) * F_bol_SSP_val2[i_Z + 1][t_inv + 1]));
    }
}

val LineFlux(std::size_t age, const val& F_line_val, const val& SFR_lum_val, const val& alpha_val,
             const val& beta_val, const val2& n_Lym_val2, double f, const valt& i_Z_val,
             const valt& t_inv_val) noexcept {
    const auto n_line = F_line_val.size();

    auto F_line_tot_val = val(n_line);

    for (auto i = std::size_t(0); i < n_line; ++i) {
        F_line_tot_val[i] = 0;
        for (auto k = std::size_t(std::fmax(age - 50, 0)); k <= age; ++k) {
            const auto q = age - k;
            // alpha: Z に関連したビン補正係数
            // beta: t に関連したビンの補正係数
            F_line_tot_val[i] += SFR_lum_val[k]
                                 * (beta_val[q]
                                    * (alpha_val[k] * f * n_Lym_val2[i_Z_val[k]][t_inv_val[q]] *
                                       F_line_val[i]
                                       + (1.0 - alpha_val[k]) * f *
                                         n_Lym_val2[i_Z_val[k] + 1][t_inv_val[q]] * F_line_val[i])
                                    + (1.0 - beta_val[q])
                                      * (
                                      alpha_val[k] * f * n_Lym_val2[i_Z_val[k]][t_inv_val[q] + 1] *
                                      F_line_val[i]
                                      + (1.0 - alpha_val[k]) * f *
                                        n_Lym_val2[i_Z_val[k] + 1][t_inv_val[q] + 1] *
                                        F_line_val[i]));
        }
    }
    return F_line_tot_val;
}

val2 YoungStellarFraction(std::size_t n_lambda, const val& lambda_angstrom_val, std::size_t age_gal,
                          const val2& F_gal_val2, const val& F_gal_val) noexcept {
    auto f_young_star_val2 = val2(val(n_lambda), 3);

    for (auto li = std::size_t(0); li < n_lambda; ++li) {
        if (lambda_angstrom_val[li] >= 900) {
            auto F_gal_numerator = 0.0;

            for (auto age_star = std::size_t(0); age_star <= age_gal; ++age_star) {
                const auto t = age_gal - age_star;
                F_gal_numerator += F_gal_val2[li][t];
                if (age_star == 10) { // young star threshold = 10 Myr
                    f_young_star_val2[0][li] = F_gal_numerator / F_gal_val[li];
                }
                if (age_star == 100) { // young star threshold = 100 Myr
                    f_young_star_val2[1][li] = F_gal_numerator / F_gal_val[li];
                }
                if (age_star == 1000) { // young star threshold = 1 Gyr
                    f_young_star_val2[2][li] = F_gal_numerator / F_gal_val[li];
                }
            }

            if (age_gal % 10 == 0 && age_gal != 0) {
                if (age_gal <= 10) {
                    f_young_star_val2[0][li] = 1.0;
                }
                if (age_gal <= 100) {
                    f_young_star_val2[1][li] = 1.0;
                }
                if (age_gal <= 1000) {
                    f_young_star_val2[2][li] = 1.0;
                }
            }
        } else {
            f_young_star_val2[0][li] = f_young_star_val2[1][li] = f_young_star_val2[2][li] = 0;
        }
    }
    return f_young_star_val2;
}

double NumberOfLymanContinuumPhoton(std::size_t age, const valt& i_Z_val, const val& SFR_lum_val,
                                    const val& alpha_val, const val2& n_Lym_val2) {
    //Calculation of the number of Lyman continuum photons.
    auto      n_Ly_tot = 0.0;
    for (auto i        = std::size_t(std::fmax(0, age - 50)); i <= age; ++i) {
        n_Ly_tot += SFR_lum_val[i]
                    * (alpha_val[i] * n_Lym_val2[i_Z_val[i]][age - i]
                       + (1.0 - alpha_val[i]) * n_Lym_val2[i_Z_val[i] + 1][age - i]);
    }
    return n_Ly_tot;
}

double NumberOfSNII(std::size_t age, const valt& i_Z_val, const val& SFR_lum_val,
                    const val& alpha_val, const val2& n_SNII_val2) {
    //Calculation of the number of SNII
    auto      n_SNII_tot = 0.0;
    for (auto i          = std::size_t(std::fmax(0, age - 100)); i <= age; ++i) {
        n_SNII_tot += SFR_lum_val[i]
                      * (alpha_val[i] * n_SNII_val2[i_Z_val[i]][age - i]
                         + (1.0 - alpha_val[i]) * n_SNII_val2[i_Z_val[i] + 1][age - i]);
    }
    return n_SNII_tot;
}

double NumberOfSNIa(std::size_t age, const valt& i_Z_val, const val& SFR_lum_val,
                    const val& alpha_val, const val2& n_SNIa_val2) {
    // Calculation of the number of SNIa.
    auto      n_SNIa_tot = 0.0;
    for (auto i          = std::size_t(0); i <= age; ++i) {
        n_SNIa_tot += SFR_lum_val[i]
                      * (alpha_val[i] * n_SNIa_val2[i_Z_val[i]][age - i]
                         + (1.0 - alpha_val[i]) * n_SNIa_val2[i_Z_val[i] + 1][age - i]);
    }
    return n_SNIa_tot;
}

[[maybe_unused]] void
Extinction([[maybe_unused]] std::size_t age, [[maybe_unused]] std::size_t i912, double F_bol,
           const val& Z_ext_val, const val& f_gra_val, const val& Z_gas_val,
           const val& lambda_angstrom_val,
           val& F_gal_val, val& F_line_tot_val) {
    //Extinction
    auto       F_ext       = 0.0;
    auto       coeff_ext_Z = val(2);
    const auto code_ext    = 0;
    if (code_ext != 0) {
        DustComposition(Z_ext_val, f_gra_val, Z_gas_val[age], coeff_ext_Z);

        // ***********If there is some extinction, the Lyman continuum photons
        // ***********not absorbed by the gas are absorbed by the dust.
        if (Z_gas_val[age] > 0.0) {
            for (auto i = std::size_t(0); i < i912; ++i) {
                F_ext += (F_gal_val[i] + F_gal_val[i + 1])
                         * (lambda_angstrom_val[i + 1] - lambda_angstrom_val[i]) * 0.5;
                //F_gal_val[i] = 0.0;
            }
            F_ext += (912.0 - lambda_angstrom_val[i912]) * F_gal_val[i912];
            F_gal_val[i912] = 0.0;
            // ***********Idem for Lyman alpha.
            F_ext += F_line_tot_val[30];
            F_line_tot_val[30] = 0.0;
        }
        // ***********Ellipticals.
    }

    if (F_bol > 0.0) {
        F_ext /= F_bol;
    } else {
        F_ext = 0.0;
    }
}

void
WriteContinuumFile(const val& lambda_cm_val, const val& F_gal_val, const val& lambda_cm_output_val,
                   FILE* ofp) {
    const auto F_gal_result_val = InterPolateVarray(lambda_cm_val, lambda_cm_output_val, F_gal_val);

    for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
        fprintf(ofp, "%5e %5e\n", lambda_cm_output_val[li], F_gal_result_val[li] * 1e8 * L_sun);
    }
}

void WriteYoungStellarFractionFile(const val& lambda_cm_output_val, const val2& f_young_star_val2,
                                   const val&& lambda_cm_val, std::ofstream& ofs) {
    auto f_young_star_output_val2 = val2(val(N_LAMBDA), 3);

    for (auto i = std::size_t(0); i < 3; ++i) {
        f_young_star_output_val2[i] = InterPolateVarray(lambda_cm_val, lambda_cm_output_val,
                                                        f_young_star_val2[i]);
        // When lambda_um_output > lambda_um, f_young may be larger than 1.
        std::for_each(std::begin(f_young_star_output_val2[i]),
                      std::end(f_young_star_output_val2[i]),
                      [](auto&& f_young) { f_young = std::clamp(f_young, 0.0, 1.0); });
        f_young_star_output_val2[i] = MoveAverage(f_young_star_output_val2[i]);
    }

    for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
        ofs << lambda_cm_output_val[li] << " " << f_young_star_output_val2[0][li] << " "
            << f_young_star_output_val2[1][li] << " " << f_young_star_output_val2[2][li]
            << std::endl;
    }
}
}
#endif
