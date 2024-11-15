#ifndef MKSED_RADIATIVE_TRANSFER_H
#define MKSED_RADIATIVE_TRANSFER_H

#include <utility>
#include <valarray>

#include <cfloat>

#include <calculation_util.h>
#include <constants.h>
#include <SED_file_util.h>
#include <SED_setting_util.h>

namespace radiative_transfer {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using SEDSetting = my_util::SED_util::SEDSettingUtil;
using SEDFile = my_util::SED_util::SEDFileUtil;

class RadiativeTransfer {
  private:
    const val lambda_cm_val_ = SEDSetting::LambdaCm(); // wavelength in cm

    double D_; // dust mass / gas mass
    val    f_young_val_; // 銀河の中の星のある基準よりも若い星の割合 stellar_spectrum で計算される
    val    kappa_eff_val_; // effective extinction coefficient per unit length
    val    g_eff_val_; // effective asymmetry parameter
    val    omega_eff_val_; // effective albedo
    val    Pesc_val_;
    val    tau_val_; // effective optical depth
    val    mu_val_; // mu = cos(theta)

    /**
     * Probability emitted externally emitted isotropic photons will interact with
     * and get absorbed by a homogeneous medium contained in a sphere.
     * @param tau Optical depth
     * @return val Probability
     */
    [[nodiscard]] static inline val PintVal(const val& tau) noexcept {
        const auto n_tau = std::size(tau);
        auto       p_int = val(n_tau);
        for (auto  i     = std::size_t(0); i < n_tau; ++i) {
            // TODO: この式の正当性を確かめる
            //  tau が非常に小さい時は計算が違う
            if (tau[i] < 1e-2) //p_int[i] = 4 * tau[i] / 3;
                p_int[i] = tau[i] / 0.75 - tau[i] * tau[i];
            else {
                const auto A = 1 / (2 * tau[i] * tau[i]);
                p_int[i] = 1 - A + (1 / tau[i] + A) * std::exp(-2 * tau[i]);
            }
        }
        return p_int;
    };

    // The scattering escape probability that internally emitted isotropic photons
    [[nodiscard]] static inline val ProbabilityEscape(const val& tau_val, const val& omega_val,
                                                      const val& Pint_val) noexcept {
        const auto Pe_val = 0.75 / tau_val * Pint_val;
        return Pe_val / (1 - omega_val * (1 - Pe_val));
    }


    [[nodiscard]] static inline val
    GClumpVal(const val& tau_clump, const val& g_val, const val& omega) noexcept {
        auto g_clump_val = val(N_LAMBDA);

        for (auto i = std::size_t(0); i < N_LAMBDA; ++i) {
            auto g = g_val[i];
            if (g < 0.2) {
                g_clump_val[i] = g;
                continue;
            } else if (g < 0) g = 0.0;

            const auto A = 1.5 + 4 * g * g * g + 2 * omega[i] * std::sqrt(g) * std::exp(-5 * g);
            const auto B = 2 - g * (1 - g) - 2 * omega[i] * g;
            const auto C = 1 / (3 - std::sqrt(2 * g) - 2 * omega[i] * g * (1 - g));

            g_clump_val[i] = g - C * (1 - (1 + std::exp(-B / A))
                                          / (1 + std::exp((tau_clump[i] - B) / A)));
            if (std::isnan(g_clump_val[i]))
                std::cout << "i = " << i << ", g[i] = " << g_val[i]
                          << ", omega = " << omega[i] << ", tau = " << tau_clump[i]
                          << ", A = " << A << ", B = " << B << ", C = " << C
                          << std::endl;
        }
        return g_clump_val;
    }

    /**
     * Calculate effective dust optical properties based on Inoue 2005
     * @param k_val Extinction coefficient per unit mass
     * @param omega_val Scattering albedo
     * @param g_val Asymmetry parameter
     */
    inline void
    MegaGrainApproximationOneAge(const val& k_val, const val& omega_val,
                                 const val& g_val, double n0_cnm) noexcept {
        // extern double n0_cnm_;  // Declare it as an external variable
        const auto     n0_CNM    = n0_cnm; // default: 1e3 1e4 1e5(kano)
        // The optical depth of a clump relative to the interclump medium
        // [g / cm3] Regard CNM as clump
        const auto     n_CNM     = n0_CNM * std::pow(PKB / T_CNM_MGA, 1.0 / INDEX_CNM);
        const auto RHO_CLUMP = MEAN_ATOMIC_WEIGHT * m_PROTON * n_CNM;
        const auto FRACTION_CLUMP = (n_HOMOGENEOUS - n_WNM) / (n_CNM - n_WNM);

        const auto tau_clump_val   = (RHO_CLUMP - RHO_ICM) * k_val * D_ * R_CLUMP;
        std::cout << "Rho_cl = " << RHO_CLUMP << ", Rho_ICM = " << RHO_ICM << ", R_CLUMP = " << R_CLUMP << ", f_cl = " << FRACTION_CLUMP
        << std::endl;
        // Interaction probability against the parallel light by a sphere with an optical depth
        const auto Pint_clump_val  = PintVal(tau_clump_val);
        // The extinction coefficient per unit length of the medium by clumps (mega-grains)
        const auto kappa_clump_val = 0.75 * FRACTION_CLUMP * Pint_clump_val / R_CLUMP;
        // The effective extinction coefficient per unit length on the interclump medium
        const auto kappa_val       = k_val * D_ * RHO_ICM;
        // The effective extinction coefficient per unit length in the two-phase medium
        kappa_eff_val_ = kappa_clump_val + kappa_val;

        // The photon escape probability from a sphere in which isotropic sources and dust whose scattering is isotropic distribute uniformly
        Pesc_val_ = ProbabilityEscape(tau_clump_val, omega_val, Pint_clump_val);
        // The albedo of a clump
        const auto omega_clump_val = omega_val * Pesc_val_;
        // The effective albedo
        omega_eff_val_ =
        (omega_clump_val * kappa_clump_val + omega_val * kappa_val) / kappa_eff_val_;

        // The asymmetry parameter of a clump
        const auto g_clump_val = GClumpVal(tau_clump_val, g_val, omega_val);
        // The effective asymmetry parameter
        g_eff_val_ = (g_clump_val * kappa_clump_val + g_val * kappa_val) / kappa_eff_val_;
    }

    // 計算で使う mu = cos(theta) の配列を求める 1から-1
    [[nodiscard]] static inline val MuVal() noexcept {
        auto      mu_val = val(N_ANGLE);
        for (auto i      = std::size_t(0); i < N_ANGLE; ++i) {
            mu_val[i] = 1 - static_cast<double>(i) * dMU;
        }
        return mu_val;
    }

    [[nodiscard]] static val Thetaval() noexcept {
        auto      theta_val = val(N_ANGLE);
        for (auto i         = std::size_t(0); i < N_ANGLE; ++i) {
            const auto mu = 1 - static_cast<double>(i) * dMU;
            theta_val[i] = std::acos(mu);
        }
        return theta_val;
    }

    // Henyey & Greenstein 1941
    [[nodiscard]] inline val2 PhaseFunctionVal(double g) const noexcept {
        auto phase_val2 = val2(val(N_ANGLE), N_ANGLE);

        for (auto i = std::size_t(0); i < N_ANGLE; ++i) { //
            for (auto j = std::size_t(0); j < N_ANGLE; ++j) { // 入射角
                auto mu = mu_val_[static_cast<size_t>(std::abs(
                static_cast<int>(j) - static_cast<int>(i)))];

                mu = mu_val_[i] * mu_val_[j] + std::sqrt(1 - mu_val_[i] * mu_val_[i]) *
                                               std::sqrt(1 - mu_val_[j] * mu_val_[j]);
                const auto denominator_inner_value = 1 + g * g - 2 * g * mu;
                phase_val2[i][j] =
                0.5 * (1 - g * g) / (denominator_inner_value * std::sqrt(denominator_inner_value));
            }

            auto      sum = 0.0;
            for (auto j   = std::size_t(1); j < N_ANGLE; ++j) {
                sum += 0.5 * dMU * (phase_val2[i][j] + phase_val2[i][j - 1]);
            }
            // std::cout << "i = " << i << ", sum = " << sum << std::endl;
            // i = 0 のときの合計が1にならないので規格化しておく
            if (std::abs(sum - 1) > DBL_MIN) {
                phase_val2[i] = phase_val2[i] / sum;
            }
        }
        return phase_val2;
    }


    [[nodiscard]] static inline val Zval(const FreeParameter& free_params) noexcept {
        auto      z_val = val(N_Z);
        for (auto i     = std::size_t(0); i < N_Z; ++i) {
            z_val[i] = free_params.h_dust_ - free_params.dz_ * static_cast<double>(i);
        }
        return z_val;
    }

    // Source function の星による放射部分を計算
    [[nodiscard]] static inline val
    SourceFunctionFromStar(double f_young_star, double kappa, double Pesc,
                           const FreeParameter& free_params) noexcept {
        if (f_young_star < DBL_MIN) return val(N_TAU);
        const auto z_val      = Zval(free_params);
        const auto h_dust     = free_params.h_dust_;
        const auto old_star   = XI * std::exp(-XI * z_val / h_dust) / (2 * h_dust);
        const auto young_star = Pesc / (2 * h_dust);
        return ((1 - f_young_star) * old_star + f_young_star * young_star) / kappa;
    }

    [[nodiscard]] inline double IStar(double f_young_star) const noexcept {
        // f_young_star == 0 のときはまだ星が作られていない
        if (f_young_star < DBL_MIN) return 0.0;
        return (1 - f_young_star) * 0.5 * std::exp(-XI);
    }

    // scattering term in source function
    [[nodiscard]] inline val ScatteringEmissivity(const val2& I_val2,
                                                  const val& phase_function_val) const noexcept {
        auto scattering_emissivity_val = val(N_TAU);

        for (auto i = std::size_t(1); i < N_ANGLE; ++i) {
            scattering_emissivity_val += 0.5 * dMU * (I_val2[i] * phase_function_val[i]
                                                      + I_val2[i - 1] * phase_function_val[i - 1]);
        }
        return scattering_emissivity_val;
    }

    [[nodiscard]] inline double
    IntensityOneStep(double d_tau_mu, double I0, double S0, double S1) const noexcept {
        const auto exp_d_tau_mu = std::exp(-d_tau_mu);
        const auto inv_d_tau_mu = 1 / d_tau_mu;
        return d_tau_mu > 1e-3 ? I0 * exp_d_tau_mu
                                 + (inv_d_tau_mu - exp_d_tau_mu * inv_d_tau_mu - exp_d_tau_mu) * S0
                                 + (1 - inv_d_tau_mu + exp_d_tau_mu * inv_d_tau_mu) * S1
               // small tau case
                               : I0 * exp_d_tau_mu
                                 + (0.5 * d_tau_mu - d_tau_mu * d_tau_mu / 3) * S0
                                 + (0.5 * d_tau_mu - d_tau_mu * d_tau_mu / 6) * S1;
    }

    // 輻射輸送の計算が終わるかどうかの判断をする
    // Source function がすべての角度と optical depth においてほとんど変化しなくっていたら終了
    [[nodiscard]] inline bool
    IsConvergence(const val2& old_S_val2, const val2& S_val2) const noexcept {
        for (auto i = std::size_t(0); i < N_ANGLE; ++i) {
            const auto diff_val = std::abs((S_val2[i] - old_S_val2[i]) / S_val2[i]);
            if (diff_val.max() > Ng_CRITERION) return false;
        }
        return true;
    }

    [[nodiscard]] inline val2 InitializeIVal2(double I_star) const noexcept {
        auto      I_val2 = val2(val(N_TAU), N_ANGLE);
        // 初期条件を設定
        for (auto i      = std::size_t(0); i < N_ANGLE; ++i) {
            // 表面での下向きの放射は外側にある古い星からの放射のみ
            if (mu_val_[i] < 0.0) I_val2[i][0] = I_star / (-mu_val_[i]);
        }
        return I_val2;
    }

    [[nodiscard]] inline val2 InitializeSourceFunction(const val& S_star_val) const noexcept {
        auto      S_val2 = val2(val(N_TAU), N_ANGLE);
        // 初期条件を設定
        for (auto i      = std::size_t(0); i < N_ANGLE; ++i) {
            // 散乱を入れない星からの放射を最初の source function にする
            S_val2[i] = S_star_val;
        }
        return S_val2;
    }

    [[nodiscard]] inline val2 RadiativeTransferSolver(double d_tau, double omega, double I_star,
                                                      const val& S_star_val,
                                                      const val2& phi_val2) const noexcept {
        auto I_val2 = InitializeIVal2(I_star);
        auto S_val2 = InitializeSourceFunction(S_star_val);

        auto Ng_count  = std::size_t(0);
        auto Ng_S_val3 = val3(4);

        for (auto n = std::size_t(0); n < N_Ng_MAX; ++n) {
            // Down stream
            for (auto i = N_ANGLE - 1, i_end = N_ANGLE / 2; i >= i_end; --i) {
                const auto d_tau_mu = d_tau / (-mu_val_[i]);
                for (auto  ti       = std::size_t(1); ti < N_TAU; ++ti) {
                    I_val2[i][ti] = IntensityOneStep(d_tau_mu, I_val2[i][ti - 1],
                                                     S_val2[i][ti - 1], S_val2[i][ti]);
                }
            }

            // Up stream
            for (int i = N_ANGLE / 2 - 1; i >= 0; --i) {
                const auto i_angle = static_cast<std::size_t >(i);
                // Mirror boundary condition
                I_val2[i_angle][N_TAU - 1] = I_val2[N_ANGLE - 1 - i_angle][N_TAU - 1];
                const auto d_tau_mu = d_tau / mu_val_[i_angle];

                for (int ti = N_TAU - 2; ti >= 0; --ti) {
                    const auto i_tau = static_cast<std::size_t>(ti);
                    I_val2[i_angle][i_tau] = IntensityOneStep(d_tau_mu, I_val2[i_angle][i_tau + 1],
                                                              S_val2[i_angle][i_tau + 1],
                                                              S_val2[i_angle][i_tau]);
                }
            }

            // calculate source function
            const auto old_S_val2 = S_val2;

            for (auto i = std::size_t(0); i < N_ANGLE; ++i) {
                S_val2[i] = S_star_val + 0.5 * omega * ScatteringEmissivity(I_val2, phi_val2[i]);
            }

            if (IsConvergence(old_S_val2, S_val2)) return I_val2;

            Ng_S_val3[Ng_count] = S_val2;
            ++Ng_count;
            if (Ng_count == 4) {
                S_val2 = NgAcc2(Ng_S_val3);
                Ng_S_val3[0] = S_val2;
                Ng_count = 1;
            }
        }

        // N_Ng_MAX 回のループを経て収束しなかった場合
        return I_val2;
    }

    [[nodiscard]] inline val2 NgAcc2(const val3& Ng_S_val) const noexcept {
        auto      a1 = 0.0;
        auto      b1 = 0.0;
        auto      b2 = 0.0;
        auto      c1 = 0.0;
        auto      c2 = 0.0;
        for (auto i  = std::size_t(0); i < N_ANGLE; ++i) {
            const auto d0 = Ng_S_val[3][i] - Ng_S_val[2][i];
            const auto d1 = Ng_S_val[3][i] - 2 * Ng_S_val[2][i] + Ng_S_val[1][i];
            const auto d2 = Ng_S_val[3][i] - Ng_S_val[2][i] - Ng_S_val[1][i] + Ng_S_val[0][i];

            a1 += (d1 * d1).sum();
            b1 += (d1 * d2).sum();
            b2 += (d2 * d2).sum();
            c1 += (d0 * d1).sum();
            c2 += (d0 * d2).sum();
        }

        const auto a      = (b2 * c1 - b1 * c2) / (b2 * a1 - b1 * b1);
        const auto b      = (a1 * c2 - b1 * c1) / (b2 * a1 - b1 * b1);
        auto       result = val2(N_ANGLE);
        for (auto  i      = std::size_t(0); i < N_ANGLE; ++i) {
            result[i] = (1 - a - b) * Ng_S_val[3][i] + a * Ng_S_val[2][i] + b * Ng_S_val[1][i];
        }
        return result;
    }

    [[nodiscard]] inline auto TransmissionRate(double I_star, const val2& I_val2) const noexcept {
        // transmission rate val
        auto      T_val = val(N_ANGLE);
        for (auto i     = std::size_t(0); i < N_ANGLE; ++i) {
            if (mu_val_[i] > 0) {
                T_val[i] = mu_val_[i] * I_val2[i][0] + I_star;
            }
        }
        return T_val;
    }

  public:
    /**
     * Read total_mass_*, stellar_spectrum_*, dust_parameter_*
     * @param D Dust to gas mass ratio
     * @param dz Delta z-axis
     * @param f_young_val Young stellar fraction
     * @param k_val Extinction coefficient per unit mass
     * @param omega_val Scattering albedo
     * @param g_val Asymmetry parameter
     */
    inline RadiativeTransfer(double D, double dz, val f_young_val, const val& k_val,
                             const val& omega_val, const val& g_val, double n0_cnm) noexcept:
    D_(D),
    f_young_val_(std::move(f_young_val)),
    mu_val_(MuVal()) {
        MegaGrainApproximationOneAge(k_val, omega_val, g_val, n0_cnm);
        tau_val_ = kappa_eff_val_ * dz;
    }

    ~RadiativeTransfer() = default;

    inline void WriteMGAParameter(const FreeParameter& free_params) {
        // Debug: MGA で計算した effective なダストのパラメータを書き出す
        auto ofs = SEDFile::OfsOpen(MGAParameterFileName(free_params));
        ofs << "Wavelength tau_eff kappa_eff omega_eff g_eff extinction" << std::endl; // header

        const auto extinction_val = tau_val_ / tau_val_[822]; // Divide by V 0.55um

        for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
            ofs << lambda_cm_val_[li] << " " << tau_val_[li] << " " << kappa_eff_val_[li] << " "
                << omega_eff_val_[li] << " " << g_eff_val_[li] << " ";
            if (std::isnan(extinction_val[li])) ofs << 0.0 << std::endl;
            else ofs << extinction_val[li] << std::endl;
        }
        ofs.close();
    }

    [[nodiscard]] inline auto Calculate(const FreeParameter& free_params) const noexcept {
        auto transmission_rate_val2 = val2(val(N_ANGLE), N_LAMBDA);
        //auto transmission_rate_val = val(N_LAMBDA);

        //for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
        for (auto li = INDEX_LYMAN_BREAK; li < N_LAMBDA; ++li) {
            const auto d_tau      = kappa_eff_val_[li] * free_params.dz_;
            const auto S_star_val = SourceFunctionFromStar(f_young_val_[li], kappa_eff_val_[li],
                                                           Pesc_val_[li], free_params);
            const auto phi_val2   = PhaseFunctionVal(g_eff_val_[li]);


            const auto I_star = IStar(f_young_val_[li]);

            const auto I_val2 = RadiativeTransferSolver(d_tau, omega_eff_val_[li],
                                                        I_star, S_star_val, phi_val2);
            transmission_rate_val2[li] = TransmissionRate(I_star, I_val2);
        }
        return transmission_rate_val2;
    }

    inline val
    WriteTransmissionRate(const val2& T_val, const FreeParameter& free_params) {
        auto ofs = my_util::FileUtil::OfsOpen(TransmissionRateFileName(free_params));
        ofs << "Wavelength(cm) T T_average" << std::endl;

        auto T_average_val = val(N_LAMBDA);

        for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
            for (auto i_a       = std::size_t(0); i_a < T_val[li].size() - 1; ++i_a) {
                T_average_val[li] += 0.5 * (T_val[li][i_a] + T_val[li][i_a + 1]) * dMU;
            }
            ofs << lambda_cm_val_[li] << " " << T_val[li][0] << " " << T_average_val[li] << std::endl;
        }
        ofs.close();
        return T_average_val;
    }

};
}
#endif //MKSED_RADIATIVE_TRANSFER_H
