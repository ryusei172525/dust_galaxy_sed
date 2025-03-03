//
// Created by 西田和樹 on 2020/09/11.
//
#include "dust_evolution.h"
#include "read_file.h"

namespace asano_model {
    
val2 DustEvolution::Accretion(const val& M_X, const val& M_dust_X, double M_ISM,
                              const val2& rho_val2) const noexcept {
    auto drho_val2 = val2(val(N_MAX_DUST_RADIUS), 2);
    if (!ISM_params_.isACC_) return drho_val2;

    const auto rho_ISM = MU * m_H * ISM_params_.n_; // Mass density of ISM

    for (auto o = std::size_t(0); o < 2; ++o) { // C or Si
        if (M_X[o] < DBL_MIN) continue; // C or Si がないときには飛ばして次にいく
        // Equ.(5.13) in Asano Dron
        // mu_dot = 3*a_dot/a
        const auto mu_dot_val = ALPHA[o] * rho_ISM / (F_GM[o] * RHO_CARSIL[o])
                                * std::sqrt(kb * ISM_params_.T_ / (2.0 * M_PI * M_ATOM[o] / N_A))
                                * M_X[o] / M_ISM * (1.0 - F_GM[o] * M_dust_X[o] / M_X[o])
                                * 3.0 / dust_radius_;

        for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius) {
            drho_val2[o][i_radius] -= HH * YR / DM
                                      * (mu_dot_val[i_radius] * rho_val2[o][i_radius]
                                         - mu_dot_val[i_radius - 1] * rho_val2[o][i_radius - 1]);
        }
    }
    return drho_val2;
}


val2 DustEvolution::DeltaMByAccretion(double M_ISM, const val& M_X_val, const val& volume,
                                      const val& m_sil_ratio_val, const val2& m_grain_dm_val2,
                                      const val2& m_dust_val2,
                                      const val2& m_carsil_val2) const noexcept {
    // total dust mass of key element X
    const auto M_dust_X = {m_carsil_val2[0].sum(), m_carsil_val2[1].sum()};

    const auto drho_carsil_val2 = Accretion(M_X_val, M_dust_X, M_ISM,
                                            RhoCarsil(volume, m_dust_val2));

    auto      dm_acc_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
    for (auto i_radius    = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius) {
        dm_acc_val2[i_radius][0] = drho_carsil_val2[0][i_radius] * m_grain_dm_val2[i_radius][0];
        for (auto i_species = std::size_t(1); i_species < N_DUST_SPECIES - 1; ++i_species) {
            dm_acc_val2[i_radius][i_species] =
            drho_carsil_val2[1][i_radius] * m_grain_dm_val2[i_radius][i_species] *
            m_sil_ratio_val[i_species];
        }
    }
    return dm_acc_val2;
}


val2 DustEvolution::Coagulation(const val2& m_carsil_val2) const noexcept {
    auto dm_carsil_val2 = val2(val(N_MAX_DUST_RADIUS), 2);
    if (ISM_params_.isCoag_ == 0) return dm_carsil_val2;

    static const auto coag_i       = GetCoagI();
    static const auto g_surface = val({75.0, 25.0}); // surface energy per unit area
    static const auto E_young   = val({1.0 * std::pow(10.0, 11.0),
	  5.4 * std::pow(10.0, 11.0)});

    [[maybe_unused]] static constexpr auto v_shat = std::array<double, 2>{1.2e5, 2.7e5};

    for (auto o = std::size_t(0); o < 2; ++o) { // carbon or silicate
        for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i) {
            for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j) {
                const auto m_ij = m_carsil_val2[o][i] * m_carsil_val2[o][j];

                for (auto w = std::size_t(0); w < 4; ++w) {
                    /* switch(media){//CNMで起こるプロセスのダスト速度をYan et al.(2004)かHirashita (2012)を選択している
                     case 0:
                     v_kj[carsil] = relative_vel(omega,grain_velocity[k_radius][carsil],grain_velocity[j_radius][carsil]);
                     break;
                     default:
                     v_kj[carsil] = relative_vel(omega,coag_vel(dust_radius[k_radius],RHO_CARSIL[carsil]),coag_vel(dust_radius[j_radius],RHO_CARSIL[carsil]));
                     break;
                     }
                     */
                    const auto a_reduce = dust_radius_[i] * dust_radius_[j] /
                                          (dust_radius_[i] + dust_radius_[j]);

                    [[maybe_unused]] auto v_coag = 2.141 * F_stick
                                                   * std::sqrt(
                    dust_radius_[i] * dust_radius_[i] - dust_radius_[i] * dust_radius_[j] +
                    dust_radius_[j] * dust_radius_[j]) /
                                                   (dust_radius_[i] + dust_radius_[j]) *
                                                   std::pow(g_surface[o] / a_reduce,
                                                            5.0 / 6.0) /
                                                   std::pow(E_young[o], 1.0 / 3.0) /
                                                   std::sqrt(RHO_CARSIL[o]);

                    //if (v_ij_[o][w][i][j] < std::min(v_coag, v_shat[o])) {
                    //if (v_ij_[o][w][i][j] < v_shat[o]) {
                    dm_carsil_val2[o][i] -= m_ij * dec_ratio_[o][w][i][j];
                    //以下m_newの値がどのサイズビンに入るかを計算(29/11/2011)
                    dm_carsil_val2[o][coag_i[o][i][j]] += m_ij * dec_ratio_[o][w][i][j];
                    //}//v_coagのif文を消すと、全ダスト粒子がcoagulationすることになる。
                }
            }
        }
    }
    return dm_carsil_val2;
}

// メンバ変数アクセス箇所にロック（std::mutex）を追加し、競合を防止する。
// mutable std::mutex shattering_mutex;

// std::mutex shattering_mutex;
// もっとも計算に時間がかかるため、メンバ変数を多用して計算を軽くしている
val2 DustEvolution::Shattering(const val2& m_carsil_val2) const noexcept {
    std::lock_guard<std::mutex> lock(shattering_mutex); // 排他制御を開始


    if (ISM_params_.isShat_ == 0) return val2(val(N_MAX_DUST_RADIUS), 2);
    // shattering が発生する相対速度の閾値
    static constexpr auto v_shat = std::array<double, 2>{1.2e5, 2.7e5};

    // carbon part
    auto dm_car_val = val(N_MAX_DUST_RADIUS);
    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i) {
        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j) {
            const auto m_ij = m_carsil_val2[0][i] * m_carsil_val2[0][j];
            for (auto  w    = std::size_t(0); w < 4; ++w) {
                //shattering occurs if the relative velocity is more than threshold velocity
                if (v_ij_[0][w][i][j] > v_shat[0]) {
                    dm_car_val[i] -= m_ij * dec_ratio_[0][w][i][j];
                    dm_car_val[j_val_[0][w][i][j]] += m_ij * inc_ratio_[0][w][i][j];
                    // 細かく砕けたダストを分配するベクトル演算をしている
                    dm_car_val += m_ij * distri_inc_ratio_[0][w][i][j];
                }
            }
        }
    }

    auto dm_sil_val = val(N_MAX_DUST_RADIUS);
    // silicate part
    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i) {
        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j) {
            const auto m_ij = m_carsil_val2[1][i] * m_carsil_val2[1][j];
            for (auto  w    = std::size_t(0); w < 4; ++w) {
                //shattering occurs if the relative velocity is more than threshold velocity
                if (v_ij_[1][w][i][j] > v_shat[1]) {
                    dm_sil_val[i] -= m_ij * dec_ratio_[1][w][i][j];
                    dm_sil_val[j_val_[1][w][i][j]] += m_ij * inc_ratio_[1][w][i][j];
                    // 細かく砕けたダストを分配するベクトル演算をしている
                    dm_sil_val += m_ij * distri_inc_ratio_[1][w][i][j];
                }
            }
        }
    }
    return {dm_car_val, dm_sil_val};
}

val2 DustEvolution::DeltaMByGrainGrainCollision(double M_ISM, const val& m_sil_ratio_val,
                                                const val2& m_carsil_val2) const noexcept {
    // なぜ n_gas をかけるのかはわからない。
    const auto n_gas = MU * m_H * ISM_params_.n_ / M_ISM;

    // carbon or silicate しか衝突のパラメータなどがわからないので、その二つに分ける。
    // また、carbon と silicate は衝突しないとこにしている
    const auto dm_coag_carsil_val2 = Coagulation(m_carsil_val2);
    const auto dm_shat_carsil_val2 = Shattering(m_carsil_val2);
    const auto dm_collision_carsil_val2 = dm_coag_carsil_val2 + dm_shat_carsil_val2;

    auto dm_collision_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
    for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius) {
        dm_collision_val2[i_radius][0] = dm_collision_carsil_val2[0][i_radius] * n_gas;
        for (auto i_species = std::size_t(1); i_species < N_DUST_SPECIES - 1; ++i_species) {
            dm_collision_val2[i_radius][i_species] =
            dm_collision_carsil_val2[1][i_radius] * m_sil_ratio_val[i_species] * n_gas;
        }
    }
    return dm_collision_val2;
}

val2 DustEvolution::Calculation(double M_ISM, const val& M_X_val, const val& volume,
                                const val2& m_grain_dm_val2,
                                const val2& m_dust_val2) const noexcept {
    const auto m_carsil_val2     = DustMassCarSil(m_dust_val2);
    const auto m_sil_ratio_val   = SilicateMassDistribution(m_dust_val2);
    auto       dm_acc_val2       = val2();
    auto       dm_collision_val2 = val2();
    // #pragma omp parallel sections
    //     {
    // #pragma omp section
    // dm_acc_val2 = DeltaMByAccretion(M_ISM, M_X_val, volume, m_sil_ratio_val,
    //                                     m_grain_dm_val2, m_dust_val2, m_carsil_val2);
    // #pragma omp section
    // dm_collision_val2 = DeltaMByGrainGrainCollision(M_ISM, m_sil_ratio_val,
	//                                               m_carsil_val2);

    // }

    // インスタンスをスレッドごとに分離
    // 各スレッドが DustEvolution の別々のインスタンスを使用するようにする。
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // DustEvolution dust_evolution_instance1 = *this; // スレッド専用のインスタンス
            DustEvolution dust_evolution_instance1(*this); // コピーコンストラクタを使用
            dm_acc_val2 = dust_evolution_instance1.DeltaMByAccretion(M_ISM, M_X_val, volume, m_sil_ratio_val,
                                                                    m_grain_dm_val2, m_dust_val2, m_carsil_val2);
        }
        #pragma omp section
        {
            // DustEvolution dust_evolution_instance2 = *this; // スレッド専用のインスタンス
            DustEvolution dust_evolution_instance2(*this); // コピーコンストラクタを使用
            dm_collision_val2 = dust_evolution_instance2.DeltaMByGrainGrainCollision(M_ISM, m_sil_ratio_val,
                                                                                    m_carsil_val2);
        }
    }
    // これにより、メンバ変数が各スレッドで独立するため競合が解消されます。

    // ISM ごとの fraction をかけてから戻す
   return val(ISM_params_.f_, N_MAX_DUST_RADIUS) * (dm_acc_val2 + dm_collision_val2);
}


DustEvolution::DustEvolution(const ISMParams& ISM_params, const val& dust_radius, val2 m_grain)
: ISM_params_(ISM_params),
  dust_radius_(dust_radius), m_grain_(std::move(m_grain)),
  grain_velocity_(ReadGrainVelocityFile(ISM_params.name_)),
  m_bin_max_(GetBinMaxMass(dust_radius)), v_ij_(RelativeVelocityVarray()),
  M_shocked_(MShoked()) {
    const auto alpha_ij = AlphaIJ();
    auto       a_fmin   = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
    auto       a_fmax   = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
    auto       m_rem    = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
    auto       m_shat   = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
    MShatMRem(a_fmin, a_fmax, m_shat, m_rem);
    j_val_ = Jval(m_rem);
    const auto distri = Distri(a_fmin, a_fmax);
    const auto sum    = GetSum(distri);
    dec_ratio_        = DecrementRatio(alpha_ij);
    inc_ratio_        = IncrementRatio(alpha_ij, m_rem);
    distri_inc_ratio_ = DistributionIncrementRatio(alpha_ij, m_shat, sum, distri);
}
}
