#ifndef DUST_EVOLUTION_H
#define DUST_EVOLUTION_H

#include <valarray>

#include <constants.h>
#include <iostream>
#include <numeric>
#include <cfloat>
#include <mutex>

#include "function.h"
#include "function.h"

namespace asano_model
{
    using val = std::valarray<double>;
    using val2 = std::valarray<val>;
    using val3 = std::valarray<val2>;
    using val4 = std::valarray<val3>;
    using val5 = std::valarray<val4>;

    using vast = std::valarray<std::size_t>;
    using vast2 = std::valarray<vast>;
    using vast3 = std::valarray<vast2>;
    using vast4 = std::valarray<vast3>;

    class DustEvolution
    {
    public:
    // コピーコンストラクタ
    DustEvolution(const DustEvolution& other)
        : ISM_params_(other.ISM_params_),
          dust_radius_(other.dust_radius_),
          m_grain_(other.m_grain_),
          grain_velocity_(other.grain_velocity_),
          m_bin_max_(other.m_bin_max_),
          v_ij_(other.v_ij_),
          M_shocked_(other.M_shocked_),
          j_val_(other.j_val_),
          dec_ratio_(other.dec_ratio_),
          inc_ratio_(other.inc_ratio_),
          distri_inc_ratio_(other.distri_inc_ratio_) {
        // std::mutex はコピーしない
    }
    private:
        ISMParams ISM_params_ = ISMParams("CNM");
        val dust_radius_;
        val2 m_grain_;
        val2 grain_velocity_;
        std::array<std::array<double, N_MAX_DUST_RADIUS>, 2> m_bin_max_{};
        val4 v_ij_;
        val4 M_shocked_;
        vast4 j_val_;
        val4 dec_ratio_;
        val4 inc_ratio_;
        val5 distri_inc_ratio_;
        mutable std::mutex shattering_mutex; // クラス内で定義

        inline double RelativeVelocity(std::size_t omega, double v_i, double v_j) const noexcept
        {
            switch (omega)
            {
            case 0:
                return (v_i + v_j); // front collision
            case 1:
                return (std::fabs(v_i - v_j)); // back collision
            case 2:
                return (v_i); // sideway collision
            case 3:
                return (v_j); // sideway collision
            default:
                return (std::sqrt(v_j * v_j + v_i * v_i - v_j * v_i));
            }
        }

        inline val4 RelativeVelocityVarray() const noexcept
        {
            auto v_ij = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
            {
                for (auto omega = std::size_t(0); omega < 4; ++omega)
                {
                    for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
                    {
                        for (auto j_radius = N_RADIUS_MIN; j_radius < N_MAX_DUST_RADIUS; ++j_radius)
                        {
                            v_ij[carsil][omega][i_radius][j_radius] = RelativeVelocity(omega,
                                                                                       grain_velocity_[i_radius][carsil],
                                                                                       grain_velocity_[j_radius][carsil]);
                        }
                    }
                }
            }
            return v_ij;
        }

        inline val4 AlphaIJ() const noexcept
        {
            auto alpha_ij = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
            {
                for (auto omega = std::size_t(0); omega < 4; ++omega)
                {
                    for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
                    {
                        for (auto j_radius = N_RADIUS_MIN; j_radius < N_MAX_DUST_RADIUS; ++j_radius)
                        {
                            alpha_ij[carsil][omega][i_radius][j_radius] = DT / 4 * M_PI * std::pow(dust_radius_[i_radius] + dust_radius_[j_radius], 2) * v_ij_[carsil][omega][i_radius][j_radius] /
                                                                          (m_grain_[carsil][i_radius] *
                                                                           m_grain_[carsil][j_radius]);
                        }
                    }
                }
            }
            return alpha_ij;
        }

        inline val4 MShoked() const noexcept
        {
            const auto M_1 = GetM1();
            const auto sigma_1 = GetSigma1(M_1);
            auto M_shocked = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
            {
                for (auto omega = std::size_t(0); omega < 4; ++omega)
                {
                    for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
                    {
                        for (auto j_radius = N_RADIUS_MIN; j_radius < N_MAX_DUST_RADIUS; ++j_radius)
                        {
                            // The ratio of the relative velocity of grains and the sound speeds
                            const auto M_r = v_ij_[carsil][omega][i_radius][j_radius] / c0[carsil];
                            // Function of the shock paramters
                            const auto sigma_r = 0.3 * std::pow(s_shat[carsil] + 1 / M_r * (1 + R_COLL) - 0.11, 1.3) / (s_shat[carsil] + 1 / M_r * (1 + R_COLL) - 1);
                            M_shocked[carsil][omega][i_radius][j_radius] =
                                0.5 * m_grain_[carsil][j_radius] * (1.0 + 2.0 * R_COLL) / (std::pow(1.0 + R_COLL, 1.778) * std::pow(sigma_r, 0.1111)) * std::pow(M_r * M_r / (sigma_1[carsil] * M_1[carsil] * M_1[carsil]), 0.8889);
                        }
                    }
                }
            }
            return M_shocked;
        }

        // 質量 bin ごとの境界値の大きい方を求めている
        // 10^0.05倍しているのは、半径のbinの大きい方を考慮することで、m_remがどのmass binに入っているかを説明するため
        inline std::array<std::array<double, N_MAX_DUST_RADIUS>, 2>
        GetBinMaxMass(const val &dust_radius) const noexcept
        {
            auto m_bin_max = std::array<std::array<double, N_MAX_DUST_RADIUS>, 2>{{{0}}};
            for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
            {
                for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
                {
                    m_bin_max[carsil][i] = Volume(
                                               dust_radius[i] * std::pow(10.0, 0.05)) *
                                           RHO_CARSIL[carsil];
                }
            }
            return m_bin_max;
        }

        inline std::array<double, 2> GetM1() const noexcept
        {
            auto M_1 = std::array<double,
                                  2>{0};
            for (auto i = std::size_t(0); i < 2; ++i)
            {
                const auto phi = P1[i] / (RHO_CARSIL[i] * c0[i] * c0[i]);
                M_1[i] = 2 * phi / (1 + std::sqrt(1 + 4 * s_shat[i] * phi));
            }
            return M_1;
        }

        inline std::array<double, 2>
        GetSigma1(const std::array<double, 2> &M_1) const noexcept
        {
            auto sigma_1 = std::array<double,
                                      2>{0};
            for (auto i = std::size_t(0); i < 2; ++i)
            {
                sigma_1[i] = 0.3 * std::pow(s_shat[i] + 1.0 / M_1[i] - 0.11, 1.3) / (s_shat[i] + 1.0 / M_1[i] - 1.0);
            }
            return sigma_1;
        }

        // - がついているのは、元が降順なのでサーチが面倒になることから計算のテクニックとして - にしている
        inline std::array<std::array<double, N_MAX_DUST_RADIUS>, 2>
        GetILeft() const noexcept
        {
            auto i_left = std::array<std::array<double, N_MAX_DUST_RADIUS>, 2>{{{0}}};
            for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
            {
                const auto tmp = -Volume(
                    std::pow(10.0,
                             -SHAT_RADIUS_INDEX_MIN - (0.1 * static_cast<double>(i))) *
                    std::pow(
                        10.0,
                        -0.05));
                for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
                {
                    i_left[carsil][i] = tmp * RHO_CARSIL[carsil];
                }
            }
            return i_left;
        }

        // "2" represents carbon or silicate respectively.
        inline std::array<std::array<std::array<double, N_MAX_DUST_RADIUS>, N_MAX_DUST_RADIUS>, 2>

        GetJLeft() const noexcept
        {
            auto j_left = std::array<std::array<std::array<double, N_MAX_DUST_RADIUS>, N_MAX_DUST_RADIUS>, 2>{{{{{0}}}}};
            for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
            {
                for (auto j = std::size_t(0); j < N_MAX_DUST_RADIUS; ++j)
                {
                    const auto tmp = Volume(
                        std::pow(10, -SHAT_RADIUS_INDEX_MIN - (0.1 * static_cast<double>(i)) + (0.1 * static_cast<double>(j))) * std::pow(10, 0.05));
                    for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
                    {
                        j_left[carsil][i][j] = tmp * RHO_CARSIL[carsil];
                    }
                }
            }
            return j_left;
        }

        inline std::array<std::array<std::array<double, N_MAX_DUST_RADIUS>, N_MAX_DUST_RADIUS>, 2>

        GetMShat() const noexcept
        {
            auto m_shat = std::array<std::array<std::array<double, N_MAX_DUST_RADIUS>, N_MAX_DUST_RADIUS>, 2>{{{{{0}}}}};
            for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
            {
                for (auto k = std::size_t(0); k < N_MAX_DUST_RADIUS; ++k)
                {
                    const auto tmp = Volume(
                        std::pow(10, -SHAT_RADIUS_INDEX_MIN - (0.1 * static_cast<double>(i)) + 0.1 * static_cast<double>(k)));
                    for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
                    {
                        const auto m = tmp * RHO_CARSIL[carsil];
                        m_shat[carsil][i][k] = m * m * std::pow(m, -(X_SHAT + 2.0) / 3.0) * DM;
                    }
                }
            }
            return m_shat;
        }

        // shatteringで形成されたfragmentの最小値が3A以下の場合、ダスト破壊とするための計算
        inline double massratio(double smin, double smax, std::size_t carsil) const noexcept
        {
            static const auto i_left = GetILeft();
            static const auto j_left = GetJLeft();
            static const auto m_shat = GetMShat();

            const auto m_min = Volume(smin) * RHO_CARSIL[carsil];
            const auto m_max = Volume(smax) * RHO_CARSIL[carsil];

            // 計算のテクニックとして、降順の配列をサーチするため、配列要素をすべて負にして、m_minも負にすることで普通にサーチできるようになっている。
            const auto i = static_cast<std::size_t>(std::lower_bound(
                                                        std::begin(i_left[carsil]), std::end(i_left[carsil]), -m_min) -
                                                    std::begin(i_left[carsil]));
            // j_left は昇順なので普通にサーチ
            const auto j = static_cast<std::size_t>(std::lower_bound(
                                                        std::begin(j_left[carsil][i]), std::end(j_left[carsil][i]), m_max) -
                                                    std::begin(j_left[carsil][i]));

            if (i > N_MAX_DUST_RADIUS)
                std::cout << "i = " << i << ", j = " << j << std::endl;
            if (j > N_MAX_DUST_RADIUS)
                std::cout << "i = " << i << ", j = " << j << std::endl;
            if (j < i)
                return 1; // つまり、maximum shattered fragment size < 3A

            const auto m_i = Volume(
                                 std::pow(10.0,
                                          -SHAT_RADIUS_INDEX_MIN - (0.1 * static_cast<double>(i)))) *
                             RHO_CARSIL[carsil];
            // m_shatのうち、3A以下のダスト半径の質量を計算
            auto sum1 = m_i * std::pow(m_i, -(X_SHAT + 2) / 3) * (Volume(std::pow(10, -SHAT_RADIUS_INDEX_MIN - (0.1 * static_cast<double>(i)) + 0.05)) * RHO_CARSIL[carsil] - m_min);
            if (i > 1)
            {
                sum1 += std::accumulate(std::begin(m_shat[carsil][i]) + 1,
                                        std::begin(m_shat[carsil][i]) + i, 0.0);
            }

            auto sum = m_i * std::pow(m_min, -(X_SHAT + 2.0) / 3.0) * (Volume(std::pow(10.0, -SHAT_RADIUS_INDEX_MIN - (0.1 * static_cast<double>(i)) + 0.05)) * RHO_CARSIL[carsil] - m_min);
            sum += std::accumulate(std::begin(m_shat[carsil][i]) + 1,
                                   std::begin(m_shat[carsil][i]) + j, 0.0);

            const auto m_j = Volume(std::pow(10.0,
                                             -SHAT_RADIUS_INDEX_MIN + (0.1 * static_cast<double>(j)))) *
                             RHO_CARSIL[carsil];
            sum += m_j * std::pow(m_j, -(X_SHAT + 2.0) / 3.0) * (m_max - Volume(std::pow(10.0, -SHAT_RADIUS_INDEX_MIN + (0.1 * static_cast<double>(j)) - 0.05)) * RHO_CARSIL[carsil]);
            return sum1 / sum;
        }

        inline void
        MShatMRem(val4 &a_fmin, val4 &a_fmax, val4 &m_shat, val4 &m_rem) noexcept
        {
            for (auto o = std::size_t(0); o < 2; ++o)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            // switching catastrophic and cratering
                            auto ak_fmax = 0.0;
                            auto ak_fmin = 0.0;
                            // auto m_shat = 0.0;
                            if (m_grain_[o][i] < 2.0 * M_shocked_[o][w][i][j])
                            { // catastrophic
                                m_shat[o][w][i][j] = f_catastrophic * m_grain_[o][i];
                                ak_fmax = 0.22 * dust_radius_[i];
                                // ak_fmin = ak_fmax*pow(P1/Pv,1.47);
                            }
                            else
                            { // cratering
                                // m_shat[carsil] = 0.4000777294*M_shocked;//0.4000777294s stand for F_M in Hirashita & Kobayashi (2013)
                                m_shat[o][w][i][j] = F_M * M_shocked_[o][w][i][j];
                                ak_fmax = std::pow(
                                              m_shat[o][w][i][j] * 3.0 * (Z_COLL + 1.0) /
                                                  (RHO_CARSIL[o] * 16.0 * 3.142 * (Z_COLL - 2.0)),
                                              0.3333) /
                                          Z_COLL;
                                if (ak_fmax > dust_radius_[i])
                                    std::cout << ak_fmax << " " << dust_radius_[i] << std::endl;
                                // ak_fmin = ak_fmax*pow(P1/Pv,1.47);
                            }
                            ak_fmin = ak_fmax * std::pow(P1[o] / Pv[o], (Z_COLL + 1.0) / Z_COLL);
                            m_rem[o][w][i][j] = m_grain_[o][i] - m_shat[o][w][i][j];

                            // fragment のサイズのすり合わせ
                            if (ak_fmin < 3.e-8)
                            { // 計算するダストサイズの範囲によって変更が必要
                                // if(ak_fmin[carsil] < 5.e-7){
                                a_fmin[o][w][i][j] = 3.e-8;
                                // 全m_shatのうち、何割が3A以下の半径のダストが担っているかを計算。このサイズのダストはダストではないので破壊されたとする。(2012/11/11)
                                m_shat[o][w][i][j] *= 1.0 - massratio(ak_fmin, ak_fmax, o);
                            }
                            else
                            {
                                a_fmin[o][w][i][j] = ak_fmin;
                            }
                            if (ak_fmax < 3.e-8)
                            {
                                // if(ak_fmax[carsil] < 5.e-7){
                                a_fmax[o][w][i][j] = 3.01e-8;
                            }
                            else
                            {
                                a_fmax[o][w][i][j] = ak_fmax;
                            }
                        }
                    }
                }
            }
        }

        inline vast4 Jval(const val4 &m_rem) const noexcept
        {
            auto j_val = vast4(vast3(vast2(vast(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto carsil = std::size_t(0); carsil < 2; ++carsil)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            j_val[carsil][w][i][j] = static_cast<std::size_t>(std::lower_bound(
                                                                                  std::begin(m_bin_max_[carsil]) + N_RADIUS_MIN,
                                                                                  std::end(m_bin_max_[carsil]),
                                                                                  m_rem[carsil][w][i][j]) -
                                                                              std::begin(m_bin_max_[carsil]));
                            if (j_val[carsil][w][i][j] > N_MAX_DUST_RADIUS - 1)
                                j_val[carsil][w][i][j] = N_MAX_DUST_RADIUS - 1;
                        }
                    }
                }
            }
            return j_val;
        }

        inline val5 Distri(const val4 &a_fmin, const val4 &a_fmax) const noexcept
        {
            auto distri = val5(val4(val3(
                                        val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS),
                                    4),
                               2);
            for (auto cs = std::size_t(0); cs < 2; ++cs)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            // 以下m_min,m_maxの値がどのサイズビンに入るかを計算
                            const auto m_min = Volume(a_fmin[cs][w][i][j]) * RHO_CARSIL[cs];
                            const auto m_max = Volume(a_fmax[cs][w][i][j]) * RHO_CARSIL[cs];
                            auto f_min = static_cast<std::size_t>(std::lower_bound(
                                                                      std::begin(m_bin_max_[cs]) + N_RADIUS_MIN,
                                                                      std::end(m_bin_max_[cs]), m_min) -
                                                                  std::begin(m_bin_max_[cs]));
                            if (j_val_[cs][w][i][j] > N_MAX_DUST_RADIUS - 1)
                                f_min = N_MAX_DUST_RADIUS - 1;
                            auto f_max = static_cast<std::size_t>(std::lower_bound(
                                                                      std::begin(m_bin_max_[cs]) + N_RADIUS_MIN,
                                                                      std::end(m_bin_max_[cs]), m_max) -
                                                                  std::begin(m_bin_max_[cs]));
                            if (j_val_[cs][w][i][j] > N_MAX_DUST_RADIUS - 1)
                                f_max = N_MAX_DUST_RADIUS - 1;

                            distri[cs][w][i][j][f_min] = m_grain_[cs][f_min] *
                                                         std::pow(m_grain_[cs][f_min],
                                                                  -(X_SHAT + 2.0) / 3.0) *
                                                         (m_bin_max_[cs][f_min] - m_min);
                            for (auto m = f_min + 1; m < f_max; ++m)
                            {
                                distri[cs][w][i][j][m] = m_grain_[cs][m] * std::pow(m_grain_[cs][m], -(X_SHAT + 2.0) / 3.0) * Volume(dust_radius_[m]) * RHO_CARSIL[cs] * DM;
                            }
                            distri[cs][w][i][j][f_max] = m_grain_[cs][f_max] * std::pow(m_grain_[cs][f_max], -(X_SHAT + 2.0) / 3.0) * (m_max - Volume(dust_radius_[f_max] * std::pow(10.0, -0.05)) * RHO_CARSIL[cs]);
                        }
                    }
                }
            }
            return distri;
        }

        inline val4 GetSum(const val5 &distri) const noexcept
        {
            auto sum = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto cs = std::size_t(0); cs < 2; ++cs)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            sum[cs][w][i][j] = distri[cs][w][i][j].sum();
                        }
                    }
                }
            }
            return sum;
        }

        inline val4 DecrementRatio(const val4 &alpha_ij) const noexcept
        {
            auto dec_ratio = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto o = std::size_t(0); o < 2; ++o)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            dec_ratio[o][w][i][j] = alpha_ij[o][w][i][j] * m_grain_[o][i];
                        }
                    }
                }
            }
            return dec_ratio;
        }

        inline val4 IncrementRatio(const val4 &alpha_ij, const val4 &m_rem) const noexcept
        {
            auto inc_ratio = val4(val3(val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 4), 2);
            for (auto o = std::size_t(0); o < 2; ++o)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            inc_ratio[o][w][i][j] = alpha_ij[o][w][i][j] * m_rem[o][w][i][j];
                        }
                    }
                }
            }
            return inc_ratio;
        }

        inline val5 DistributionIncrementRatio(const val4 &alpha_ij, const val4 &m_shat,
                                               const val4 &sum,
                                               const val5 &distri) const noexcept
        {
            auto distri_inc_ratio = val5(val4(val3(
                                                  val2(val(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS),
                                              4),
                                         2);

            for (auto o = std::size_t(0); o < 2; ++o)
            {
                for (auto w = std::size_t(0); w < 4; ++w)
                {
                    for (auto i = N_RADIUS_MIN; i < N_MAX_DUST_RADIUS; ++i)
                    {
                        for (auto j = N_RADIUS_MIN; j < N_MAX_DUST_RADIUS; ++j)
                        {
                            for (auto k = std::size_t(0); k < N_MAX_DUST_RADIUS; ++k)
                            {
                                distri_inc_ratio[o][w][i][j][k] =
                                    alpha_ij[o][w][i][j] * m_shat[o][w][i][j] * distri[o][w][i][j][k] / sum[o][w][i][j];
                            }
                        }
                    }
                }
            }
            return distri_inc_ratio;
        }

        inline vast3 GetCoagI() const noexcept
        {
            auto m_new_left = val2(val(N_MAX_DUST_RADIUS), 2);

            for (auto o = std::size_t(0); o < 2; ++o)
            {
                for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
                {
                    m_new_left[o][i] = Volume(dust_radius_[i] * std::pow(10.0, 0.05)) * RHO_CARSIL[o];
                }
            }
            auto coag_i_vast3 = vast3(vast2(vast(N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS), 2);

            for (auto o = std::size_t(0); o < 2; ++o)
            {
                for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
                {
                    for (auto j = std::size_t(0); j < N_MAX_DUST_RADIUS; ++j)
                    {
                        const auto m_new = m_grain_[o][i] + m_grain_[o][j];
                        const auto coag_i = static_cast<std::size_t>(std::lower_bound(std::begin(m_new_left[o]), std::end(m_new_left[o]), m_new) - std::begin(m_new_left[o]));
                        coag_i_vast3[o][i][j] = std::clamp(coag_i, N_RADIUS_MIN, N_MAX_DUST_RADIUS - 1);
                    }
                }
            }
            return coag_i_vast3;
        }

        val SilicateNumberDistribution(const val &volume, const val2 &m_dust_val2) const noexcept
        {
            const auto n_sil_per_species = SilicateNumberPerSpecies(volume, m_dust_val2);
            const auto total_n_sil = n_sil_per_species.sum();

            if (total_n_sil < DBL_MIN)
                return val(N_DUST_SPECIES - 1);
            return n_sil_per_species / total_n_sil;
        }

        val SilicateMassDistribution(const val2 &m_dust_val2) const noexcept
        {
            const auto m_sil_per_species = SilicateMassPerSpecies(m_dust_val2);
            const auto total_m_sil = m_sil_per_species.sum();
            // total が 0 のときは全ての半径で 0 を戻す
            if (total_m_sil < DBL_MIN)
                return val(N_DUST_SPECIES - 1);
            return m_sil_per_species / total_m_sil;
        }

        val2 Accretion(const val &M_X, const val &M_dust_X, double M_ISM,
                       const val2 &rho_val2) const noexcept;

        val2 DeltaMByAccretion(double M_ISM, const val &M_X_val, const val &volume,
                               const val &m_sil_ratio_val, const val2 &m_grain_dm_val2,
                               const val2 &m_dust_val2, const val2 &m_carsil_val2) const noexcept;

        val2 Coagulation(const val2 &m_carsil_val2) const noexcept;
        
        
        val2 Shattering(const val2 &m_carsil_val2) const noexcept;

        val2 DeltaMByGrainGrainCollision(double M_ISM, const val &m_sil_ratio_val,
                                         const val2 &m_carsil_val2) const noexcept;

    public:
        DustEvolution(const ISMParams &ISM_params, const val &dust_radius, val2 m_grain);

        val2 Calculation(double M_ISM, const val &M_X_val, const val &volume,
                         const val2 &m_grain_dm_val2, const val2 &m_dust_val2) const noexcept;
    };
}
#endif // DUST_EVOLUTION_H
