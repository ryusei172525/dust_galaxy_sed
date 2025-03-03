#ifndef ASANO_MODEL_FUNCTION_H
#define ASANO_MODEL_FUNCTION_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <numeric>
#include <tuple>

#include <constants.h>
#include <star_formation_history.h>

namespace asano_model
{
    using val = std::valarray<double>;
    using val2 = std::valarray<val>;
    using val3 = std::valarray<val2>;
    using val4 = std::valarray<val3>;

    using valt = std::valarray<std::size_t>;
    using valt2 = std::valarray<valt>;

    // python, Javascript などの他の言語の map 関数に相当する
    // 非破壊的に関数を適用したコンテナを戻す
    template <class Container, class Transform>
    static inline Container Map(const Container &c, Transform f) noexcept
    {
        auto result = Container(c.size());
        std::transform(begin(c), end(c), std::begin(result), f);
        return result;
    }

    inline double Volume(double a) noexcept
    {
        return 4.0 / 3.0 * M_PI * a * a * a;
    }

    inline val StellarLifeTime() noexcept
    {
        // metallicityは時間に依存ではなく、常にsolar metallicityとする。(2010/12/29)
        // これは、今回採用している星の寿命の式において、使用できる範囲を越えると、問題がある。
        // 同じ星の寿命でも一桁以上変わってしまう。この点はまだ理論的にも不明なので使用しない。
        // 使用できる範囲でならばmetallicity依存はほとんどないので、solar metallicityを使用する。
        const auto lZ = std::log10(0.02);

        const auto a0 = 10.13 + 0.07547 * lZ - 0.008084 * lZ * lZ;
        const auto a1 = -4.424 - 0.7939 * lZ - 0.1187 * lZ * lZ;
        const auto a2 = 1.262 + 0.3385 * lZ + 0.05417 * lZ * lZ;

        auto tau_star_val = val(N_M_STAR);

        for (auto i = std::size_t(0); i < N_M_STAR; ++i)
        {
            const auto m = 1.0 + 0.01 * static_cast<double>(i);
            tau_star_val[i] = std::pow(10, a0 + a1 * std::log10(m + 0.005) +
                                               a2 * std::log10(m + 0.005) * std::log10(m + 0.005));
        }
        return tau_star_val;
    }

    // inline std::valarray<size_t>
    // SearchMetallicityIndex(std::size_t n, std::valarray<size_t> tau, const val& M_Z_val,
    //                        const val& M_gas_val) noexcept {
    //      maximum index number of metallicity is 1000 and minimum one is 5
    //     auto      Z = std::valarray<size_t>(N_M_STAR);
    //     for (auto i = std::size_t(0); i < N_M_STAR; ++i) {
    //         Z[i] = static_cast<std::size_t>(M_Z_val[n - tau[i]] / M_gas_val[n - tau[i]] / Zsun * 1000);
    //         Z[i] = std::clamp(Z[i], std::size_t(5), std::size_t(1000));
    //     }
    //     return Z;
    // }
    //
    // size_t SearchMetallicityIndex(std::size_t n, std::size_t k, const val& m_metal,
    //                               const val& m_gas) noexcept;

    // 時刻がi_ageのときに超新星爆発できる星の質量の下限値のインデックスを探す
    // m_min(t) in Equ. (4.7)-(4.10) of Nishida D thesis
    inline std::size_t MinimumStellarMassIndex(std::size_t i_age, const val &tau_star_val) noexcept
    {
        const auto t = static_cast<double>(i_age) * TIME_BIN_DUST_MODEL;
        for (auto i = std::size_t(0); i < N_M_STAR; ++i)
        {
            if (t - tau_star_val[i] > 0)
                return i;
        }
        return N_M_STAR; // no star can explode a supernova
    }

    // Index of (t-tau_star) on Equ. (4.7)-(4.10) of Nishida D thesis
    inline valt
    StellarBirthIndex(std::size_t i_age, std::size_t i_m_min, const val &tau_star_val) noexcept
    {
        auto i_star_birth_valt = valt(N_M_STAR);

        for (auto i = i_m_min; i < N_M_STAR; ++i)
        {
            i_star_birth_valt[i] =
                i_age - static_cast<std::size_t>(tau_star_val[i] / TIME_BIN_DUST_MODEL);
        }
        return i_star_birth_valt;
    }

    // Index of Z(t-tau_star) in Equ. (4.7)-(4.10) of Nishida D thesis
    inline valt
    StellarMetallicityIndex(std::size_t i_m_min, const valt &i_star_birth_valt, const val &Z_val)
    {
        // maximum index number of metallicity is 1000 and minimum one is 5
        auto i_Z_valt = valt(N_M_STAR);

        for (auto i = i_m_min; i < N_M_STAR; ++i)
        {
            const auto i_Z = static_cast<std::size_t>(Z_val[i_star_birth_valt[i]] / Zsun * 1000);
            i_Z_valt[i] = std::clamp(i_Z, std::size_t(5), std::size_t(1000));
        }
        return i_Z_valt;
    }

    // Equ. (4.7) of Nishida D-thesis
    inline double
    ReturnMass(std::size_t i_m_min, const valt &i_birth_valt, const valt &i_Z_valt, const val &IMF_val,
               const val &SFR_val, const val2 &M_rem_val2) noexcept
    {
        auto R = 0.0;
        const auto dm = 0.01;

        for (auto i = i_m_min; i < N_M_STAR - 1; ++i)
        {
            const auto m_star = 1.0 + dm * static_cast<double>(i); // Msun

            const auto rem1 = M_rem_val2[i_Z_valt[i]][i];
            const auto rem2 = M_rem_val2[i_Z_valt[i + 1]][i + 1];

            R += ((m_star - rem1) * IMF_val[i] * SFR_val[i_birth_valt[i]] + (m_star + dm - rem2) * IMF_val[i + 1] * SFR_val[i_birth_valt[i + 1]]);
        }
        R *= dm * 0.5;
        return R;
    }

    // Equ. (4.9) of Nishida D-thesis
    inline std::tuple<double, double, double>
    MetalYield(std::size_t i_m_min, const valt &i_birth_valt, const valt &i_Z_valt, const val &IMF_val,
               const val &SFR_val, const val3 &M_Z_per_star_val3) noexcept
    {
        auto Y_Z_val = val(3);
        const auto dm = 0.01;

        for (auto i_species = std::size_t(0); i_species < 3; ++i_species)
        {
            for (auto i = i_m_min; i < N_M_STAR - 1; ++i)
            {
                const auto M_Z1 = M_Z_per_star_val3[i_Z_valt[i]][i][i_species];
                const auto M_Z2 = M_Z_per_star_val3[i_Z_valt[i + 1]][i + 1][i_species];

                Y_Z_val[i_species] += (M_Z1 * IMF_val[i] * SFR_val[i_birth_valt[i]] + M_Z2 * IMF_val[i + 1] * SFR_val[i_birth_valt[i + 1]]);
            }
        }
        Y_Z_val *= dm * 0.5;
        return std::make_tuple(Y_Z_val[0], Y_Z_val[1], Y_Z_val[2]);
    }

    inline double GetStarBurst(std::size_t i_age, double m_burst_rem) noexcept
    {
        if (!isBURST)
            return 0;
        const auto age = static_cast<double>(i_age) * TIME_BIN_DUST_MODEL;
        if (age < T_BURST)
            return 0;
        return m_burst_rem * std::exp(-(age - T_BURST) / T_E_FOLDING) / T_E_FOLDING *
               TIME_BIN_DUST_MODEL;
    }

    inline double MetalMass(double M_gas, double SFR, double M_Z, double Y_Z) noexcept
    {
        return M_Z + TIME_BIN_DUST_MODEL * (-M_Z / M_gas * SFR + Y_Z);
    }

    inline std::size_t AGBMetalicityIteration(double Z) noexcept
    {
        if (AGB_DATA_TYPE == 0)
        {
            if (Z <= 0.05)
                return 0;
            else if (Z <= 0.1)
                return 1;
            else if (Z <= 0.2)
                return 2;
            else if (Z <= 0.4)
                return 3;
            else if (Z <= 0.75)
                return 4;
            else
                return 5;
        }
        else
        {
            if (Z <= 0.2)
                return 0;
            else
                return 1;
        }
    }

    inline double SNR(std::size_t i_m_min, const valt &i_birth_valt, const val &IMF_val,
                      const val &SFR_val) noexcept
    {
        auto R = 0.0;

        // 8Msun-40Msun
        for (auto i = std::max(i_m_min, std::size_t(700)); i < N_M_STAR - 1; ++i)
        {
            R += IMF_val[i] * SFR_val[i_birth_valt[i]] + IMF_val[i + 1] * SFR_val[i_birth_valt[i + 1]];
        }
        R *= 0.01 / 2;
        return R;
    }

    inline val DeltaRadius(const val &dust_radius) noexcept
    {
        auto da = val(N_MAX_DUST_RADIUS);
        for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
        {
            da[i] = std::pow(10.0, std::log10(dust_radius[i])) * (std::pow(10.0, 0.05) - std::pow(10.0, -0.05));
        }
        return da;
    }

    inline std::array<std::size_t, 2>
    SNStellarMassIteration(std::size_t i_stellar_mass) noexcept
    {
        if (i_stellar_mass < 1399)
        { // SN part
            return {0, 0};
        }
        else if (i_stellar_mass == 1399)
        { // 14.99Msun - Ms - 15.0Msun
            return {0, 1};
        }
        else if (i_stellar_mass < 2099)
        {
            return {1, 1};
        }
        else if (i_stellar_mass == 2099)
        { // 21.99Msun - Ms - 22.0Msun
            return {1, 2};
        }
        else if (i_stellar_mass < 2899)
        {
            return {2, 2};
        }
        else if (i_stellar_mass == 2899)
        { // 29.99Msun - Ms - 30.0Msun
            return {2, 3};
        }
        else if (i_stellar_mass <= 3899)
        {
            return {3, 3};
        }
        return {0, 0};
    }

    inline valt2
    AGBMetallicityIterationVal2(std::size_t i_m_min, const valt &i_birth_valt, const val &Z_val)
    {
        auto i_z_AGB_val2 = valt2(valt(2), 700);

        for (auto i_m = i_m_min; i_m < 700; ++i_m)
        { // AGB:0-699
            i_z_AGB_val2[i_m][0] = AGBMetalicityIteration(Z_val[i_birth_valt[i_m] + 1]);
            i_z_AGB_val2[i_m][1] = AGBMetalicityIteration(Z_val[i_birth_valt[i_m + 1] + 1]);
        }
        return i_z_AGB_val2;
    }

    // 星質量ごとに使用するデータが変わるのでその変わり目で、iteration number を変えていく
    inline std::array<std::size_t, 2> AGBStellarMassIteration(std::size_t i_m) noexcept
    {
        auto i_AGB = std::array<std::size_t, 2>();
        if (AGB_DATA_TYPE == 0)
        { // Zhukovska et al.(2008)
            if (i_m < 5)
            { // 1.0Msunのデータを、星質量<1.05Msunの星が形成するダストの代表値としている
                i_AGB[0] = 0;
                i_AGB[1] = 0;
            }
            else if (i_m == 5)
            {
                i_AGB[0] = 0;
                i_AGB[1] = 1;
            }
            else if (i_m < 15)
            {
                i_AGB[0] = 1;
                i_AGB[1] = 1;
            }
            else if (i_m == 15)
            {
                i_AGB[0] = 1;
                i_AGB[1] = 2;
            }
            else if (i_m < 22)
            {
                i_AGB[0] = 2;
                i_AGB[1] = 2;
            }
            else if (i_m == 22)
            {
                i_AGB[0] = 2;
                i_AGB[1] = 3;
            }
            else if (i_m < 27)
            {
                i_AGB[0] = 3;
                i_AGB[1] = 3;
            }
            else if (i_m == 27)
            {
                i_AGB[0] = 3;
                i_AGB[1] = 4;
            }
            else if (i_m < 35)
            {
                i_AGB[0] = 4;
                i_AGB[1] = 4;
            }
            else if (i_m == 35)
            {
                i_AGB[0] = 4;
                i_AGB[1] = 5;
            }
            else if (i_m < 45)
            {
                i_AGB[0] = 5;
                i_AGB[1] = 5;
            }
            else if (i_m == 45)
            {
                i_AGB[0] = 5;
                i_AGB[1] = 6;
            }
            else if (i_m < 55)
            {
                i_AGB[0] = 6;
                i_AGB[1] = 6;
            }
            else if (i_m == 55)
            {
                i_AGB[0] = 6;
                i_AGB[1] = 7;
            }
            else if (i_m < 65)
            {
                i_AGB[0] = 7;
                i_AGB[1] = 7;
            }
            else if (i_m == 65)
            {
                i_AGB[0] = 7;
                i_AGB[1] = 8;
            }
            else if (i_m < 75)
            {
                i_AGB[0] = 8;
                i_AGB[1] = 8;
            }
            else if (i_m == 75)
            {
                i_AGB[0] = 8;
                i_AGB[1] = 9;
            }
            else if (i_m < 85)
            {
                i_AGB[0] = 9;
                i_AGB[1] = 9;
            }
            else if (i_m == 85)
            {
                i_AGB[0] = 9;
                i_AGB[1] = 10;
            }
            else if (i_m < 95)
            {
                i_AGB[0] = 10;
                i_AGB[1] = 10;
            }
            else if (i_m == 95)
            {
                i_AGB[0] = 10;
                i_AGB[1] = 11;
            }
            else if (i_m < 105)
            {
                i_AGB[0] = 11;
                i_AGB[1] = 11;
            }
            else if (i_m == 105)
            {
                i_AGB[0] = 11;
                i_AGB[1] = 12;
            }
            else if (i_m < 115)
            {
                i_AGB[0] = 12;
                i_AGB[1] = 12;
            }
            else if (i_m == 115)
            {
                i_AGB[0] = 12;
                i_AGB[1] = 13;
            }
            else if (i_m < 125)
            {
                i_AGB[0] = 13;
                i_AGB[1] = 13;
            }
            else if (i_m == 125)
            {
                i_AGB[0] = 13;
                i_AGB[1] = 14;
            }
            else if (i_m < 135)
            {
                i_AGB[0] = 14;
                i_AGB[1] = 14;
            }
            else if (i_m == 135)
            {
                i_AGB[0] = 14;
                i_AGB[1] = 15;
            }
            else if (i_m < 145)
            {
                i_AGB[0] = 15;
                i_AGB[1] = 15;
            }
            else if (i_m == 145)
            {
                i_AGB[0] = 15;
                i_AGB[1] = 16;
            }
            else if (i_m < 175)
            {
                i_AGB[0] = 16;
                i_AGB[1] = 16;
            }
            else if (i_m == 175)
            {
                i_AGB[0] = 16;
                i_AGB[1] = 17;
            }
            else if (i_m < 225)
            {
                i_AGB[0] = 17;
                i_AGB[1] = 17;
            }
            else if (i_m == 225)
            {
                i_AGB[0] = 17;
                i_AGB[1] = 18;
            }
            else if (i_m < 275)
            {
                i_AGB[0] = 18;
                i_AGB[1] = 18;
            }
            else if (i_m == 275)
            {
                i_AGB[0] = 18;
                i_AGB[1] = 19;
            }
            else if (i_m < 301)
            {
                i_AGB[0] = 19;
                i_AGB[1] = 19;
            }
            else if (i_m == 301)
            {
                i_AGB[0] = 19;
                i_AGB[1] = 20;
            }
            else if (i_m < 325)
            {
                i_AGB[0] = 20;
                i_AGB[1] = 20;
            }
            else if (i_m == 325)
            {
                i_AGB[0] = 20;
                i_AGB[1] = 21;
            }
            else if (i_m < 375)
            {
                i_AGB[0] = 21;
                i_AGB[1] = 21;
            }
            else if (i_m == 375)
            {
                i_AGB[0] = 21;
                i_AGB[1] = 22;
            }
            else if (i_m < 425)
            {
                i_AGB[0] = 22;
                i_AGB[1] = 22;
            }
            else if (i_m == 425)
            {
                i_AGB[0] = 22;
                i_AGB[1] = 23;
            }
            else if (i_m < 475)
            {
                i_AGB[0] = 23;
                i_AGB[1] = 23;
            }
            else if (i_m == 475)
            {
                i_AGB[0] = 23;
                i_AGB[1] = 24;
            }
            else if (i_m < 525)
            {
                i_AGB[0] = 24;
                i_AGB[1] = 24;
            }
            else if (i_m == 525)
            {
                i_AGB[0] = 24;
                i_AGB[1] = 25;
            }
            else if (i_m < 575)
            {
                i_AGB[0] = 25;
                i_AGB[1] = 25;
            }
            else if (i_m == 575)
            {
                i_AGB[0] = 25;
                i_AGB[1] = 26;
            }
            else if (i_m < 700)
            {
                i_AGB[0] = 26;
                i_AGB[1] = 26;
            }
        }
        else
        {
            if (i_m < 75)
            { // 1.5Msunのデータを、星質量<1.75Msunの星が形成するダストの代表値としている
                i_AGB[0] = 0;
                i_AGB[1] = 0;
            }
            else if (i_m == 75)
            {
                i_AGB[0] = 0;
                i_AGB[1] = 1;
            }
            else if (i_m < 125)
            {
                i_AGB[0] = 1;
                i_AGB[1] = 1;
            }
            else if (i_m == 125)
            {
                i_AGB[0] = 1;
                i_AGB[1] = 2;
            }
            else if (i_m < 175)
            {
                i_AGB[0] = 2;
                i_AGB[1] = 2;
            }
            else if (i_m == 175)
            {
                i_AGB[0] = 2;
                i_AGB[1] = 3;
            }
            else if (i_m < 225)
            {
                i_AGB[0] = 3;
                i_AGB[1] = 3;
            }
            else if (i_m == 225)
            {
                i_AGB[0] = 3;
                i_AGB[1] = 4;
            }
            else if (i_m < 275)
            {
                i_AGB[0] = 4;
                i_AGB[1] = 4;
            }
            else if (i_m == 275)
            {
                i_AGB[0] = 4;
                i_AGB[1] = 5;
            }
            else if (i_m < 325)
            {
                i_AGB[0] = 5;
                i_AGB[1] = 5;
            }
            else if (i_m == 325)
            {
                i_AGB[0] = 5;
                i_AGB[1] = 6;
            }
            else if (i_m < 375)
            {
                i_AGB[0] = 6;
                i_AGB[1] = 6;
            }
            else if (i_m == 375)
            {
                i_AGB[0] = 6;
                i_AGB[1] = 7;
            }
            else if (i_m < 425)
            {
                i_AGB[0] = 7;
                i_AGB[1] = 7;
            }
            else if (i_m == 425)
            {
                i_AGB[0] = 7;
                i_AGB[1] = 8;
            }
            else if (i_m < 475)
            {
                i_AGB[0] = 8;
                i_AGB[1] = 8;
            }
            else if (i_m == 475)
            {
                i_AGB[0] = 8;
                i_AGB[1] = 9;
            }
            else if (i_m < 525)
            {
                i_AGB[0] = 9;
                i_AGB[1] = 9;
            }
            else if (i_m == 525)
            {
                i_AGB[0] = 9;
                i_AGB[1] = 10;
            }
            else if (i_m < 575)
            {
                i_AGB[0] = 10;
                i_AGB[1] = 10;
            }
            else if (i_m == 575)
            {
                i_AGB[0] = 10;
                i_AGB[1] = 11;
            }
            else if (i_m < 625)
            {
                i_AGB[0] = 11;
                i_AGB[1] = 11;
            }
            else if (i_m == 625)
            {
                i_AGB[0] = 11;
                i_AGB[1] = 12;
            }
            else if (i_m < 675)
            {
                i_AGB[0] = 12;
                i_AGB[1] = 12;
            }
            else if (i_m == 675)
            {
                i_AGB[0] = 12;
                i_AGB[1] = 13;
            }
            else if (i_m < 700)
            {
                i_AGB[0] = 13;
                i_AGB[1] = 13;
            }
        }
        return i_AGB;
    }

    inline std::array<std::array<std::size_t, 2>, 700>
    AGBStellarMassIterationArr2(std::size_t i_m_min) noexcept
    {
        auto i_m_AGB_arr2 = std::array<std::array<std::size_t, 2>, 700>();

        for (auto i_m = i_m_min; i_m < 700; ++i_m) // AGB:0-699
            i_m_AGB_arr2[i_m] = AGBStellarMassIteration(i_m);
        return i_m_AGB_arr2;
    }

    inline std::array<std::array<std::size_t, 2>, N_M_STAR>
    SNStellarMassIterationArr2(std::size_t i_m_min) noexcept
    {
        auto i_m_SN_arr2 = std::array<std::array<std::size_t, 2>, N_M_STAR>();

        for (auto i_m = std::max(i_m_min, std::size_t(700)); i_m < N_M_STAR; ++i_m) // SN:
            i_m_SN_arr2[i_m] = SNStellarMassIteration(i_m);
        return i_m_SN_arr2;
    }

    inline double
    AGBYield(std::size_t i_m_min, std::size_t i_s, std::size_t i_a, const valt &i_birth_valt,
             const valt2 &i_z_AGB_val2, const std::array<std::array<std::size_t, 2>, 700> &i_m_AGB_ar2,
             const val &IMF_val, const val &SFR_val, const val4 &M_AGB_val4) noexcept
    {
        auto Y_AGB = 0.0;
        const auto dm = 0.01;

        for (auto i_m = i_m_min; i_m < 700; ++i_m)
        { // AGB:0-699
            const auto dust1 = M_AGB_val4[i_z_AGB_val2[i_m][0]][i_m_AGB_ar2[i_m][0]][i_s][i_a];
            const auto dust2 = M_AGB_val4[i_z_AGB_val2[i_m][1]][i_m_AGB_ar2[i_m][1]][i_s][i_a];

            Y_AGB += (dust1 * IMF_val[i_m] * SFR_val[i_birth_valt[i_m]] + dust2 * IMF_val[i_m + 1] * SFR_val[i_birth_valt[i_m + 1]]);
        }
        Y_AGB *= dm * 0.5;
        return Y_AGB;
    }

    inline double
    SNYield(std::size_t i_m_min, std::size_t i_s, std::size_t i_a, const valt &i_birth_valt,
            const std::array<std::array<std::size_t, 2>, N_M_STAR> &i_m_SN_arr2, const val &IMF_val,
            const val &SFR_val, const val3 &M_SN_val3) noexcept
    {
        auto Y_SN = 0.0;
        const auto dm = 0.01;

        for (auto i_m = std::max(i_m_min, std::size_t(700)); i_m < N_M_STAR - 1; ++i_m)
        { // SN:700-3899
            const auto dust1 = M_SN_val3[i_m_SN_arr2[i_m][0]][i_s][i_a];
            const auto dust2 = M_SN_val3[i_m_SN_arr2[i_m][1]][i_s][i_a];
            Y_SN += (dust1 * IMF_val[i_m] * SFR_val[i_birth_valt[i_m]] +
                     dust2 * IMF_val[i_m + 1] * SFR_val[i_birth_valt[i_m + 1]]);
        }
        Y_SN *= dm * 0.5;
        return Y_SN;
    }

    inline val2
    DustYield(std::size_t i_m_min, const valt &i_birth_valt, const val &m_gas, const val &m_metal,
              const val &IMF_val, const val &SFR_val, const val3 &M_SN_val3,
              const val4 &M_AGB_val4) noexcept
    {
        auto Y_d_val2 = val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES);
        const auto Z_val = m_metal / m_gas / Zsun;

        const auto i_m_SN_arr2 = SNStellarMassIterationArr2(i_m_min);
        const auto i_z_AGB_val2 = AGBMetallicityIterationVal2(i_m_min, i_birth_valt, Z_val);
        const auto i_m_AGB_arr2 = AGBStellarMassIterationArr2(i_m_min);

        // 最後のFe3O4はmixed modelの場合のみのため、除外
        for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
        {
            for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
            {
                Y_d_val2[i_s][i_a] = AGBYield(i_m_min, i_s, i_a, i_birth_valt, i_z_AGB_val2,
                                              i_m_AGB_arr2, IMF_val, SFR_val, M_AGB_val4);
                Y_d_val2[i_s][i_a] += SNYield(i_m_min, i_s, i_a, i_birth_valt, i_m_SN_arr2,
                                              IMF_val, SFR_val, M_SN_val3);
            }
        }
        return Y_d_val2;
    }

    inline val2 GrainMass(const val &volume) noexcept
    {
        auto m_grain = val2(val(N_MAX_DUST_RADIUS), 2);
        for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
        {
            m_grain[0][i_radius] = volume[i_radius] * RHO_GRAIN[0];  // carbon dust
            m_grain[1][i_radius] = volume[i_radius] * RHO_GRAIN[10]; // silicate
        }
        return m_grain;
    }

    inline val2 GrainMassDm(const val &volume) noexcept
    {
        auto m_grain_dm = val2(val(N_DUST_SPECIES - 1), N_MAX_DUST_RADIUS);
        for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
        {
            m_grain_dm[i_radius] = volume[i_radius] * RHO_GRAIN * DM;
        }
        return m_grain_dm;
    }

    static inline val2 DustMassCarSil(const val2 &m_dust_val2) noexcept
    {
        auto m_carsil_val2 = val2(val(N_MAX_DUST_RADIUS), 2);

        for (auto i_radius = std::size_t(0); i_radius < N_MAX_DUST_RADIUS; ++i_radius)
        {
            m_carsil_val2[0][i_radius] = m_dust_val2[i_radius][0];
            m_carsil_val2[1][i_radius] = std::accumulate(std::begin(m_dust_val2[i_radius]) + 1,
                                                         std::end(m_dust_val2[i_radius]), 0.0);
        }
        return m_carsil_val2;
    }

    static inline val2 RhoCarsil(const val &volume, const val2 &m_dust_val2) noexcept
    {
        auto rho_carsil_val2 = val2(val(N_MAX_DUST_RADIUS), 2);

        for (auto i_radius = std::size_t(0); i_radius < N_MAX_DUST_RADIUS; ++i_radius)
        {
            rho_carsil_val2[0][i_radius] =
                m_dust_val2[i_radius][0] / RHO_GRAIN[0] / volume[i_radius] / DM;
            for (auto i_species = std::size_t(1); i_species < N_DUST_SPECIES - 1; ++i_species)
            {
                rho_carsil_val2[1][i_radius] +=
                    m_dust_val2[i_radius][i_species] / RHO_GRAIN[i_species] / volume[i_radius] / DM;
            }
        }
        return rho_carsil_val2;
    }

    static inline val SilicateMassPerSpecies(const val2 &m_dust) noexcept
    {
        auto m_sil_per_species = val(N_DUST_SPECIES - 1);
        for (auto i_species = std::size_t(1); i_species < N_DUST_SPECIES - 1; ++i_species)
        {
            for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
            {
                m_sil_per_species[i_species] += m_dust[i_radius][i_species];
            }
        }
        return m_sil_per_species;
    }

    static inline val
    SilicateNumberPerSpecies(const val &volume, const val2 &m_dust) noexcept
    {
        auto n_per_species = val(N_DUST_SPECIES - 1);
        for (auto i_species = std::size_t(1); i_species < N_DUST_SPECIES - 1; ++i_species)
        {
            for (auto i_radius = N_RADIUS_MIN; i_radius < N_MAX_DUST_RADIUS; ++i_radius)
            {
                n_per_species[i_species] += m_dust[i_radius][i_species] * Msun / (volume[i_radius] * RHO_GRAIN[i_species]);
            }
        }
        return n_per_species;
    }

    inline double swept(double gas, double metal) noexcept
    {
        auto metallicity = metal / gas;
        if (metal / gas < pow(10.0, -4.0) * Zsun)
            metallicity = pow(10.0, -4.0) * Zsun;
        else if (metal / gas > Zsun)
            metallicity = Zsun;
        return 1535.0 * pow(n_H, -0.202) * pow((metallicity / Zsun) + 0.039, -0.298);
    }

    inline double
    GasMass(bool is_infall, size_t n, double m_gas_pre, double m_return, double schmidt_index,
            double tau_SF, double tau_infall) noexcept
    {
        const auto m_gas = m_gas_pre +
                           TIME_BIN_DUST_MODEL *
                               (-SFH::SFR(m_gas_pre, schmidt_index, tau_SF, static_cast<double>(n) * TIME_BIN_DUST_MODEL) + m_return);
        if (is_infall == false)
            return m_gas;
        else
            return m_gas +
                   SFH::InfallMass(static_cast<double>(n) * TIME_BIN_DUST_MODEL, 1.0, tau_infall) *
                       TIME_BIN_DUST_MODEL;
    }

    inline double StellarMass(double M_star_pre, double M_gas, double M_return, double schmidt_index,
                              double tau_SF, size_t n) noexcept
    {
        return M_star_pre + TIME_BIN_DUST_MODEL * (SFH::SFR(M_gas, schmidt_index, tau_SF, static_cast<double>(n) * TIME_BIN_DUST_MODEL) - M_return);
    }

    // SN によるダストの破壊を計算
    inline val2
    SNDestruction(double M_gas, double M_Z, double SN_rate, const val &a_val, const val2 &M_dust_val2,
                  const val4 &dest_eta_val4) noexcept
    {
        // ダストの破壊量を格納するための2次元配列を初期化
        auto M_SND_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
        
        // 超新星によるダストの破壊が無効な場合、空の配列を返す
        if (isSN_DEST == false)
            return M_SND_val2;

        // 金属量に基づいて破壊効率のインデックスを決定
        auto l = std::size_t(3);
        if (M_Z / M_gas / Zsun < 5.e-3)
            l = 0; // 低金属量
        else if (M_Z / M_gas / Zsun < 5.e-1)
            l = 1; // 中金属量
        else if (M_Z / M_gas / Zsun < 0.5)
            l = 2; // 高金属量

        // 各ダスト種に対してループ
        for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
        {
            // 各ダスト半径に対してループ
            for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
            {
                // ダストの破壊量を計算するための変数を初期化
                // dest_etaの配列の左から1番目はダスト種、2番目は、破壊前のダストサイズ、3番目は破壊後のダストサイズ、4番目はmetallicityをそれぞれ表している。
                // 単純に破壊によって減少している分を計算する際はpre=i_radiusの場合のみなので、for文で計算を回す必要がないはず(24/5/2013)
                auto dest = 0.0;
                for (auto pre = N_RADIUS_MIN; pre < N_MAX_DUST_RADIUS; ++pre)
                {
                    // 超新星衝撃によるダストの減少を計算
                    if (i_a != N_MAX_DUST_RADIUS - 1)
                    {
                        if (pre != N_MAX_DUST_RADIUS - 1)
                        {
                            // 破壊前のダストサイズに基づいて減少量を計算
                            dest += M_dust_val2[pre][i_s] * std::pow(a_val[i_a] / a_val[pre], 3.0) *
                                    dest_eta_val4[l][pre + 1][i_a + 1][i_s + 1];
                        }
                    }
                    else
                    {
                        // 最大ダスト半径の場合の処理
                        if (pre == N_MAX_DUST_RADIUS - 1)
                        {
                            dest += M_dust_val2[pre][i_s] * std::pow(a_val[i_a] / a_val[pre], 3.0);
                        }
                    }
                }
                // 超新星発生率に基づいてダストの破壊量を計算
                M_SND_val2[i_a][i_s] =
                    (-SN_rate * swept(M_gas, M_Z) / M_gas * (M_dust_val2[i_a][i_s] - dest)) *
                    TIME_BIN_DUST_MODEL;
            }
        }
        // 計算したダストの破壊量を返す
        return M_SND_val2;
    }

    auto MakeDustEvolutionClassVector(const val &dust_radius, const val2 &m_grain);

    val2 DustMassEvolution(double M_gas, const val &M_X_val, const val &a_val, const val &volume_val,
                           const val2 &M_dust_val2) noexcept;

    inline val2 DustMass(double M_gas, double schmidt_index, double tau_SF, const val2 &M_dust_val2,
                         const val2 &Y_d_val2, const val2 &M_SND_val2,
                         const val2 &M_evolution_val2, size_t n) noexcept
    {
        // ダスト質量の変化量を格納するための2次元配列を初期化
        auto dM_dust_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);

        // 各ダスト半径に対してループ
        for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
        {
            // 各ダスト種に対してループ
            for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
            {
                // Nishidaの論文の式(4.4)に基づいて、星化する質量 -D(t)SFR(t) を計算
                const auto M_astration = -isINJECTION * M_dust_val2[i_a][i_s] / M_gas *
                                         SFH::SFR(M_gas, schmidt_index, tau_SF, static_cast<double>(n) * TIME_BIN_DUST_MODEL);

                // ダストの質量を時間に基づいて計算
                const auto M_dust_by_star = TIME_BIN_DUST_MODEL * (M_astration + Y_d_val2[i_s][i_a]);

                // NaNチェック: ダストの質量がNaNの場合、警告を表示
                if (std::isnan(M_dust_by_star))
                    std::cout << "i_a = " << i_a << ", i_s = " << i_s << ", M_dust_by_star is nan!!"
                              << std::endl;
                if (std::isnan(M_evolution_val2[i_a][i_s]))
                    std::cout << "i_a = " << i_a << ", i_s = " << i_s << ", M_evolution is nan!!"
                              << std::endl;

                // ダストの質量の変化量を計算
                dM_dust_val2[i_a][i_s] =
                    M_dust_by_star + M_SND_val2[i_a][i_s] + M_evolution_val2[i_a][i_s];

                // ダストの質量が非常に小さい場合、ダストの質量を0に設定
                if (dM_dust_val2[i_a][i_s] + M_dust_val2[i_a][i_s] < DBL_MIN)
                {
                    dM_dust_val2[i_a][i_s] = -M_dust_val2[i_a][i_s]; // 総質量を0に設定
                }
            }
        }
        // ダストの質量の変化量を現在のダスト質量に加算して返す
        return dM_dust_val2 + M_dust_val2;
    }

} // namespace asano_model
#endif // ASANO_MODEL_FUNCTION_H
