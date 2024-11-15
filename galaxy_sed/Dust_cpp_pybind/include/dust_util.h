#ifndef SED_MODEL_DUST_UTIL_H
#define SED_MODEL_DUST_UTIL_H

#include <SED_file_util.h>
#include <SED_setting_util.h>

#include <utility>

namespace my_util::dust_util {
using SEDFile = SED_util::SEDFileUtil;

class Dust {
  private:
    std::string dust_species_;
    std::size_t n_a_;
    double      T_sublimation_;
    val         a_cm_val_;
    val         n_atom_val_;
    val         T_val_;
    val2        E_inter_val2_;
    val2        Qabs_val2_;
    val2        n_val2_;
    val2        t_cool_val2_;
    val2        t_int_cool_val2_;
    val2        epsilon_val2_;

    inline std::size_t NumberOfDustRadius() noexcept {
        auto n_dust_radius = N_DUST_RADIUS;
        if (dust_species_ == "PAHion" || dust_species_ == "PAHneu")
            n_dust_radius = N_DUST_RADIUS_PAH;
        return n_dust_radius;
    }

    [[nodiscard]] val NumberOfAtom() const noexcept {
        auto C = C_Sil;
        if (dust_species_ != "Sil") C = C_Gra;
        return C * a_cm_val_ * a_cm_val_ * a_cm_val_;
    }

    [[nodiscard]] inline double
    SublimationTemperature() const noexcept {
        if (dust_species_ == "Sil") return SILICATE_SUBLIMATION_TEMPERATURE;
        return CARBON_SUBLIMATION_TEMPERATURE;
    }

    [[nodiscard]] static inline double HCRatio(double N_atom) noexcept {
        if (N_atom <= 25) return 0.5;
        if (N_atom <= 100) return 0.5 / std::sqrt(N_atom / 25);
        return 0.25;
    }

    inline val2 ReadInternalEnergy(const std::string& dust_species) {
        auto ifs = my_util::FileUtil::IfsOpen(EnthalpyPerAtomFileName(dust_species));

        T_val_.resize(N_TEMPERATURE);

        auto      E_inter_per_atom_val = val(N_TEMPERATURE);
        for (auto i_T                  = std::size_t(0); i_T < N_TEMPERATURE; ++i_T) {
            ifs >> T_val_[i_T] >> E_inter_per_atom_val[i_T];
        }

        auto      E_inter_val2 = val2(n_a_);
        for (auto di           = std::size_t(0); di < n_a_; ++di) {
            E_inter_val2[di] = E_inter_per_atom_val * (n_atom_val_[di] - 2);
        }
        return E_inter_val2;
    }

    [[nodiscard]] inline val2 InternalEnergy() {
        // Silicate or Graphite case
        if (dust_species_ == "Sil" || dust_species_ == "Gra")
            return ReadInternalEnergy(dust_species_);

        // PAHs case
        const auto E_inter_gra_val2 = ReadInternalEnergy("Gra");
        const auto E_inter_pah_val2 = ReadInternalEnergy(dust_species_);
        auto       E_inter_val2     = val2(n_a_);

        for (auto di = std::size_t(0); di < n_a_; ++di) {
            const auto H_C_ratio = HCRatio(n_atom_val_[di]);
            E_inter_val2[di] = E_inter_gra_val2[di] + H_C_ratio * E_inter_pah_val2[di];
        }
        return E_inter_val2;
    }


    inline void ReadCoolingTimeFile() {
        t_cool_val2_.resize(n_a_);
        t_int_cool_val2_.resize(n_a_);
        epsilon_val2_.resize(n_a_);
        for (auto di = std::size_t(0); di < n_a_; ++di) {
            // file 名として使うために radius を um 単位の小数点以下３桁の指数表示に変換する
            auto dust_radius_cm_oss = std::ostringstream();
            dust_radius_cm_oss << std::scientific << std::setprecision(3) << a_cm_val_[di];
            auto ifs = my_util::FileUtil::IfsOpen(
            CoolingTimeFileName(dust_species_, dust_radius_cm_oss.str()));

            t_cool_val2_[di].resize(N_TEMPERATURE);
            t_int_cool_val2_[di].resize(N_TEMPERATURE);
            epsilon_val2_[di].resize(N_TEMPERATURE);
            SEDFile::SkipLine(ifs); // Skip header
            for (auto T = std::size_t(0); T < N_TEMPERATURE; ++T) {
                auto buff = 0.0; // temperature
                ifs >> buff >> t_cool_val2_[di][T]
                    >> t_int_cool_val2_[di][N_TEMPERATURE - T] // 昇順になるように逆から挿入する
                    >> epsilon_val2_[di][T];
            }
        }
    }

    [[nodiscard]] inline std::size_t SearchTemperatureIterationNumber(double T) const noexcept {
        // T 以上の最初の温度ビンのイテレータ
        const auto iter_T = std::lower_bound(std::begin(T_val_), std::end(T_val_), T);
        return static_cast<std::size_t>(std::distance(std::begin(T_val_), iter_T));
    }

    [[nodiscard]] static inline std::size_t
    Heating(double hit_energy_erg, const val& E_internal_val,
            double& E_internal_now) noexcept {
        E_internal_now += hit_energy_erg;
        // key より大きい最初の要素のイテレータを返す
        // key が E_internal_val の範囲外だった場合は範囲外のイテレータを戻す
        const auto heat_after_iter = std::upper_bound(std::begin(E_internal_val),
                                                      std::end(E_internal_val), E_internal_now);

        if (heat_after_iter == std::end(E_internal_val)) return N_TEMPERATURE - 1;

        // heat_after_iter が範囲よりも大きかったときの処理をかんたんにするため、-1はここでつける
        return static_cast<std::size_t>(std::distance(std::begin(E_internal_val),
                                                      heat_after_iter - 1));
    }

    [[nodiscard]] static inline std::size_t
    Cooling(std::size_t i_T, double t_delta, const val& t_integral_cool_val) noexcept {
        if (t_delta < 0) return i_T;
        const auto t_after = t_integral_cool_val[N_TEMPERATURE - i_T] + t_delta;

        // lower_bound は t_after を超える初めのイテレータを返すがその前がほしいので -1 する
        const auto iter_after = std::lower_bound(std::begin(t_integral_cool_val),
                                                 std::end(t_integral_cool_val), t_after) - 1;

        const auto i_T_after = std::distance(std::begin(t_integral_cool_val), iter_after);

        // t_integral_cool_val は昇順にするために温度とビンが逆なので N_TEMPERATURE から引いてやることで
        // 温度のビン番号になる
        return N_TEMPERATURE - static_cast<std::size_t>(i_T_after);
    }

    static inline void
    FillTemperatureHistogram(std::size_t i_T_cool_before, std::size_t i_T_cool_after,
                             double t_delta, const val& t_cool_val, const val& E_inter_val,
                             const val& epsilon_val, val& T_hist_val,
                             double& E_inter_now) noexcept {
        // 冷却によって1bin以上温度が変わる前に次の光子があたった場合
        if (i_T_cool_after == i_T_cool_before) {
            T_hist_val[i_T_cool_after] += t_delta;
            return;
        }

        // E_inter_now が i_T_cool_before での内部エネルギーを超えている分の補正
        auto t_elapsed =
             (E_inter_now - E_inter_val[i_T_cool_before]) / epsilon_val[i_T_cool_before];
        T_hist_val[i_T_cool_before + 1] += t_elapsed;

        for (auto i_T = i_T_cool_before; i_T > i_T_cool_after; --i_T) {
            T_hist_val[i_T] += t_cool_val[i_T];
            t_elapsed += t_cool_val[i_T];
            E_inter_now -= t_cool_val[i_T] * epsilon_val[i_T];
        }

        if (t_delta - t_elapsed < 0)
            std::cout << "t_delta - t_elapsed is minus!! t_delta = " <<
                      t_delta << ", t_elapsed = " << t_elapsed
                      << ", i_T_cool_before = " << i_T_cool_before << ", i_T_cool_after = "
                      << i_T_cool_after << std::endl;

        // 余った時間を冷却後の温度に足し合わせる
        T_hist_val[i_T_cool_after] += t_delta - t_elapsed;
        E_inter_now -= (t_delta - t_elapsed) * epsilon_val[i_T_cool_after];
        if (E_inter_now < 0) E_inter_now = 0;
    }


    [[nodiscard]] inline val
    CoolingAndHeating(double T_equilibrium, const val& t_hit_val, const val& E_hit_val,
                      const val& t_integral_cool_val, const val& t_cool_val, const val& E_inter_val,
                      const val& epsilon_val) const noexcept {
        // T_equilibrium 以上の最初の温度ビンのイテレーション番号
        auto i_T = SearchTemperatureIterationNumber(T_equilibrium);

//        if (t_hit_val[0] > t_cool_val[i_T] * 100) {
//            auto T_hist_val = val(N_TEMPERATURE);
//            T_hist_val[i_T] = 1.0;
//            return T_hist_val;
//        }

        auto E_inter_now = E_inter_val[i_T]; // 初期状態のエンタルピー

        auto       t_before   = 0.0; // 前回の衝突時刻
        auto       T_hist_val = val(N_TEMPERATURE); // temperature distribution
        auto       n_lost     = 0;
        const auto n_hit      = t_hit_val.size();

        auto total_time = 0.0;

        for (auto i = std::size_t(0); i < n_hit; ++i) {
            i_T = Heating(E_hit_val[i], E_inter_val, E_inter_now);
            // 昇華温度に達していた場合を判断する
            //if (i_T == N_TEMPERATURE - 1) {
            if (T_val_[i_T] > T_sublimation_) {
                i_T         = SearchTemperatureIterationNumber(T_equilibrium);
                t_before    = t_hit_val[i]; // 次の時間差を計算するためにここで記録
                E_inter_now = E_inter_val[i_T];
                ++n_lost;
                continue;
            }

            const auto i_T_before_cool = i_T; // FillTemperatureHistogram で冷却前の温度を使うので記録しておく

            // この時間で冷却が行われる
            const auto t_delta = (t_hit_val[i] - t_before);

            total_time += t_delta;

            // t_delta には E_inter_now が T での内部エネルギーを超えている分の補正を入れる
            i_T = Cooling(i_T, t_delta - (E_inter_now - E_inter_val[i_T]) / epsilon_val[i_T],
                          t_integral_cool_val);

//            std::cout << "T_heat = " << T_val_[i_T_before_cool] << ", T_cool = " << T_val_[i_T]
//                      << ", t_delta = " << t_delta << ", EquilibriumTemperature = " << T_equilibrium
//                      << std::endl;

            // 温度ごとの滞在時間をヒストグラムに詰める
            FillTemperatureHistogram(i_T_before_cool, i_T, t_delta, t_cool_val, E_inter_val,
                                     epsilon_val, T_hist_val, E_inter_now);

            t_before = t_hit_val[i]; // 次の時間差を計算するためにここで記録
        }

        if (n_lost != 0) {
            const auto ratio = 1 - static_cast<double>(n_lost) / static_cast<double>(n_hit);
            //std::cout << "n_lost = " << n_lost << ", ratio = " << ratio << std::endl;
            return T_hist_val / T_hist_val.sum() * ratio; // 最後に規格化してから戻す
        }
        return T_hist_val / T_hist_val.sum(); // 最後に規格化してから戻す
    }

    // erg/s/cm
    // Nishida et al 2022の式(86)に対応
    [[nodiscard]] inline double
    LuminosityOfOneLambda(const val& Qabs_val, const val& n_val, const val& BB_val,
                          const val2& T_hist_val2) const noexcept {
        auto              L        = 0.0;
        static const auto S_pi_val = 4 * PI * PI * a_cm_val_ * a_cm_val_;
        const auto        tmp_val  = S_pi_val * Qabs_val * n_val;
        // forループはn_a_ 回繰り返され、各繰り返しで"L"に新しい項が加算されます。
        for (auto di = std::size_t(0); di < n_a_; ++di)
            L += tmp_val[di] * (BB_val * T_hist_val2[di]).sum();
        return L;
    }

    // erg/s/cm
    inline void
    LuminosityOfOneLambdaPerRadius(const val& Qabs_val, const val& n_val, const val& BB_val,
                                   const val2& T_hist_val2, std::ofstream& ofs) const noexcept {
        static const auto S_pi_val = 4 * PI * PI * a_cm_val_ * a_cm_val_;
        const auto        tmp_val  = S_pi_val * Qabs_val * n_val;

        for (auto di = std::size_t(0); di < n_a_; ++di)
            ofs << tmp_val[di] * (BB_val * T_hist_val2[di]).sum() << " ";
        ofs << std::endl;
    }


  public:
    Dust() = delete;

    inline Dust(std::string dust_species, val a_cm_val) :
    dust_species_(std::move(dust_species)),
    n_a_(NumberOfDustRadius()),
    T_sublimation_(SublimationTemperature()),
    a_cm_val_(std::move(a_cm_val)),
    n_atom_val_(NumberOfAtom()),
    E_inter_val2_(InternalEnergy()),
    Qabs_val2_(SEDFile::ReadQabsFile(n_a_, DustPropertyFileName(dust_species_, "Qabs"))) {
        auto t_cool_val2            = val2();
        auto t_cool_integrated_val2 = val2();
        auto epsilon_val2           = val2();
        ReadCoolingTimeFile();
    }

    ~Dust() = default;

    [[nodiscard]] inline auto GetDustSpecies() const noexcept { return dust_species_; }

    [[nodiscard]] inline auto GetNumberOfDustRadius() const noexcept { return n_a_; }

    [[nodiscard]] inline auto GetRadiusCmValarray() const noexcept { return a_cm_val_; }

    [[nodiscard]] inline auto GetTemperatureBin() const noexcept { return T_val_; }

    [[nodiscard]] inline auto GetNumberDistribution() const noexcept { return n_val2_; }

    inline void SetNumberDistribution(double M_gal, const val2& n_val2) noexcept {
        n_val2_.resize(n_val2.size());
        for (auto i = std::size_t(0); i < n_val2.size(); ++i) {
            n_val2_[i] = n_val2[i] * M_gal;
        }
    }

    [[nodiscard]] inline val2
    DpDt(const val& lambda_cm_val, const val& E_photon_val, const val& u_cgs_val) const noexcept {
        auto dpdt_val2 = val2(val(N_LAMBDA), n_a_);

        for (auto di = std::size_t(0); di < n_a_; ++di) {
            for (auto li = std::size_t(1); li < N_LAMBDA; ++li) {
                const auto dE = E_photon_val[li - 1] - E_photon_val[li];
                dpdt_val2[di][li] = PI * a_cm_val_[di] * a_cm_val_[di] * Qabs_val2_[li][di]
                                    * u_cgs_val[li] * lambda_cm_val[li] *
                                    lambda_cm_val[li] * lambda_cm_val[li]
                                    / (h_P * h_P * cl) * dE;
            }
        }
        return dpdt_val2;
    }

// 最大の衝突確率が max_probability になるように dt を決定する
    [[nodiscard]] inline val
    DeltaTime(double total_probability, const val2& dpdt_val2) const noexcept {
        auto      dt_val = val(n_a_);
        for (auto di     = std::size_t(0); di < n_a_; ++di)
            dt_val[di] = total_probability / dpdt_val2[di].sum();
        return dt_val;
    }

    [[nodiscard]] inline auto
    Dp(std::size_t n_age, const val2& dt_val2, const val3& dpdt_val3) const noexcept {
        // dpdt に十分小さい時間 dt をかけて、その間に光子が衝突する確率にする
        auto      dp_val3 = val3(val2(val(N_LAMBDA), n_a_), n_age);
        for (auto ai      = std::size_t(0); ai < n_age; ++ai) {
            for (auto di = std::size_t(0); di < n_a_; ++di) {
                dp_val3[ai][di] = dt_val2[ai][di] * dpdt_val3[ai][di];
            }
        }
        return dp_val3;
    }

    [[nodiscard]] inline val2
    CalculateTemperature(const val& T_equilibrium_val, const val2& t_hit_val2,
                         const val2& E_hit_val2) const noexcept {
        auto T_hist_val2 = val2(val(N_TEMPERATURE), n_a_);

        for (auto di = std::size_t(0); di < n_a_; ++di) {
            T_hist_val2[di] = CoolingAndHeating(T_equilibrium_val[di],
                                                t_hit_val2[di], E_hit_val2[di],
                                                t_int_cool_val2_[di],
                                                t_cool_val2_[di],
                                                E_inter_val2_[di], epsilon_val2_[di]);
        }
        return T_hist_val2;
    }

    [[nodiscard]] inline val
    DustRadiation(std::size_t ai, const val2& T_hist_val2, const val2& BB_val2) const noexcept {
        auto      L_val = val(N_LAMBDA);
        for (auto li    = std::size_t(0); li < N_LAMBDA; ++li) {
            L_val[li] = LuminosityOfOneLambda(Qabs_val2_[li], n_val2_[ai], BB_val2[li],
                                              T_hist_val2);
        }
        return L_val;
    }

    inline void WriteDustRadiationPerRadius(const std::string& dust_species, std::size_t ai,
                                            const val2& T_hist_val2, const val2& BB_val2,
                                            const val& lambda_cm_val,
                                            const FreeParameter& free_params) const noexcept {
        auto ofs = SEDFile::OfsOpen(DustRadiationPerRadiusFileName(dust_species, free_params));
        ofs << "Wavelength[cm] ";
        for (auto i = std::size_t(0); i < n_a_; ++i) {
            ofs << a_cm_val_[i] << " ";
        }
        ofs << std::endl;

        for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
            ofs << lambda_cm_val[li] << " ";
            LuminosityOfOneLambdaPerRadius(Qabs_val2_[li], n_val2_[ai],
                                           BB_val2[li], T_hist_val2, ofs);
        }
    }

};

/**
 * Make dust classes vector
 * @param dust_species_vec {Sil, Gra, PAHion, PAHneu}
 * @param a_cm_val dust radius [cm]
 * @return std::vector<Dust>
 */
inline std::vector<Dust>
MakeDustClassVector(const std::vector<std::string>& dust_species_vec, const val& a_cm_val) {
    auto dust_vec = std::vector<my_util::dust_util::Dust>();
    for (const auto& dust_species : dust_species_vec)
        dust_vec.emplace_back(my_util::dust_util::Dust(dust_species, a_cm_val));
    return dust_vec;
}

// 梨本 p33, Li & Draine 2001 p3
[[nodiscard]]  inline double GraphiteFraction(double radius_cm) noexcept {
    constexpr auto a_threshold = 50e-8; // 50A
    if (radius_cm < a_threshold) return 0.01;
    //auto f_gra = 0.01 + 0.99 * (1 - std::pow(a_threshold / radius_cm, 3));
    //std::cout << "radius = " << radius_cm
    //        << ", GraphiteFraction = " << f_gra << std::endl;
    //f_gra = 1 - (1 - 0.01) * std::min(1.0, std::pow(a_threshold / radius_cm, 3));
    //std::cout << "radius = " << radius_cm
    //        << ", GraphiteFraction = " << f_gra << std::endl;
    //return 1;
    return 1 - (1 - 0.01) * std::min(1.0, std::pow(a_threshold / radius_cm, 3));
}

[[nodiscard]]  inline val GraphiteFractionVarray(const val& a_cm_val) noexcept {
    auto f_gra_val = val(N_DUST_RADIUS);

    for (auto i = std::size_t(0); i < N_DUST_RADIUS_PAH; ++i)
        f_gra_val[i] = GraphiteFraction(a_cm_val[i]);

    for (auto i = N_DUST_RADIUS_PAH; i < N_DUST_RADIUS; ++i)
        f_gra_val[i] = 1.0;
    return f_gra_val;
}

[[nodiscard]]  inline double PAHneuFraction(double a_cm) noexcept {
    const auto radius_um = a_cm * 1e4; // 式が um で書いてあるので変換する
    return (1.0 / (1.0 + 3.5e6 * 0.460 * std::pow(radius_um, 3) * 1.14
                         * std::pow(100, 0.5) / 0.045)
            + 1.0 / (1.0 + 3.5e6 * 0.460 * std::pow(radius_um, 3) * 1.14
                           * std::pow(6000, 0.5) / 0.04))
           * 0.5;
}

[[nodiscard]] inline val PAHneuFractionVarray(const val& a_cm_val) noexcept {
    auto      f_PAHion_val = val(N_DUST_RADIUS_PAH);
    for (auto i            = std::size_t(0); i < N_DUST_RADIUS_PAH; ++i)
        f_PAHion_val[i] = PAHneuFraction(a_cm_val[i]);
    return f_PAHion_val;
}

/**
 * Read dust number distribution file and divide carbonaceous grains into graphite, neutral PAH, and ionized PAH
 * @param ifn
 * @param n_age
 * @param a_cm_val
 * @return val3 [dust species][galaxy age][dust radius]
 * 0: silicate, 1: graphite, 2: neutral PAH, 3: ionized PAH
 */
inline val3 ReadAndDivideDustNumberDistribution(const std::string& ifn, std::size_t n_age,
                                                const val& a_cm_val) {

    // n_C_sil_val3[0]: Carbon dust number, n_C_sil_val3[1]: Silicate dust number
    const auto n_C_sil_val3 = SEDFile::ReadSilicateAndCarbonDustDistributionFile(ifn, n_age);

    const auto n_sil_val2 = n_C_sil_val3[0];

    // Carbonaceous grain case
    const auto f_gra_val = GraphiteFractionVarray(a_cm_val);

    // TODO: Check ion or neu
    const auto f_PAHneu_val = PAHneuFractionVarray(a_cm_val);

    auto n_gra_val2    = val2(n_age);
    auto n_PAHneu_val2 = val2(n_age);
    auto n_PAHion_val2 = val2(n_age);

    for (auto i = std::size_t(0); i < n_age; ++i) {
        n_gra_val2[i]    = n_C_sil_val3[1][i] * f_gra_val;
        n_PAHneu_val2[i] = n_C_sil_val3[1][i] * (1 - f_gra_val) * f_PAHneu_val;
        n_PAHion_val2[i] = n_C_sil_val3[1][i] * (1 - f_gra_val) * (1 - f_PAHneu_val);
    }
    return {n_sil_val2, n_gra_val2, n_PAHneu_val2, n_PAHion_val2};
}
} // my_util
#endif //SED_MODEL_DUST_UTIL_H
