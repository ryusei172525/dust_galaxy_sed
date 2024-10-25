#ifndef DUST_TEMPERATURE_H
#define DUST_TEMPERATURE_H

#include <cfloat>

#include <SED_file_util.h>
#include <SED_setting_util.h>
#include <file_name.h>
#include <constants.h>

namespace dust_temperature {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using SEDFile = my_util::SED_util::SEDFileUtil;

std::vector<std::string> DustSpeciesVector(int argc, char** argv) noexcept {
    // コマンドライン引数が無いときは全てのダスト種を返す
    if (argc == 1) return {"Sil", "Gra", "PAHion", "PAHneu"};

    auto dust_species_vector = std::vector<std::string>();
    dust_species_vector.reserve(static_cast<std::size_t>(argc - 1));
    for (auto i = 1; i < argc; ++i) {
        // argv[i] と "Sil" などを比較するのに便利なように std::string に変換する
        auto dust_species = std::string(argv[i]);
        if (dust_species == "Sil" || dust_species == "Gra" || dust_species == "PAHion" ||
            dust_species == "PAHneu") {
            dust_species_vector.push_back(dust_species);
        } else {
            std::cout << "Error: \"" << dust_species
                      << "\" can not accepted dust species (accepted dust are only Sil, Gra, PAHion or PAHneu) !"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    return dust_species_vector;
}

std::size_t NDustRadius(const std::string& dust_species) noexcept {
    if (dust_species == "Sil" || dust_species == "Gra") return N_DUST_RADIUS;
    return N_DUST_RADIUS_PAH; // PAHion or  PAHneu
}

// hit counts is not constant
inline void ReadHitTimeFile(std::size_t n_a, const std::string& ifn, val2& t_hit_val2,
                            val2& E_hit_val2) {
    auto ifs = my_util::FileUtil::IfsOpen(ifn);

    t_hit_val2.resize(n_a);
    E_hit_val2.resize(n_a);

    for (auto di = std::size_t(0); di < n_a; ++di) {
        auto n_hit = std::size_t(0);
        ifs >> n_hit;
        t_hit_val2[di].resize(n_hit);
        E_hit_val2[di].resize(n_hit);
        for (auto i = std::size_t(0); i < n_hit; ++i) ifs >> t_hit_val2[di][i] >> E_hit_val2[di][i];
    }
}


[[nodiscard]] inline double
IntegrateEnergyDensity(const val& lambda_cm_val, const val& u_cgs_val) noexcept {
    // ファイルから波長とルミノシティを読み込む
    auto sum_u_cgs = 0.0;

    for (auto li = std::size_t(1); li < N_LAMBDA; ++li) {
        sum_u_cgs +=
        0.5 * (u_cgs_val[li - 1] + u_cgs_val[li]) * (lambda_cm_val[li] - lambda_cm_val[li - 1]);
    }
    return sum_u_cgs;
}

// ダストの平衡温度を求める. Takeuchi 2003 equ.8 より
[[nodiscard]] inline val EquilibriumTemperature(const std::string& dust_species, double sum_u_cgs,
                                                const val& a_val) noexcept {
    constexpr auto a = h_P * cl / (PI * k_B);
    auto           A = A_Sil; // A is constant value in Takeuchi 2003
    if (dust_species != "Sil") A = A_C;
    const auto b_val = 945 * sum_u_cgs / (960 * PI * (2 * PI * A * a_val) * h_P * cl);
    return a * std::pow(b_val, 1.0 / 6);
}

[[nodiscard]] inline val2
TemperatureHistogram(double sum_u_cgs, const val& a_cm_val, const val2& t_hit_val2,
                     const val2& E_hit_val2, const my_util::dust_util::Dust& dust) {
    [[maybe_unused]] const auto start_time = my_util::MyTime::GetTime();

    const auto T_equilibrium_val = EquilibriumTemperature(dust.GetDustSpecies(),
                                                          sum_u_cgs, a_cm_val);

    const auto T_hist_val2 = dust.CalculateTemperature(T_equilibrium_val, t_hit_val2, E_hit_val2);
//    std::cout << "Dust temperature: ";
//    my_util::MyTime::PrintElapsedTime(start_time, "milli");

    return T_hist_val2;
}

inline void
WriteHist(const std::string& dust_species, const FreeParameter& free_params, const val& a_cm_val,
          const val& T_val, const val2& T_hist_val2) {

    auto ofp = my_util::FileUtil::OfpOpen(HistogramFileName(dust_species, free_params));
    // Write header
    std::fprintf(ofp, "Temperature");
    for (const auto& a : a_cm_val) std::fprintf(ofp, " %.2ecm", a);
    std::fprintf(ofp, "\n");

    for (auto i_T = std::size_t(0); i_T < N_TEMPERATURE; ++i_T) {
        std::fprintf(ofp, "%e", T_val[i_T]);
        for (const auto& T_hist_val : T_hist_val2) {
            std::fprintf(ofp, " %le", T_hist_val[i_T]);
        }
        std::fprintf(ofp, "\n");
    }
    std::fclose(ofp);
}
}
#endif // DUST_TEMPERATURE_H
