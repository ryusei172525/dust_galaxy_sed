#ifndef MKSED_MC_DUST_HIT_FILE_H
#define MKSED_MC_DUST_HIT_FILE_H

#include <valarray>
#include <random>

#include <SED_file_util.h>
#include <constants.h>
#include <file_name.h>

namespace mc_dust_hit {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using SEDFile = my_util::SED_util::SEDFileUtil;

std::vector<std::string> DustSpeciesVector(int argc, char** argv) {
    // コマンドライン引数が無いときは全てのダスト種を返す
    if (argc == 1) return std::vector<std::string>{"Sil", "Gra", "PAHion", "PAHneu"};

    auto dust_species_vector = std::vector<std::string>();
    dust_species_vector.reserve(static_cast<unsigned long>(argc - 1));
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


// すべての時間ですべての波長に関して判定を行う。
// 同時刻で衝突してしまうときは、エネルギーを足す
// 終了条件は合計のヒット数
//[[nodiscard]] inline auto
//MCSimulation(double dt, const val& dp_val, const val& E_photon_val) noexcept {
//    auto rnd = std::random_device();     // 非決定的な乱数生成器を生成
//    auto mt  = std::mt19937(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
//    auto mt = std::mt19937(0);     //  シードを追試可能なように固定しておく
//
//    auto hit_time_energy_map = std::map<double, double>(); //output
//
//     dp に対応したベルヌーイ分布を作成
//    auto bernoulli_vec = std::vector<std::bernoulli_distribution>(N_LAMBDA);
//
//    for (auto i = INDEX_LYMAN_BREAK; i < INDEX_LYMAN_BREAK + 1000; ++i) {
//        bernoulli_vec[i] = std::bernoulli_distribution(dp_val[i]);
//    }
//
//    auto t_elapsed = 0.0;
//
//    for (auto hit_count = std::size_t(0); hit_count < static_cast<std::size_t>(1e3);) {
//        t_elapsed += dt;
//         上限値をprobabilityから決めるほうがよい
//        for (auto li = INDEX_LYMAN_BREAK; li < INDEX_LYMAN_BREAK + 1000; ++li) {
//            if (bernoulli_vec[li](mt)) {
//                hit_time_energy_map[t_elapsed] += E_photon_val[li];
//                ++hit_count;
//            }
//        }
//    }
//    return hit_time_energy_map;
//}

[[nodiscard]] inline auto
MCSimulationAtOneRadius(double dt, double total_probability, const val& dpdt_val,
                        const val& E_photon_val) noexcept {
    const auto hit_max             = std::size_t(1e4);
    auto       hit_time_energy_val = val2(val(hit_max), 2); // return value

//    auto rnd = std::random_device();     // 非決定的な乱数生成器を生成
//    auto mt  = std::mt19937(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    auto mt = std::mt19937();     //  シードを追試可能なように固定しておく

    auto bernoulli_distribution = std::bernoulli_distribution(total_probability);
    auto discrete_distribution  = std::discrete_distribution<std::size_t>(std::begin(dpdt_val),
                                                                          std::end(dpdt_val));

    auto t_elapsed = 0.0;

    for (auto hit_count = std::size_t(0); hit_count < hit_max;) {
        t_elapsed += dt;
        if (bernoulli_distribution(mt)) {
            hit_time_energy_val[0][hit_count] = t_elapsed;
            hit_time_energy_val[1][hit_count] = E_photon_val[discrete_distribution(mt)];
            ++hit_count;
        }
    }
    return hit_time_energy_val;
}

[[nodiscard]] inline auto
MCSimulation(double total_probability, const val& E_photon_val, const val& dt_val,
             const val2& dpdt_val2) noexcept {
    [[maybe_unused]] const auto start_time = my_util::MyTime::GetTime();

    const auto n_a                  = dt_val.size();
    auto       hit_time_energy_val3 = val3(val2(n_a), 2);

    for (auto di = std::size_t(0); di < n_a; ++di) {
        const auto hit_time_energy_val2 = MCSimulationAtOneRadius(dt_val[di], total_probability,
                                                                  dpdt_val2[di], E_photon_val);
        hit_time_energy_val3[0][di] = hit_time_energy_val2[0];
        hit_time_energy_val3[1][di] = hit_time_energy_val2[1];
    }

//    std::cout << "MC simulation: ";
//    my_util::MyTime::PrintElapsedTime(start_time, "milli");

    return hit_time_energy_val3;
}

void WriteHitTimeEnergy(const std::string& dust_species, const FreeParameter& free_params,
                        const val3& hit_time_energy_val3) {
    auto ofp = my_util::FileUtil::OfpOpen(HitTimeFileName(dust_species, free_params));

    const auto n_a = hit_time_energy_val3[0].size();

    for (auto di = std::size_t(0); di < n_a; ++di) {
        const auto hit_count = hit_time_energy_val3[0][di].size();
        std::fprintf(ofp, "%lu\n", hit_count);
        for (auto i = std::size_t(0); i < hit_count; ++i) {
            std::fprintf(ofp, "%le %le\n", hit_time_energy_val3[0][di][i],
                         hit_time_energy_val3[1][di][i]);
        }
        std::fclose(ofp);
    }
}
}
#endif //MKSED_MC_DUST_HIT_FILE_H
