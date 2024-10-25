#ifndef MKSED_DUST_RADIATION_H
#define MKSED_DUST_RADIATION_H

#include <vector>
#include <string>
#include <valarray>
#include <iostream>

#include <file_name.h>
#include <file_util.h>
#include <constants.h>
#include <SED_file_util.h>

namespace dust_radiation {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using SEDFile = my_util::SED_util::SEDFileUtil;

std::vector<std::string> DustSpeciesVector(int argc, char** argv) {
    if (argc == 1) return std::vector<std::string>{"Sil", "Gra", "PAHion", "PAHneu"};

    auto dust_species_vector = std::vector<std::string>();

    for (auto i = 1; i < argc; ++i) {
        const auto dust_species = argv[i];
        if (dust_species == std::string("Sil") || dust_species == std::string("Gra") ||
            dust_species == std::string("PAHion") || dust_species == std::string("PAHneu")) {
            dust_species_vector.emplace_back(argv[i]);
        } else {
            std::cout << "Error: \"" << argv[i]
                      << "\" can not accepted dust species (accepted dust are only Sil, Gra, PAHion or PAHneu) !"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    return dust_species_vector;
}

// 梨本 p33, Li & Draine 2001 p3
double GraphiteFraction(double radius_cm) noexcept {
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

[[nodiscard]] inline val GraphiteFractionVarray(const val& a_cm_val) noexcept {
    auto f_gra_val = val(N_DUST_RADIUS);

    for (auto i = std::size_t(0); i < N_DUST_RADIUS_PAH; ++i)
        f_gra_val[i] = GraphiteFraction(a_cm_val[i]);

    for (auto i = N_DUST_RADIUS_PAH; i < N_DUST_RADIUS; ++i)
        f_gra_val[i] = 1.0;
    return f_gra_val;
}

[[nodiscard]] inline double PAHionFraction(double a_cm) noexcept {
    const auto radius_um = a_cm * 1e4; // 式が um で書いてあるので変換する
    return (1.0 / (1.0 + 3.5e6 * 0.460 * std::pow(radius_um, 3) * 1.14
                         * std::pow(100, 0.5) / 0.045)
            + 1.0 / (1.0 + 3.5e6 * 0.460 * std::pow(radius_um, 3) * 1.14
                           * std::pow(6000, 0.5) / 0.04))
           * 0.5;
}

[[nodiscard]] inline val PAHionFractionVarray(const val& a_cm_val) noexcept {
    auto      f_PAHion_val = val(N_DUST_RADIUS_PAH);
    for (auto i            = std::size_t(0); i < N_DUST_RADIUS_PAH; ++i)
        f_PAHion_val[i] = PAHionFraction(a_cm_val[i]);
    return f_PAHion_val;
}

// Calculate fraction of carbonaceous dust grains
[[nodiscard]] inline val2
DustNumberDistribution(const std::string& dust_species, std::size_t n_age, std::size_t n_a,
                       const val& a_cm_val, const val2& n_car_val2,
                       const val2& n_sil_val2) noexcept {
    if (dust_species == "Sil") return n_sil_val2; // sil の場合はそのまま戻す

    const auto f_gra_val = GraphiteFractionVarray(a_cm_val);
    auto       n_val2    = val2(val(n_a), n_age);

    if (dust_species == "Gra") {
        for (auto i = std::size_t(0); i < n_age; ++i) {
            n_val2[i] = n_car_val2[i] * f_gra_val;
        }
        return n_val2;
    }

    const auto f_PAHion_val = PAHionFractionVarray(a_cm_val);
    if (dust_species == "PAHion") {
        for (auto i = std::size_t(0); i < n_age; ++i)
            n_val2[i] = n_car_val2[i] * (1 - f_gra_val) * f_PAHion_val;
        return n_val2;
    }

    // PAHneu
    for (auto i = std::size_t(0); i < n_age; ++i)
        n_val2[i] = n_car_val2[i] * (1 - f_gra_val) * (1 - f_PAHion_val);
    return n_val2;
}

[[nodiscard]] inline double BlackBody(double T, double lambda_cm) noexcept {
    const auto in_exp = h_P * cl / (lambda_cm * k_B * static_cast<double>(T));
    return 2 * h_P * cl * cl / std::pow(lambda_cm, 5) / (std::exp(in_exp) - 1);
}

[[nodiscard]] inline val2 BlackBodyVarray2(const val& T_val, const val& lambda_cm_val) noexcept {
    auto black_body_val2 = val2(val(N_TEMPERATURE), N_LAMBDA);

    for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
        for (auto Ti = std::size_t(0); Ti < N_TEMPERATURE; ++Ti) {
            black_body_val2[li][Ti] = BlackBody(T_val[Ti], lambda_cm_val[li]);
        }
    }
    return black_body_val2;
}

val2 ReadTemperatureHistogram(const std::string& ifn, std::size_t n_a) {
    auto ifs = my_util::FileUtil::IfsOpen(ifn);

    auto T_hist_val2 = val2(val(N_TEMPERATURE), n_a);

    SEDFile::SkipLine(ifs); // skip header

    // 各温度に対するデータを読み取り、配列に格納する
    for (auto Ti = std::size_t(0); Ti < N_TEMPERATURE; ++Ti) {
        auto T = 0.0;
        ifs >> T; // ファイルから温度の値を読み取る
        std::cout << "T = " << T << std::endl;
        for (auto di = std::size_t(0); di < n_a; ++di) {
            ifs >> T_hist_val2[di][Ti]; // ファイルから読み取った値を二次元配列 T_hist_val2 に格納する
        }
    }

    // 二次元の配列を返す
    return T_hist_val2;
}


void WriteLuminosity(const std::string& dust_species, const FreeParameter& free_params,
                     const val& lambda_cm_val, const val& L_val) {
    auto ofp = SEDFile::OfpOpen(DustRadiationFileName(dust_species, free_params));
    fprintf(ofp, "wavelength(cm) Luminosity(erg/s/cm)\n"); // Write header
    for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
        fprintf(ofp, "%e %e\n", lambda_cm_val[li], L_val[li]);
    }
    std::fclose(ofp);
}
}
#endif //MKSED_DUST_RADIATION_H
