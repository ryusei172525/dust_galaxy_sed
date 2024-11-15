#ifndef MKSED_MAKE_SED_H
#define MKSED_MAKE_SED_H

#include <vector>
#include <valarray>
#include <string>
#include <iostream>

#include <file_util.h>
#include <file_name.h>
#include <constants.h>
#include <dust_util.h>
#include <calculation_util.h>

namespace make_SED {
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using SEDFile = my_util::SED_util::SEDFileUtil;

template<class T>
void Print(const std::string& comment, const T& container) {
    std::cout << comment;
    for (const auto& item : container) {
        std::cout << ", " << item;
    }
    std::cout << std::endl;
}

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

[[nodiscard]] inline val ReadTransmissionRateFile(const FreeParameter& free_params) {
    auto ifs = my_util::FileUtil::IfsOpen(TransmissionRateFileName(free_params));
    SEDFile::SkipLine(ifs);

    auto transmission_rate_val = val(N_LAMBDA);
    for (auto&& transmission_rate : transmission_rate_val) {
        auto buff = 0.0;
        ifs >> buff >> transmission_rate >> buff;
    }
    return transmission_rate_val;
}

val ReadLuminosityOfDust(const std::string& dust_species, const FreeParameter& free_params) {
    auto       ifs  = SEDFile::IfsOpen(DustRadiationFileName(dust_species, free_params));
    const auto head = SEDFile::SkipLine(ifs);

    auto L_dust_val = val(N_LAMBDA);
    for (auto&& L_dust : L_dust_val) {
        auto buff = 0.0; // lambda
        ifs >> buff >> L_dust;
    }
    return L_dust_val;
}

void WriteSED(const FreeParameter& free_params, const val& lambda_cm_val,
              const val& transmission_rate_val, const val& L_star_val, const val2& L_dust_val2) {
    const auto L_star_attenuated_val = L_star_val * transmission_rate_val;

    const auto n_species = L_dust_val2.size();

    auto L_dust_attenuated_val2 = val2(val(N_LAMBDA), n_species);

    for (auto i = std::size_t(0); i < n_species; ++i) {
        L_dust_attenuated_val2[i] = L_dust_val2[i] * transmission_rate_val;
    }

    const auto SED_val2 = L_star_attenuated_val + L_dust_attenuated_val2[0]
                          + L_dust_attenuated_val2[1] + L_dust_attenuated_val2[2] +
                          L_dust_attenuated_val2[3];

    auto ofp = SEDFile::OfpOpen(SEDFileName(free_params), true);
    fprintf(ofp, "Wavelength(cm)");
    fprintf(ofp, " Total Stellar_int Stellar Sil Gra PAHion PAHneu\n");

    for (auto li = std::size_t(0); li < N_LAMBDA; ++li) {
        fprintf(ofp, "%le", lambda_cm_val[li]);
        fprintf(ofp, " %le %le %le", SED_val2[li], L_star_val[li],
                L_star_attenuated_val[li]);
        for (auto i = std::size_t(0); i < n_species; ++i) {
            fprintf(ofp, " %le", L_dust_attenuated_val2[i][li]);
        }
        fprintf(ofp, "\n");
    }
    std::fclose(ofp);
}

void CombineSEDs(const FreeParameter& free_params, const std::vector<std::string>& dust_species_vec,
                 const val& lambda_cm_val, const val& L_star_val) {
    const auto start_time = my_util::MyTime::GetTime();
// まずは星の放射と attenuation factor を読み込む
    const auto        transmission_rate_val = SEDFile::ReadTransmissionRateFile(free_params);
// Read dust radiation from file
    static const auto n_species             = dust_species_vec.size();
    auto              L_dust_val2           = val2(val(N_LAMBDA), n_species);
    for (auto         i                     = std::size_t(0); i < n_species; ++i)
        L_dust_val2[i] = make_SED::ReadLuminosityOfDust(dust_species_vec[i], free_params);
    make_SED::WriteSED(free_params, lambda_cm_val, transmission_rate_val, L_star_val, L_dust_val2);

    std::cout << "Combine SEDs: ";
    my_util::MyTime::PrintElapsedTime(start_time, "sec");

}

}
#endif //MKSED_MAKE_SED_H
