#ifndef SED_MODEL_PROJECTS_SED_MODEL_H
#define SED_MODEL_PROJECTS_SED_MODEL_H

#include <SED_file_util.h>
#include <dust_util.h>
#include <calculation_util.h>
#include <stellar_spectrum.h>
#include <dust_temperature.h>
#include <radiative_transfer.h>
#include <mc_dust_hit.h>
#include <dust_radiation.h>
#include <combine_SED.h>

namespace SED_model {
using SEDFile = my_util::SED_util::SEDFileUtil;
using Time = my_util::MyTime;

using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using vec = std::vector<double>;

auto MakePairVector(const vec& vec1, const vec& vec2) {
    auto tau_pair_vec = std::vector<std::pair<double, double>>();
    for (const auto& value1 : vec1) {
        for (const auto& value2 : vec2) {
            tau_pair_vec.emplace_back(std::make_pair(value1, value2));
        }
    }
    return tau_pair_vec;
}

inline void AgeDependent(const std::vector<std::string>& dust_species_vec, double M_gal,
                         const val& lambda_cm_val, const val& a_cm_val, const val& E_photon_val,
                         double D, const val& k_val, const val& omega_val, const val& g_val,
                         std::vector<my_util::dust_util::Dust>& dust_vec,
                         FreeParameter& free_params,
                         const std::vector<std::pair<double, double>>& geometry_pair_vec,
                         const val& L_star_val, const val& f_young_val, const val& Cabs_sum_val, double n0_cnm) {

    // first: h_dust, second: R_gal
    for (const auto& geometry_pair : geometry_pair_vec) {
        free_params.SetDustScaleHeight(geometry_pair.first);
        free_params.SetGalaxyRadius(geometry_pair.second);

        /////////////////////// Radiative transfer ////////////////////////////////////////////////////
        const auto ai = free_params.ai_;

        auto RT = radiative_transfer::RadiativeTransfer(D, free_params.dz_, f_young_val,
                                                        k_val, omega_val, g_val, n0_cnm);

        // RT.WriteMGAParameter(free_params);
        const auto transmission_rate_val2 = RT.Calculate(free_params);
        const auto T_average_val         = RT.WriteTransmissionRate(transmission_rate_val2,
                                                                    free_params);

        auto transmission_rate_val = val(N_LAMBDA);
        for (auto i = std::size_t(0); i < N_LAMBDA; ++i) {
            transmission_rate_val[i] =transmission_rate_val2[i][0];
        }

        const auto u_cgs_val = SEDFile::CalcEnergyDensity(L_star_val, Cabs_sum_val,
                                                          M_gal, transmission_rate_val);
        const auto sum_u_cgs = dust_temperature::IntegrateEnergyDensity(lambda_cm_val, u_cgs_val);

        auto       dpdt_map       = std::map<std::string, val2>();
        auto       dt_map         = std::map<std::string, val>();
        const auto probability_MC = 1e-2;

        for (const auto& dust : dust_vec) {
            const auto dust_species = dust.GetDustSpecies();
            dpdt_map.emplace(dust_species, dust.DpDt(lambda_cm_val, E_photon_val, u_cgs_val));
            dt_map.emplace(dust_species,
                           dust.DeltaTime(probability_MC, dpdt_map[dust_species]));
        }


        /////// MC dust hit //////////////////////////////////////////////////////////////////////////
        const auto start_time_dust_species_loop = Time::GetTime();

#pragma omp parallel for
        for (const auto& dust : dust_vec) {
            const auto dust_species         = dust.GetDustSpecies();
            const auto hit_time_energy_val3 = mc_dust_hit::MCSimulation(probability_MC,
                                                                        E_photon_val,
                                                                        dt_map[dust_species],
                                                                        dpdt_map[dust_species]);
//                mc_dust_hit::WriteHitTimeEnergy(dust_species, free_params, hit_time_energy_val3);

            // Dust temperature
            const auto T_hist_val2 = dust_temperature::TemperatureHistogram(sum_u_cgs,
                                                                            a_cm_val,
                                                                            hit_time_energy_val3[0],
                                                                            hit_time_energy_val3[1],
                                                                            dust);

            const auto T_val = dust.GetTemperatureBin();
            dust_temperature::WriteHist(dust_species, free_params, a_cm_val, T_val, T_hist_val2);

            // Dust radiation
            static const auto BB_val2 = dust_radiation::BlackBodyVarray2(T_val, lambda_cm_val);

            const auto L_dust_val = dust.DustRadiation(ai, T_hist_val2, BB_val2);


            dust.WriteDustRadiationPerRadius(dust_species, ai, T_hist_val2, BB_val2, lambda_cm_val,
                                             free_params);

            dust_radiation::WriteLuminosity(dust_species, free_params, lambda_cm_val,
                                            L_dust_val);
        }

        my_util::MyTime::PrintElapsedTime(start_time_dust_species_loop, "sec",
                                          "Dust species loop: ");

        make_SED::CombineSEDs(free_params, dust_species_vec, lambda_cm_val, L_star_val);
        std::cout << std::endl;
    }
}

} // namespace SED_model
#endif //SED_MODEL_PROJECTS_SED_MODEL_H
