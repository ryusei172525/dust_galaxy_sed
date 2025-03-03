//
// Created by 西田和樹 on 2020/09/11.
//

#include "function.h"

#include <cfloat>
#include <vector>

#include "dust_evolution.h"

namespace asano_model
{
    auto MakeDustEvolutionClassVector(const val &dust_radius, const val2 &m_grain)
    {
        const auto ISM_name_vector = std::vector<std::string>{"CNM", "MC", "WNM"};

        auto evolution_vector = std::vector<DustEvolution>();
        for (const auto &ISM_name : ISM_name_vector)
        {
            evolution_vector.emplace_back(DustEvolution(ISMParams(ISM_name), dust_radius, m_grain));
        }
        return evolution_vector;
    }

    val2 DustMassEvolution(double M_gas, const val &M_X_val, const val &a_val, const val &volume_val,
                           const val2 &M_dust_val2) noexcept
    {
        // The timestep is required CFL condition (dt < da/dot{a}) (23/8/2011)
        static const auto m_grain_val2 = GrainMass(volume_val);
        static const auto m_grain_dm_val2 = GrainMassDm(volume_val);

        static auto evolution_vector = MakeDustEvolutionClassVector(a_val, m_grain_val2);

        // auto M_evolution_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);

        //    auto dM_in_MC = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
        //    auto dM_in_WNM = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
        //    auto M_evo3 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
        //    auto M_evo4 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS);
        auto dM_in_CNM = val2();
        auto dM_in_MC = val2();
        auto dM_in_WNM = val2();
        #pragma omp parallel sections
        {
            #pragma omp section
                        dM_in_CNM = evolution_vector[0].Calculation(M_gas, M_X_val, volume_val, m_grain_dm_val2,
                                                                    M_dust_val2);
            #pragma omp section
                        dM_in_MC = evolution_vector[1].Calculation(M_gas, M_X_val, volume_val, m_grain_dm_val2,
                                                                M_dust_val2);
            #pragma omp section
                        dM_in_WNM = evolution_vector[2].Calculation(M_gas, M_X_val, volume_val, m_grain_dm_val2,
                                                                M_dust_val2);
        }
        //    for (auto&& evolution : evolution_vector) {
        //        M_evolution_val2 += evolution.Calculation(M_gas, M_X_val, volume_val,//                                                  m_grain_dm_val2, M_dust_val2);
        //    }
        return dM_in_CNM + dM_in_MC + dM_in_WNM;
    }

}