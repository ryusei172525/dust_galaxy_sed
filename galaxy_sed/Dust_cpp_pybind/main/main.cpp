/**
 * @file SED_mode.cpp
 * @brief Calculate SED model
 */

#include <SED_model.h>

#include <SED_setting_util.h>
#include <SED_file_util.h>
#include <dust_util.h>
#include <thread>

#include <asano_model.h>
#include <averaged_dust_propery.h>

#include <pybind11/pybind11.h>

using valt = std::valarray<std::size_t>;
using val = std::valarray<double>;
using val2 = std::valarray<val>;
using val3 = std::valarray<val2>;

using vec = std::vector<double>;

using SEDSetting = my_util::SED_util::SEDSettingUtil;
using SEDFile = my_util::SED_util::SEDFileUtil;
using DustProperty = dust_property::AveragedDustProperty;
using Time = my_util::MyTime;

// 既存のmain関数の内容を新しい関数に移動
void SED_calculator(int galaxy_age) {
    const auto start_time = Time::GetTime();
    const auto dust_species_vec = std::vector<std::string>{"Sil", "Gra", "PAHneu", "PAHion"};

    const auto i_age_valt = valt{galaxy_age};      // galaxy age,0=100Myr,4=500Myr,9=1Gyr,49=5Gyr,99=10Gyr,129=13Gyr
    const auto M_gal_val = val{1.7e11};   //! Total galaxy mass list [Msun]  default:1.0e11   2e10 2.9e12(kawamoto distant assumption) 3.1e11(kano)　1.7e11(kano)
    const auto tau_SF_vec = vec{1e8};     //! Star formation timescale list  [yr] default:3e9  5e8(kawamoto distant assumption) 1e8(kano)
    const auto tau_infall_vec = vec{5e8}; //! Infall timescale [yr]  default:15e9   1e9(kawamoto distant assumption) 5e8(kano)
    const auto tau_pair_vec = SED_model::MakePairVector(tau_SF_vec, tau_infall_vec);

    const auto h_dust_vec = vec{150 * PC}; //! Dust scale height [cm] default:150PC
    const auto R_gal_vec = vec{10e3 * PC}; //! Galaxy radius [cm]    default:10e3PC(10kPC)
    const auto geometry_pair_vec = SED_model::MakePairVector(h_dust_vec, R_gal_vec);

    auto free_params = FreeParameter();
    free_params.SetIsInfall(true); //  false:closed-box model, true:infall model
    free_params.SetAgeMax((static_cast<double>(i_age_valt.max() + 1)) * TIME_BIN);

    const auto a_cm_val = SEDSetting::DustRadiusCm();
    const auto lambda_cm_val = SEDSetting::LambdaCm();
    const auto E_photon_val = cl * h_P / lambda_cm_val;

    //! dust class vector
    auto dust_vec = my_util::dust_util::MakeDustClassVector(dust_species_vec, a_cm_val);

    const auto asano_model = asano_model::AsanoModel();
    const auto stellar_spectrum = stellar_spectrum::StellarSpectrum(lambda_cm_val);

    // first: tau_SF, second: tau_infall
    for (const auto &tau_pair : tau_pair_vec)
    {
        free_params.SetStarFormationTimescale(tau_pair.first);
        free_params.SetInfallTimescale(tau_pair.second);

        const auto fn_m_total = TotalDustMassFileName(free_params);
        const auto fn_m = DustMassDistributionFileName(free_params);
        const auto fn_n = DustNumberDistributionFileName(free_params);
        asano_model.Calculate(free_params, fn_m_total, fn_m, fn_n);

        // used in radiative transfer
        const auto D_val = SEDFile::ReadDustToGasMassRatio(fn_m_total, free_params.n_age_);

        const auto M_total_gal_val = SEDFile::ReadTotalGalaxyMassRatio(fn_m_total,
                                                                       free_params.n_age_);

        // n_val3[0]: silicate, n_val3[1]: graphite, n_val3[2]: neutral PAH, and n_val3[3]: ionized PAH
        const auto n_val3 = my_util::dust_util::ReadAndDivideDustNumberDistribution(
            fn_n, free_params.n_age_, a_cm_val);

        // dust extinction per unit mass, scattering albedo, and asymmetry parameter
        const auto fn_averaged_params = AveragedDustPropertyFileName(free_params);
        const auto dust_params_val3 = DustProperty::Calculate(a_cm_val, free_params, n_val3,
                                                              fn_averaged_params);

        const auto Cabs_sum_val2 = SEDFile::ReadSumAbsorptionCoefficient(fn_averaged_params,
                                                                         free_params.n_age_);

        for (const auto &i_age : i_age_valt)
        {
            free_params.SetIndexOfAge(i_age);

            stellar_spectrum.Calculate(free_params);
            const auto f_young_val = SEDFile::ReadYoungStellarFraction(FYoungFileName(free_params));
            const auto L_star_val = SEDFile::ReadStellarContinuumFile(
                StellarContinuumFileName(free_params));

            for (const auto &M_gal : M_gal_val)
            {
                free_params.SetTotalGalaxyMass(M_gal);

                for (auto i = std::size_t(0); i < dust_vec.size(); ++i)
                    dust_vec[i].SetNumberDistribution(M_gal, n_val3[i]);

                SED_model::AgeDependent(dust_species_vec, M_gal,
                                        lambda_cm_val, a_cm_val, E_photon_val, D_val[i_age],
                                        dust_params_val3[0][i_age],
                                        dust_params_val3[1][i_age], dust_params_val3[2][i_age],
                                        dust_vec, free_params, geometry_pair_vec,
                                        L_star_val * M_gal, f_young_val, Cabs_sum_val2[i_age]);
            }
        }
    }

    my_util::MyTime::PrintElapsedTime(start_time, "sec", "Total: ");
    std::cout << std::endl;
}

// int main(int argc, char *argv[]) {
//     SED_calculator(galaxy_age);
// }

// pybind11モジュールの定義
namespace py = pybind11;

PYBIND11_MODULE(sed_module, m) {
    m.def("SED_calculator", &SED_calculator, "Run the sed calculation");
}