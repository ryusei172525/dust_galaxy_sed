#ifndef FILE_NAME_H
#define FILE_NAME_H

#include <string>
#include <iomanip>
#include <sstream>

#include <constants.h>
#include <cfloat>

// Suffix of file names
inline static auto Suffix(const FreeParameter& free_params) noexcept {
    const auto M_GAS_INI_EXP = static_cast<int>(std::floor(
    std::log10(free_params.M_total_gal_)));
    const auto M_GAS_INI_SIG = free_params.M_total_gal_ / std::pow(10, M_GAS_INI_EXP);
    const auto TAU_SF_EXP    = static_cast<int>(std::floor(std::log10(free_params.tau_SF_)));
    const auto TAU_SF_SIG    = free_params.tau_SF_ / std::pow(10, TAU_SF_EXP);

    auto ss = std::stringstream();

    ss << std::fixed << std::setprecision(2)
       << "M=" << M_GAS_INI_SIG << "e" << M_GAS_INI_EXP
       << "Msun_R=" << free_params.R_gal_ / PC / 1e3 << "kpc_h=" << free_params.h_dust_ / PC
       << "pc_t_sf=" << TAU_SF_SIG << "e" << TAU_SF_EXP
       << "yr_wnm" << F_WNM << "_cnm" << F_CNM << "_mc" << F_MC;


    if (free_params.is_infall_) {
        const auto M_INFALL_EXP = static_cast<int>(std::floor(
        std::log10(free_params.M_total_infall_)));
        const auto M_INFALL_SIG = free_params.M_total_infall_ / std::pow(10, M_INFALL_EXP);

        const auto TAU_INFALL_EXP = static_cast<int>(std::floor(
        std::log10(free_params.tau_infall_)));
        const auto TAU_INFALL_SIG = free_params.tau_infall_ / std::pow(10, TAU_INFALL_EXP);

        ss << "_M_infall=" << M_INFALL_SIG << "e" << M_INFALL_EXP
           << "Msun_t_infall=" << TAU_INFALL_SIG << "e" << TAU_INFALL_EXP;
    }
    return ss.str();
}

/**
 * Suffix for Asano model result files
 * @param free_params
 * @return
 */
inline static auto SuffixIndependentOfGeometry(const FreeParameter& free_params) noexcept {
    const auto TAU_SF_EXP = static_cast<int>(std::floor(std::log10(free_params.tau_SF_)));
    const auto TAU_SF_SIG = free_params.tau_SF_ / std::pow(10, TAU_SF_EXP);

    auto ss = std::stringstream();

    ss << std::fixed << std::setprecision(2)
       << "t_sf=" << TAU_SF_SIG << "e" << TAU_SF_EXP
       << "yr_wnm" << F_WNM << "_cnm" << F_CNM << "_mc" << F_MC;

    if (free_params.is_infall_) {
        const auto TAU_INFALL_EXP = static_cast<int>(std::floor(
        std::log10(free_params.tau_infall_)));
        const auto TAU_INFALL_SIG = free_params.tau_infall_ / std::pow(10, TAU_INFALL_EXP);

        ss << "_t_infall=" << TAU_INFALL_SIG << "e" << TAU_INFALL_EXP;
    }
    return ss.str();
}

// Depend only on star formation timescale
inline static auto SuffixStellarSpectrum(const FreeParameter& free_params) noexcept {
    const auto M_GAS_INI_EXP = static_cast<int>(std::floor(std::log10(free_params.M_total_gal_)));
    const auto M_GAS_INI_SIG = free_params.M_total_gal_ / std::pow(10, M_GAS_INI_EXP);
    const auto TAU_SF_EXP    = static_cast<int>(std::floor(std::log10(free_params.tau_SF_)));
    const auto TAU_SF_SIG    = free_params.tau_SF_ / std::pow(10, TAU_SF_EXP);

    auto ss = std::stringstream();

    ss << std::fixed << std::setprecision(2)
       << "M=" << M_GAS_INI_SIG << "e" << M_GAS_INI_EXP << "Msun_t_sf=" << TAU_SF_SIG << "e"
       << TAU_SF_EXP << "yr";

    if (free_params.is_infall_) {
        const auto M_INFALL_EXP = static_cast<int>(std::floor(
        std::log10(free_params.M_total_infall_)));
        const auto M_INFALL_SIG = free_params.M_total_infall_ / std::pow(10, M_INFALL_EXP);

        const auto TAU_INFALL_EXP = static_cast<int>(std::floor(
        std::log10(free_params.tau_infall_)));
        const auto TAU_INFALL_SIG = free_params.tau_infall_ / std::pow(10, TAU_INFALL_EXP);

        ss << "_M_infall=" << M_INFALL_SIG << "e" << M_INFALL_EXP
           << "Msun_t_infall=" << TAU_INFALL_SIG << "e" << TAU_INFALL_EXP;
    }
    return ss.str();
}

// ################ Directories ####################################################################
// SEDモデルのディレクトリ
inline static const auto SED_MODEL_DIR = std::string(SED_MODEL_ROOT_DIR);  // CMakeで定義された値を使用
inline static const auto DATA_FILE_DIR = SED_MODEL_DIR + std::string("data/");
inline static const auto CALCULATION_RESULT_DIR = SED_MODEL_DIR + std::string("calculation_result/");
inline static const auto DUST_MODEL_RESULT_DIR = SED_MODEL_DIR + std::string("calculation_result/dust_model/");

// Dust model
inline static const auto DUST_MODEL_DIR                    = SED_MODEL_DIR + std::string("dust_model/");
inline static const auto AVERAGED_DUST_PROPERTY_RESULT_DIR = DUST_MODEL_RESULT_DIR + std::string("averaged_dust_property/");
inline static const auto ASANO_MODEL_RESULT_DIR            = DUST_MODEL_RESULT_DIR + std::string("asano_model/");
inline static const auto DEBUG_ASANO_MODEL_RESULT_DIR      = DUST_MODEL_DIR + std::string("asano_model_new2/result/");
inline static const auto ASANO_MODEL_REF_RESULT_DIR        = DUST_MODEL_DIR + std::string("asano_model_ref/result/");

// Stellar spectrum
inline static const auto STELLAR_SPECTRUM_DIR        = DATA_FILE_DIR + std::string("stellar_spectrum/");
inline static const auto STELLAR_SPECTRUM_RESULT_DIR = CALCULATION_RESULT_DIR + std::string("stellar_spectrum/");
inline static const auto STELLAR_SPECTRUM_DATA_DIR   = STELLAR_SPECTRUM_DIR + std::string("ssps_assumption/");

// mkSED
inline static const auto MKSED_DIR               = DATA_FILE_DIR + std::string("makeSED/");
inline static const auto MKSED_RESULT_DIR        = CALCULATION_RESULT_DIR + std::string("makeSED/");
inline static const auto DEBYE_RESULT_DIR        = MKSED_DIR + std::string("debye_capacity/");
inline static const auto COOLING_TIME_RESULT_DIR = MKSED_DIR + std::string("cooling_time/");

// ################ File names #####################################################################
// Dust model input files
inline std::string TotalMetalMassFileName(const std::string& metalicity) noexcept {
    return DATA_FILE_DIR + "totalmetal/newyield/totalmetal/totalmetalmass_z"
           + metalicity + ".dat";
}

inline std::string TotalRemnantMassFileName(const std::string& metalicity) noexcept {
    return DATA_FILE_DIR + std::string("totalremnant/newyield/totalremnant/totalremnantmass_z")
           + metalicity + ".dat";
}

//! Figure.13 in Nozawa et al. (2003) and Table.3 in Nozawa et al. (2007)
inline std::string SNDustMassFileName(const std::string& stellar_mass) noexcept {
    return DATA_FILE_DIR + std::string("SN_dust_mass/SN_dust_mass_m") + stellar_mass + ".dat";
}


inline std::string SNDustSizeDistributionFileName(const std::string& stellar_mass) noexcept {
    return DATA_FILE_DIR + std::string("size_distribution/shockeddust/linear/num_den")
           + std::to_string(n_H) + "/" + stellar_mass + "n" + std::to_string(n_H) + "_linear.dat";
}

inline std::string SNdustSizeDistributionBsFileName() {
    return DATA_FILE_DIR +
           std::string("size_distribution/shockeddust/linear/num_den1/snsize_BS07.dat");
}

inline std::string DestEtaFileName(const std::string& metalicity) noexcept {
    return DATA_FILE_DIR + std::string(
    "dest_eta/adjust_eta/num.n1.z") //ISM densityを変えるときはここの数字も変更する必要有 TODO:どういうこと？
           + metalicity + ".dat";
}

inline std::string AGBDustSizeDistributionFileName() {
    return DATA_FILE_DIR + std::string(
    "size_distribution/agbdust/agb_dsdis_yasudakozasa_resize_carbon.dat"); //Yasuda & Kozasa (2012)
}

inline std::string AGBDustMassFileName(const std::string& metalicity) noexcept {
    return DATA_FILE_DIR + std::string("AGB_dust_mass/AGB_dust_Z=") + metalicity + ".dat";
}

inline std::string AGBDustMassVenturaFileName(const std::string& metalicity) noexcept {
    return DATA_FILE_DIR + std::string("AGB_dust_mass/AGB_dust_ventura_Z=") + metalicity + ".dat";
}

inline static std::string
GrainVelocityFileName(const std::string& ISM_name, const std::string& dust_type) noexcept {
    //if(CNMsilvel == 2) strcpy(name,"dc2");
    //auto ISM_name = std::string();
    //if (ISM_type == 0) {
    //  ISM_name = "wnm";
    //strcpy(name,"wim_n1");//wimの速度をテストするときはwnmからwim_n1へ変更(26/11/2012)
    //} else if (ISM_type == 1) {
    //  ISM_name = "cnm";
    //strcpy(name,"dc2");//dc2の速度をテストするときはcnmからdc2へ変更(27/11/2012)
    //strcpy(name,"mc");//mcの速度をテストするときはcnmからmcへ変更(1/12/2012)
    //strcpy(name,"dc1");//dc1の速度をテストするときはcnmからdc1へ変更(1/12/2012)
    //} else if (ISM_type == 2) {
    //  ISM_name = "wim";
    //} else {
    //        ISM_name = "mc";
    //  }
    return DATA_FILE_DIR + std::string("turbulence_velocity/extrapolate_") + ISM_name + "_" +
           dust_type + "_vel.dat";
}

// Dust model output files
inline static auto
TotalDustMassFileName(const FreeParameter& free_params) noexcept {
    return ASANO_MODEL_RESULT_DIR + "total_dust_mass_" + SuffixIndependentOfGeometry(free_params) +
           ".dat";
}

inline static auto DustMassDistributionFileName(const FreeParameter& free_params) noexcept {
    return ASANO_MODEL_RESULT_DIR + "mass_dis_" + SuffixIndependentOfGeometry(free_params) + ".dat";

}

inline static auto DustNumberDistributionFileName(const FreeParameter& free_params) noexcept {
    return ASANO_MODEL_RESULT_DIR + "num_dis_" + SuffixIndependentOfGeometry(free_params) + ".dat";
}

// Dust property
inline static std::string DustPropertyDraineFileName(const std::string& dust_species) noexcept {
    return DATA_FILE_DIR + "dust_property/draine/" + dust_species + ".dat";
}

inline static std::string
DustPropertyFileName(const std::string& dust_species, const std::string& param_name) noexcept {
    return MKSED_DIR + "dust_property/" + param_name + "_" + dust_species + "_hokan.dat";
}

inline static auto AveragedDustPropertyFileName(const FreeParameter& free_params) noexcept {
    return AVERAGED_DUST_PROPERTY_RESULT_DIR + "averaged_dust_property_" + Suffix(free_params) +
           ".dat";
}

// Stellar spectrum
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "190702_SSPs.dat";
inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=1_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=0.5_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=0.1_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=0_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=-0.1_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=-0.5_SSPs.dat";
//inline static const auto SSPs_FILE_NAME = STELLAR_SPECTRUM_DATA_DIR + "Chabrier_x=-1_SSPs.dat";

inline static auto AgeString(std::size_t ai) noexcept { return "_age" + std::to_string(ai); }

[[nodiscard]] inline static auto
StellarContinuumFileName(const FreeParameter& free_params) noexcept {
    return STELLAR_SPECTRUM_RESULT_DIR + "stellar_continuum_" + SuffixStellarSpectrum(free_params)
           + AgeString(free_params.ai_) + ".dat";
}

inline static auto FYoungFileName(const FreeParameter& free_params) noexcept {
    return STELLAR_SPECTRUM_RESULT_DIR + "f_young_" + SuffixStellarSpectrum(free_params)
           + AgeString(free_params.ai_) + ".dat";
}

inline static auto SuffixMkSED(const FreeParameter& free_params) noexcept {
    return Suffix(free_params) + AgeString(free_params.ai_) + ".dat";
}

inline static auto TransmissionRateFileName(const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "transmission_rate_" + SuffixMkSED(free_params);
}

inline static auto MGAParameterFileName(const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "mga_parameter" + SuffixMkSED(free_params);
}

inline static std::string
HitTimeFileName(const std::string& dust_species, const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "hit_time_" + dust_species + "_" + SuffixMkSED(free_params);
}

inline static std::string EnthalpyPerAtomFileName(const std::string& dust_species) noexcept {
    if (dust_species == "Sil") {
        return DEBYE_RESULT_DIR + "enthalpy_per_atom_Sil.dat";
    } else if (dust_species == "Gra") {
        return DEBYE_RESULT_DIR + "enthalpy_per_atom_Gra.dat";
    } else {
        return DEBYE_RESULT_DIR + "enthalpy_per_atom_PAH.dat";
    }
}

inline static std::string CoolingTimeFileName(const std::string& dust_species,
                                              const std::string& dust_radius_string) noexcept {
    return COOLING_TIME_RESULT_DIR + dust_species + "/" + dust_radius_string + ".dat";
}

inline static std::string
HistogramFileName(const std::string& dust_species, const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "temperature_histogram_" + dust_species + "_" +
           SuffixMkSED(free_params);
}

inline static std::string
DustRadiationFileName(const std::string& dust_species, const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "dust_radiation_" + dust_species + "_" + SuffixMkSED(free_params);
}

inline static std::string
DustRadiationPerRadiusFileName(const std::string& dust_species,
                               const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "dust_radiation_per_radius_" + dust_species + "_" +
           SuffixMkSED(free_params);
}

[[nodiscard]] inline static auto SEDFileName(const FreeParameter& free_params) noexcept {
    return MKSED_RESULT_DIR + "SED_" + SuffixMkSED(free_params);
}

inline static const auto INOUE_FYOUNG_FILE_NAME =
                         MKSED_RESULT_DIR + "young_star_fraction_inoue2005.dat";
inline static const auto FYOUNG_FILE_NAME       = MKSED_RESULT_DIR + "young_star_fraction.dat";

#endif // FILE_NAME_H
