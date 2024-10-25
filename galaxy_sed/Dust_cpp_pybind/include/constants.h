#ifndef CONSTANT_H
#define CONSTANT_H
#define _USE_MATH_DEFINES

#include <array>
#include <valarray>
#include <cmath>
#include <string>
#include <iostream>
#include <free_parameter.h>

// Free parameters
inline static constexpr auto TIME_BIN = 1e8; // bin scale of galaxy time in output files[yr] 1e7=10Myr  default:1e8
//inline static constexpr auto N_MAX_AGE = static_cast<std::size_t>(MAX_AGE / TIME_BIN);

//inline static constexpr auto N_AGE = CALCULATE_AGE + 1;

inline static constexpr auto TIME_BIN_DUST_MODEL = 1e6; // When SF timescale smaller than 1e8, change this value to 1e6; default:1e7 1e6

//inline static constexpr auto N_AGE_ =
//                             static_cast<std::size_t>(MAX_AGE / TIME_BIN_DUST_MODEL) + 1; // [yr]

inline static constexpr auto INITIAL_Z     = 0.0; // initial metalicity
inline static constexpr auto F_SUB_STELLAR = 0.0; // Mass fraction of substellar objects

// Regarding to Radiative transfer
inline static constexpr auto H_DUST_MW = 150 * PC; // from Inoue 2005, h_d = 150 pc

//inline static constexpr auto H_DUST     = 150 * PC; // tamura_san
inline static constexpr auto R_GAL_MW = 10e3 * PC; // MW: Fixed galaxy radius

//inline static constexpr auto R_GAL      = 10e3 * PC; // tamura_san
inline static constexpr auto Volume_MW = R_GAL_MW * R_GAL_MW * 3.14 * H_DUST_MW;

//inline static constexpr auto M_TOTAL_GAL = 2.5e9; // LMC
//inline static constexpr auto R_GAL         = 3.65e3 * PC; // LMCr
//inline static constexpr auto H_STAR      = 180 * PC; // LMC

// Regarding to Asano model
inline static constexpr auto F_WNM = 0.5; // mass fraction of the WNM where shattering and/or coagulation occur default:0.5  0.3
inline static constexpr auto F_CNM = 0.3; // mass fraction of the CNM where accretion and/or coagulation and/or shattering occur default:0.3 0.0
inline static constexpr auto F_MC  = 0.2; // mass fraction of the MC where accretion and/or shattering and/or coagulation occur default:0.2  0.7
//inline static constexpr auto F_WNM = 0.0; // mass fraction of the WNM where shattering and/or coagulation occur
//inline static constexpr auto F_CNM = 0.0; // mass fraction of the CNM where accretion and/or coagulation and/or shattering occur
//inline static constexpr auto F_MC  = 1.0; // mass fraction of the MC where accretion and/or shattering and/or coagulation occur

// ############################################################################################################
// Physical constants
inline static constexpr auto Msun     = 2.0e33; // solar mass
inline static constexpr auto Zsun     = 0.02; // solar metallicity
inline static constexpr auto m_H      = 1.673e-24; //hydrogen atom mass g
inline static constexpr auto MU       = 1.33; // mean atomic weight (corresponding to Y = 0.25) 2012/10/22
inline static constexpr auto kb       = 1.38e-16; //boltzmann constant
inline static constexpr auto PI       = 3.141592654;
inline static constexpr auto cl       = 2.99792458e10; // [cm / s]; Speed of light
inline static constexpr auto h_P      = 6.626e-27; //[erg s]; Planck's constant
inline static constexpr auto k_B      = 1.381e-16; //[erg / K]; Boltzmann's constant
inline static constexpr auto GRAV     = 6.67e-8; //[cm3 / s2 / g]; Gravitational constant
inline static constexpr auto m_PROTON = 1.673e-24; //[g]; proton mass

inline static constexpr auto L_sun = 3.839e33; //[erg / s]; Solar luminosity
//inline static constexpr auto M_GALAXY = M_TOTAL_GAL; // [Msun]; Total galaxy mass of baryon (stellar and gas)
inline static constexpr auto N_A   = 6.022e23; // Avogadro const.

// Numerical constants
inline static constexpr auto N_DUST_RADIUS      = std::size_t(45);
inline static constexpr auto N_PAH_RADIUS       = std::size_t(16);
inline static constexpr auto N_LAMBDA           = std::size_t(2700);
inline static constexpr auto N_SILICATE_SPECIES = std::size_t(9); // Under 10 species
inline static constexpr auto N_DUST_SPECIES     = N_SILICATE_SPECIES + 1; // Under 10 species

// Bulk density (g/cm-3) 0:C, 1:Si, 2:Fe, 3:FeS, 4:Al2O3, 5:MgSiO3, 6:Mg2SiO4, 7:SiO2, 8:MgO, 9:Fe3O4, 10:Silicate
inline static const auto RHO_GRAIN = std::valarray<double>{2.2631, 2.3244, 7.8913, 4.8358, 3.9824,
                                                           3.1813, 3.2018, 2.6428, 3.5633, 5.2073,
                                                           3.3};
// The atomic mass of [0: carbon, 1: silicate]
inline static const auto M_ATOM    = std::valarray<double>{12.0, 28.1};

//########################## mkSED constants ##############################################################
// Number of PAHion and PAHneu dust radius steps
static constexpr auto N_DUST_RADIUS_PAH = std::size_t(16);


// Number of angle. This should even number for upstream and downstream numbers are equal.
static constexpr auto N_ANGLE              = std::size_t(32);
static constexpr auto Ng_CRITERION         = 1e-10; // Criterion of Ng-iteration convergence
static constexpr auto N_Ng_MAX             = std::size_t(50); // Number of max Ng-iteration
//0:10Myr 1:100Myr 2:1Gyr // young star threshold
static constexpr auto YOUNG_STAR_THRESHOLD = std::size_t(0); // 0: 10 Myr, 1: 100 Myr 3: 1 Gyr
static constexpr auto INDEX_LYMAN_BREAK    = std::size_t(588); // Index number of 912A

static constexpr auto A_C   = 3.2e-3; //[cm] Coefficient of equation of calculating carbon equilibrium temperature
static constexpr auto A_Sil = 1.4e-3;  //[cm] Drapatz & Michel 1977 Coefficient of equation of calculating Silicate equilibrium temperature
//static constexpr auto C_Sil   = 3.59e23; // 3.59e23 [cm-3] Coefficient of atom number related in Nasimoto's master paper
static constexpr auto C_Sil = 4.38e23; // 3.59e23 [cm-3] Coefficient of atom number related in Nasimoto's master paper
//static constexpr auto C_Sil = 3.50 * 4 / 3 * PI / 172.25 * N_A; // 5.12e22
static constexpr auto C_Gra = 4.68e23; //4.68e23 [cm-3] Coefficient of atom number related in Nasimoto's master paper
//static constexpr auto C_Gra = 2.24 * 4 / 3 * PI / 12 * N_A; // 4.71e23
//static constexpr auto N_c    = 1.14e23;        //[cm^-3] number density of carbon atom
//static constexpr auto N_atom = 8.5e22;        //[cm^-3] number density of silicate atom by takeuchi et al. 2003
//static constexpr auto C_Sil   = N_atom * 4 / 3 * PI; // 3.56e23
//static constexpr auto C_Gra   = N_c * 4 / 3 * PI; // 4.78e23



static constexpr auto dMU    = 2.0 / (N_ANGLE - 1); // Distance of one step of mu-axis
static constexpr auto dTHETA = M_PI / (N_ANGLE - 1); // Distance of one step of mu-axis

// Constants of used in monte-calro simulation.
static constexpr auto N_MC_MAX = std::size_t(2e9); // Number of max MC iteration
//static constexpr auto STEP_LAMBDA_CM  = (LAST_LAMBDA_CM - FIRST_LAMBDA_CM) / N_DIV_LAMBDA;
static constexpr auto N_HIT    = 10000; //
static constexpr auto N_REBIN  = 100;

//!Dust properties//////////////////////////////////////////////////////////////////////////
static constexpr auto GAMMA = 0.5;  //!Extended MGA parameter; 0.5 for slab(Varosi & Dwek 1999)


//static constexpr auto CARBON_SUBLIMATION_TEMPERATURE   = 5000; // 2500 K
//static constexpr auto SILICATE_SUBLIMATION_TEMPERATURE = 5000; // 1500 K
static constexpr auto CARBON_SUBLIMATION_TEMPERATURE   = 2500; // 2500 K
static constexpr auto SILICATE_SUBLIMATION_TEMPERATURE = 1500; // 1500 K
static constexpr auto N_TEMPERATURE                    = 5000;
static constexpr auto T_MIN                            = 0.1;

// Constants of used in mega-grain approximation
static const auto PKB = std::pow(10, 3.5); //[K / cm3] ISM equilibrium thermal pressure (p/kB)

static constexpr auto MEAN_ATOMIC_WEIGHT = 1.4; //!Mean atomic weigh


// hydrogen number densities of the CNM and the WNM [cm-3]
static const auto     INDEX_CNM = 0.7;
static const auto     INDEX_WNM = 1.0;
static const auto     n0_CNM    = 1e5; // default: 1e3 1e4 1e5(kano)
static const auto     n0_WNM    = 1.0; // default: 1.0
//![K / cm3]; WIM temperature (Warm ionized medium)
static const double   T_CNM_MGA = std::pow(10, 4.5);
// ~ 37.2 [cm-3] hydrogen density in cold neutral medium (Inoue 2005)
static const auto     n_CNM     = n0_CNM * std::pow(PKB / T_CNM_MGA, 1.0 / INDEX_CNM);
//[K / cm3] WIM temperature (Warm ionized medium)f
static constexpr auto T_WNM_MGA = 1e+4;
//~0.316 [cm-3] hydrogen density in warm neutral medium
static const auto     n_WNM     = n0_WNM * std::pow(PKB / T_WNM_MGA, 1.0 / INDEX_WNM);

//![cm - 3]; Homogeneous gas number density
static constexpr auto n_HOMOGENEOUS = 1.0;

// [g / cm3] Regard WNM as interclump medium
static const auto RHO_ICM   = MEAN_ATOMIC_WEIGHT * m_PROTON * n_WNM;
// [g / cm3] Regard CNM as clump
static const auto RHO_CLUMP = MEAN_ATOMIC_WEIGHT * m_PROTON * n_CNM;

// volume filling fraction of clumps (~1.85e-2)
// a + b = 1, a * n_CNM + b * n_WNM = n_HOMOGENEOUS の解
static const auto FRACTION_CLUMP = (n_HOMOGENEOUS - n_WNM) / (n_CNM - n_WNM);

//[cm] Jeans length for a sphere of the CNM
// static const auto R_CLUMP = std::sqrt(15.0 * k_B * PKB / (4.0 * PI * GRAV)) /
//                             RHO_CLUMP;  // calculate R_CLUMP/10 (I want to scale 1/10)

// FIX R_CLUMP
static const auto R_CLUMP = 3.20201e+19;

//static const auto R_CLUMP = std::sqrt(15.0 * k_B * PKB / (4.0 * PI * GRAV))
//                            / RHO_CLUMP * std::pow((1.0 - FRACTION_CLUMP), GAMMA);

// ######################## Asano model ###############################################################################
// Constants of iteration numbers
inline static constexpr auto ACC_TIME_SCALE = 1e7; // accretion growth time-scale

inline static constexpr auto N_TIME_BIN_SHORT = 1; // the number of timebin for shattering and coagulation.

//inline static constexpr auto N_TIME_BIN       = 101; // the number of timebin.
//inline static constexpr auto HH = 1e4; // short timestep, 10^4yr from Hirashita-san's code
inline static constexpr auto HH = TIME_BIN_DUST_MODEL /
                                  N_TIME_BIN_SHORT; // short timestep, 10^4yr from Hirashita-san's code
inline static constexpr auto YR = 365 * 24 * 60 * 60.0; // conversion factor yr to sec
inline static constexpr auto DT = HH * YR; // delta t unit on sec

inline static constexpr auto N_METALLICITY     = std::size_t(
1001); // The number of metallicity index in metal mass files and remnant mass files
inline static constexpr auto N_M_STAR          = std::size_t(
3901); // The number of mass in metal mass files and remnant mass files
//inline static constexpr auto N_M_STAR          = std::size_t(11901); // The number of mass in metal mass files and remnant mass files
inline static constexpr auto N_M_SN_STAR       = std::size_t(4); // The number of metalicity fo SN
inline static constexpr auto N_M_AGB_STAR      = std::size_t(27); // The number of metalicity fo SN
inline static constexpr auto N_METALICITY_AGB  = std::size_t(6); // The number of metalicity of AGB
// Minimum size of grains 0:10^-7.9cm, 4:10^-7.5cm, 9:10^-7.0cm, 16:10^-6.3cm
inline static constexpr auto N_RADIUS_MIN      = std::size_t(4);
// The index number of dust radius. 49:10^-3.1cm. The number of dust radius is N_DUST_RADIUS - N_RADIUS_MIN;
inline static constexpr auto N_MAX_DUST_RADIUS = std::size_t(49);

//inline static constexpr auto SCHMIDT_INDEX = 1; // Schmidt index
inline static constexpr auto M_CH       = 0.35; // characteristic mass used in IMF function.
inline static constexpr auto A_exp      =0.158;
inline static constexpr auto mc         =0.079;
inline static constexpr auto sigma      =0.69;
inline static constexpr auto A_pow      =4.43e-2;
inline static constexpr auto n_H        = 1; // Hydrogen number density of surrounding ISM of a SN
//inline static const auto     DM            = std::pow(10.0, 0.15) - std::pow(10.0, -0.15); // The size of mass
inline static const auto     DM         = 3 * std::pow(10.0, 0.1); // たぶんこっちが正しい
inline static constexpr auto RHO_CARSIL = std::array<double, 2>{2.2631,
                                                                3.3}; // carbon and silicate
// The switch of Data type used in code
inline static constexpr auto AGB_DATA_TYPE    = 0; // 0: Zhukovska et al.(2008) 1: Ventura et al.(2012a,b)
inline static constexpr auto SN_DATA_TYPE     = 0; // 0:SN size distribution data written in Nozawa et al.(2007), 1:the data written in Bianchi & Schneider = 2007)
inline static constexpr auto DUST_SOURCE_TYPE = 0; // 0:AGB+SNeII, 1:AGB only, 2:SNeII only
inline static constexpr auto IMF_TYPE         = 2; // 0:Salpeter IMF, 1:Larson IMF, 2:Chabrier IMF, 3:Chabrier_1+x=2, 4:Chabrier_1+x=1.5, 5:Chabrier_1+x=1.1 7:Chabrier_1+x=0.9

// Used in stellar spectrum (PEGASE)
inline static constexpr auto isGALACTIC_WIND    = false;
inline static constexpr auto isNEBULAR_EMISSION = true;
inline static constexpr auto SFR_TYPE           = 3; // 3: Schmitd low
inline static constexpr auto T_WIND             = 200003;
inline static constexpr auto F_SNIa             = 0.0;

class ISMParams {
  public:
    std::string name_   = "NAN";
    bool        isACC_  = false; // Accretion
    bool        isCoag_ = false; // Coagulation
    bool        isShat_ = false; // Shattering
    double      n_      = 0.0; // Number density [cm^-3]
    double      T_      = 0.0; // temperature [K]
    double      f_      = 0.0; // Mass fraction of this phase in ISM

    ISMParams() = delete;

    ISMParams(const std::string& ISM_name) : name_(ISM_name) {
        if (name_ == "CNM") SetCNMParams();
        else if (name_ == "MC") SetMCParams();
        else if (name_ == "WNM") SetWNMParams();
        else if (name_ == "WIM") SetWIMParams();
        else {
            std::cout << "Error: name is only accepted as CNM, WNM, WIM or MC!" << std::endl;
        }
    }

    void SetCNMParams() {
        isACC_  = true;     //default:true
        isCoag_ = true;     //default:true
        isShat_ = true;     //default:true
        n_      = 30.0;
        T_      = 100.0;
        f_      = F_CNM;
    }

    void SetMCParams() {
        isACC_  = true;    //default:true
        isCoag_ = true;    //default:true
        isShat_ = true;    //default:true
        n_      = 300.0;
        T_      = 25.0;
        f_      = F_MC;
    }

    void SetWNMParams() {
        isACC_  = false;    //default:false
        isCoag_ = true;    //default:true
        isShat_ = true;    //default:true
        n_      = 0.3;
        T_      = 1e4;
        f_      = F_WNM;
    }

    void SetWIMParams() {
        isACC_  = false;
        isCoag_ = false;
        isShat_ = false;
        n_      = 0.1;
        //T        = 100.0;
        f_      = 0.0;
    }
};

//basic parameters for calculation of dust process
inline static constexpr auto isSN_ADJUST = false; // 0:original SN dust mass, 1:0.1*dust mass supply by SNe
inline static constexpr auto isINJECTION = true; // Injected by stars
inline static constexpr auto isSN_DEST   = true; // SN destruction effect　default:true
inline static constexpr auto isCNM       = true; // switch whether dust processes occurs in the CNM or not
inline static constexpr auto isMC        = true; // switch whether dust processes occurs in the MC or not
inline static constexpr auto isWIM       = false; // switch whether dust processes occurs in the WIM or not
inline static constexpr auto isWNM       = true; // switch whether dust processes occurs in the WNM or not

//accretion on/ off 切り替え true:on false:off
inline static constexpr auto isACC_CNM = true; // accretion in the CNM   default:true
inline static constexpr auto isACC_MC  = true; // accretion in the MC    default:true
inline static constexpr auto isACC_WNM = false; // accretion in the WNM  default:false
inline static constexpr auto isACC_WIM = false; // accretion in the CNM  default:false

//shattering on/ off 切り替え true:on false:off
inline static constexpr auto isSHAT_CNM = false; // switch whether shattering occurs or not in the CNM  default:true
inline static constexpr auto isSHAT_MC  = false; // switch whether shattering occurs or not in the MC   default:true
inline static constexpr auto isSHAT_WIM = false; // switch whether shattering occurs or not in the WIM default:false
inline static constexpr auto isSHAT_WNM = false; // switch whether shattering occurs or not in the WNM  default:true

//coagulation on/ off 切り替え true:on false:off
inline static constexpr auto isCOAG_CNM = false; // switch whether shattering occurs or not in the CNM  default:true
inline static constexpr auto isCOAG_MC  = false; // switch whether shattering occurs or not in the MC   default:true
inline static constexpr auto isCOAG_WIM = false; // switch whether shattering occurs or not in the WIM default:false
inline static constexpr auto isCOAG_WNM = false; // switch whether shattering occurs or not in the WNM  default:true

// ISM parameters
inline static constexpr auto F_WIM = 0.0; // mass fraction of the WIM where shattering and/or coagulation occur
inline static constexpr auto N_CNM = 30.0; // number density in the CNM

inline static constexpr auto N_WIM = 0.1; // number density in the WIM
inline static constexpr auto N_WNM = 0.3; // number density in the WNM
inline static constexpr auto T_CNM = 100.0; // temperature in the CNM
inline static constexpr auto N_MC  = 300.0; // number density in the MC
inline static constexpr auto T_MC  = 25.0; // temperature in the MC

// The coefficients of shattering and coagulation process
inline static constexpr auto SHAT_RADIUS_INDEX_MIN = 7.5; // The index of minimum size of shattered grains. 7.5:10^-7.5cm
inline static constexpr auto F_GM                  = std::array<double, 2>{1.0,
                                                                           0.166}; // Key element (C or Si) Mass fraction in grain by Hirashita & Kuo (2011)
inline static constexpr auto R_COLL                = 1.0; // For the collision between the same species
inline static constexpr auto Z_COLL                = 3.4; // The radial velocity of the cratering flow in the shattered material in JTH96
inline static constexpr auto F_M                   = 0.4; // ejection fraction of shocked mass by cratering, 0.4 is fiduciary value (Hirashita & Kobayashi 2013;
inline static constexpr auto f_catastrophic        = 1.0; // ejection fraction of shattered fragments by catastrophic, 1.0 is fiduciary value.
inline static constexpr auto X_SHAT                = (Z_COLL * 4.0 + 1.0) / (Z_COLL + 1.0); // ~3.3
inline static constexpr auto ALPHA                 = std::array<double, 2>{1.0, 1.0};
inline static constexpr auto V_SHAT                = std::array<double, 2>{1.2e5, 2.7e5};

inline static constexpr auto F_stick = 10.0;

// Used in shattering constants
inline static constexpr auto P1     = std::array<double, 2>{4.0e10,
                                                            3.0e11}; // Critical pressure (carbon and silicate) by Jones et al. 1996
inline static constexpr auto c0     = std::array<double, 2>{1.8e5,
                                                            5.0e5}; // Sound speed (carbon and silicate)
inline static constexpr auto s_shat = std::array<double, 2>{1.9,
                                                            1.2}; // Dimensionless material constants
inline static constexpr auto Pv     = std::array<double, 2>{5.8e12, 5.4e12};

#endif // CONSTANT_H
