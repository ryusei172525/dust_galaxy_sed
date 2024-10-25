#ifndef SED_MODEL_FREE_PARAMETER_H
#define SED_MODEL_FREE_PARAMETER_H

#include <cmath>

inline static constexpr auto PC = 3.09e18; //[cm]; 1 parsec

static constexpr auto XI    = 0.5;
static constexpr auto N_Z   = std::size_t(101);  // Number of z-axis steps
static constexpr auto N_TAU = N_Z; // Number of the optical depth steps 2dtau = tauhom / N_TAU less than 0.1

class FreeParameter {
  public:
    bool        is_infall_ = true;
    double      age_max_   = 13e9;
    std::size_t ai_        = std::size_t(129);
    std::size_t n_age_     = static_cast<std::size_t >(age_max_ / 1e8);

    double M_total_gal_    = 1e11;
    double tau_SF_         = 5e9; //default:3e9
    double schmidt_index_  = 1.0;
    double M_total_infall_ = M_total_gal_;
    double tau_infall_     = 5e9; //default:4e9

    // Geometrical parameter
    double h_dust_ = 150 * PC;
    double R_gal_  = 20e3 * PC;

    //double h_star_ = h_dust_ / xi_; // MW
    double S_gal_ = 2 * M_PI * R_gal_ * R_gal_ + 2 * M_PI * R_gal_ * h_dust_ * 2;
    double V_gal_ = R_gal_ * R_gal_ * M_PI * h_dust_;
    double dz_    = h_dust_ / (N_Z - 1); // Distance of one step of z-axis

    /**
     * Free parameter list class
     */
    FreeParameter() = default;

    FreeParameter(bool is_infall, double M_total_gal, double tau_SF, double schmidt_index,
                  double M_total_infall, double tau_infall, double age_max, std::size_t ai,
                  std::size_t n_age)
    : is_infall_(is_infall),
      age_max_(age_max),
      ai_(ai),
      n_age_(n_age),
      M_total_gal_(M_total_gal),
      tau_SF_(tau_SF),
      schmidt_index_(schmidt_index),
      M_total_infall_(M_total_infall),
      tau_infall_(tau_infall) {}

    ~FreeParameter() = default;

    inline void SetIsInfall(bool is_infall) noexcept { is_infall_ = is_infall; }

    inline void SetAgeMax(double age_max) noexcept {
        age_max_ = age_max;
        n_age_   = static_cast<std::size_t >(age_max_ / 1e8);
    }

    inline void SetIndexOfAge(std::size_t ai) noexcept { ai_ = ai; }

    inline void SetNumberOfAge(std::size_t n_age) noexcept { n_age_ = n_age; }

    inline void SetTotalGalaxyMass(double M_total_gal) noexcept {
        M_total_gal_    = M_total_gal;
        M_total_infall_ = M_total_gal;
    }

    inline void SetStarFormationTimescale(double tau_SF) noexcept { tau_SF_ = tau_SF; }

    inline void SetSchmidtIndex(double schmidt_index) noexcept { schmidt_index_ = schmidt_index; }

    inline void
    SetTotalInfallMass(double M_total_infall) noexcept { M_total_infall_ = M_total_infall; }

    inline void SetInfallTimescale(double tau_infall) noexcept { tau_infall_ = tau_infall; }

    inline void SetDustScaleHeight(double h_dust) noexcept {
        h_dust_ = h_dust;
        //h_star_ = h_dust_ / xi_;
        V_gal_  = R_gal_ * R_gal_ * 3.14 * h_dust;
        dz_     = h_dust / (N_Z - 1);
        S_gal_  = 2 * M_PI * R_gal_ * R_gal_ + 2 * M_PI * R_gal_ * h_dust * 2;
    }

    inline void SetGalaxyRadius(double R_gal) noexcept {
        R_gal_ = R_gal;
        V_gal_ = R_gal_ * R_gal_ * 3.14 * h_dust_;
        S_gal_ = 2 * M_PI * R_gal_ * R_gal_ + 2 * M_PI * R_gal_ * h_dust_ * 2;
    }
};


//inline static constexpr auto MAX_AGE       = 13e9; // [yr]
//inline static constexpr auto CALCULATE_AGE = std::size_t(129); // [0, 1, 2, 3, 4]


// According to star formation history
//inline static constexpr auto M_TOTAL_GAL   = 1.0e11; // [Msun]
//inline static constexpr auto SF_TIME_SCALE = 3e9; // [yr]

//inline static constexpr auto isInfallModel     = true; // 0: model, 1:infall model
//inline static constexpr auto M_INFALL_GAS      = M_TOTAL_GAL; // Total infall mass in galactic age M82
//inline static constexpr auto INFALL_TIME_SCALE = 6e9;
inline static constexpr auto Z_INFALL = 0.0; // Metalicity of infall gas

inline static constexpr auto isBURST     = false ; // 0=true:closed-box model, 1=false:infall model   default:false
inline static constexpr auto M_BURST     = 2e8; // M82
inline static constexpr auto T_BURST     = 12.95e9; // M82
inline static constexpr auto T_E_FOLDING = 5e6; // M82


#endif //SED_MODEL_FREE_PARAMETER_H
