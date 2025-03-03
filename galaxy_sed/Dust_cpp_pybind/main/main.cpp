/**
 * @file SED_mode.cpp
 * @brief SEDモデルを計算する
 */

#include <SED_model.h>

#include <SED_setting_util.h>
#include <SED_file_util.h>
#include <dust_util.h>
#include <thread>
#include <asano_model.h>
#include <averaged_dust_propery.h>
#include <string>
#include <constants.h>
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

// ハイパーパラメータを格納する構造体
struct alignas(8) Hyperparameters {
    int galaxy_age;                // 銀河の年齢 0=100Myr,4=500Myr,9=1Gyr,49=5Gyr,99=10Gyr,129=13Gyr
    double galaxy_mass;            // 銀河の質量
    double starformation_timescale; // 星形成のタイムスケール
    double gas_infall_timescale;   // ガス流入のタイムスケール
    double dust_scale_height;      // ダストのスケールハイト
    double galaxy_radius;          // 銀河の半径
    double n0_cnm;                 // 中性水素の密度
    bool is_infall;                // 流入モデルかどうか
    int imf_type;                  // 初期質量関数のタイプ
    // 他のハイパーパラメータを追加
};

// SED計算を行う関数
void SED_calculator(const std::shared_ptr<Hyperparameters>& mcmc_params) {
    const auto start_time = Time::GetTime(); // 計算開始時間を記録

    // STEP0. パラメータの初期化と物理定数の計算
    const auto i_age_valt = valt{mcmc_params->galaxy_age};      // 銀河の年齢を取得
    const auto M_gal_val = val{mcmc_params->galaxy_mass * 1e9};   //! 銀河の質量を太陽質量に変換 [Msun] 天の川銀河:1.0e11 [Msun]
    const auto tau_SF_vec = vec{mcmc_params->starformation_timescale * 1e9};     //! 星形成のタイムスケールリスト [yr] 天の川銀河:3e9[yr]
    const auto tau_infall_vec = vec{mcmc_params->gas_infall_timescale * 1e9}; //! ガス流入のタイムスケール [yr] 天の川銀河: 15e9[yr]
    const auto h_dust_vec = vec{mcmc_params->dust_scale_height * PC}; //! ダストのスケールハイト [cm]　天の川銀河: 150pc
    const auto R_gal_vec = vec{mcmc_params->galaxy_radius * PC}; //! 銀河の半径 [cm] 天の川銀河: 10e3pc
    const double n0_cnm_ = mcmc_params->n0_cnm; // 中性水素の密度を取得
    const auto IMF_TYPE = mcmc_params->imf_type; // 初期質量関数のタイプを取得

    const auto tau_pair_vec = SED_model::MakePairVector(tau_SF_vec, tau_infall_vec); // タイムスケールのペアを作成
    const auto geometry_pair_vec = SED_model::MakePairVector(h_dust_vec, R_gal_vec); // ジオメトリのペアを作成

    // フリーパラメータの初期化
    auto free_params = FreeParameter(); // フリーパラメータのインスタンスを作成
    free_params.SetIsInfall(mcmc_params->is_infall); // 流入モデルかどうかを設定 false:closed-box model, true:infall model
    free_params.SetAgeMax((static_cast<double>(i_age_valt.max() + 1)) * TIME_BIN); // 最大年齢を設定

    // ダスト半径リストと波長リストを作成
    const auto a_cm_val = SEDSetting::DustRadiusCm(); // ダスト半径をセンチメートルで取得
    const auto lambda_cm_val = SEDSetting::LambdaCm(); // 波長をセンチメートルで取得

    const auto E_photon_val = cl * h_P / lambda_cm_val; // 光子エネルギーを計算

    // ダストの種類を定義
    const auto dust_species_vec = std::vector<std::string>{"Sil", "Gra", "PAHneu", "PAHion"};

    // ダストクラスのベクトルを作成
    auto dust_vec = my_util::dust_util::MakeDustClassVector(dust_species_vec, a_cm_val);

    // Asanoモデルの初期化
    const auto asano_model = asano_model::AsanoModel(); // Asanoモデルのインスタンスを作成

    // SED構成STEP1. 星のスペクトルファイルを読み込む
    const auto stellar_spectrum = stellar_spectrum::StellarSpectrum(lambda_cm_val);

    // タイムスケールのペアに基づいて計算を行う
    for (const auto &tau_pair : tau_pair_vec) {
        free_params.SetStarFormationTimescale(tau_pair.first); // 星形成のタイムスケールを設定
        free_params.SetInfallTimescale(tau_pair.second); // 流入のタイムスケールを設定

        // Asanoモデル計算後のファイルの置き場所指定
        const auto fn_m_total = TotalDustMassFileName(free_params); // 総ダスト質量ファイル名を取得
        const auto fn_m = DustMassDistributionFileName(free_params); // ダスト質量分布ファイル名を取得
        const auto fn_n = DustNumberDistributionFileName(free_params); // ダスト数分布ファイル名を取得

        // SED構成STEP2. ダストの進化を計算
        asano_model.Calculate(free_params, fn_m_total, fn_m, fn_n, IMF_TYPE); // ダストの進化を計算

        // Asano modelの計算結果から放射伝達に使用する値を読み込む
        const auto D_val = SEDFile::ReadDustToGasMassRatio(fn_m_total, free_params.n_age_); // Asanoモデルで計算したダストとガスの質量比を読み込む
        const auto M_total_gal_val = SEDFile::ReadTotalGalaxyMassRatio(fn_m_total, free_params.n_age_); // Asanoモデルで計算した銀河の総質量比を読み込む

        // Asano modelの計算結果からダストの種類ごとの数分布を読み込む
        // n_val3[0]: silicate, n_val3[1]: graphite, n_val3[2]: neutral PAH, and n_val3[3]: ionized PAH
        const auto n_val3 = my_util::dust_util::ReadAndDivideDustNumberDistribution(fn_n, free_params.n_age_, a_cm_val);

        // ダストの特性を計算 dust extinction per unit mass, scattering albedo, and asymmetry parameter
        const auto fn_averaged_params = AveragedDustPropertyFileName(free_params); // 平均ダスト特性ファイル名を取得
        const auto dust_params_val3 = DustProperty::Calculate(a_cm_val, free_params, n_val3, fn_averaged_params); // ダスト特性を計算

        const auto Cabs_sum_val2 = SEDFile::ReadSumAbsorptionCoefficient(fn_averaged_params, free_params.n_age_); // 吸収係数を読み込む

        // 銀河の年齢に基づいて計算を行う
        for (const auto &i_age : i_age_valt) {
            free_params.SetIndexOfAge(i_age); // 年齢を設定

            stellar_spectrum.Calculate(free_params); // 星のスペクトルを計算
            const auto f_young_val = SEDFile::ReadYoungStellarFraction(FYoungFileName(free_params)); // 若い星の割合を読み込む
            const auto L_star_val = SEDFile::ReadStellarContinuumFile(StellarContinuumFileName(free_params)); // 星の連続体フラックスを読み込む

            for (const auto &M_gal : M_gal_val) {
                free_params.SetTotalGalaxyMass(M_gal); // 銀河の質量を設定

                for (auto i = std::size_t(0); i < dust_vec.size(); ++i)
                    dust_vec[i].SetNumberDistribution(M_gal, n_val3[i]); // ダストの数分布を設定

                // 年齢依存のSEDを計算
                SED_model::AgeDependent(dust_species_vec, M_gal,
                                        lambda_cm_val, a_cm_val, E_photon_val, D_val[i_age],
                                        dust_params_val3[0][i_age],
                                        dust_params_val3[1][i_age], dust_params_val3[2][i_age],
                                        dust_vec, free_params, geometry_pair_vec,
                                        L_star_val * M_gal, f_young_val, Cabs_sum_val2[i_age], n0_cnm_);
            }
        }
    }

    my_util::MyTime::PrintElapsedTime(start_time, "sec", "Total: "); // 計算にかかった時間を出力
    std::cout << std::endl;
}



// pybind11モジュールの定義
namespace py = pybind11;

// pybind11でのバインディング
PYBIND11_MODULE(sed_module, m) {
    py::class_<Hyperparameters, std::shared_ptr<Hyperparameters>>(m, "Hyperparameters")
        .def(py::init<>())
        .def_readwrite("galaxy_age", &Hyperparameters::galaxy_age)
        .def_readwrite("galaxy_mass", &Hyperparameters::galaxy_mass)
        .def_readwrite("starformation_timescale", &Hyperparameters::starformation_timescale)
        .def_readwrite("gas_infall_timescale", &Hyperparameters::gas_infall_timescale)
        .def_readwrite("dust_scale_height", &Hyperparameters::dust_scale_height)
        .def_readwrite("galaxy_radius", &Hyperparameters::galaxy_radius)
        .def_readwrite("n0_cnm", &Hyperparameters::n0_cnm)
        .def_readwrite("is_infall", &Hyperparameters::is_infall)
        .def_readwrite("imf_type", &Hyperparameters::imf_type);

    m.def("SED_calculator", &SED_calculator, "ハイパーパラメータに基づいてSEDを計算する関数");
}