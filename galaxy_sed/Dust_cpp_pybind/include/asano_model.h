#ifndef ASANO_MODEL_H
#define ASANO_MODEL_H

#include <cfloat>
#include <iostream>
#include <chrono>

#include <constants.h>
#include <file_util.h>
#include <SED_setting_util.h>
#include <calculation_util.h>

#include "function.h"
#include "read_file.h"

namespace asano_model {
using val = std::valarray<double>; // doubleの配列を表す
using val2 = std::valarray<val>; // valの配列なので、2次元配列である
using val3 = std::valarray<val2>;

/**
 int main() {
    val array1 = {1.0, 2.0, 3.0};  // double の配列
    val2 array2 = {array1, {4.0, 5.0, 6.0}};  // val の配列（2次元配列）

    // ここで array1 や array2 を使用するなどの処理が続く...

    return 0;
}
　のようにして型名を簡潔にし、配列を表すのに使う。
*/

using SEDSetting = my_util::SED_util::SEDSettingUtil;

/**
 * Calculate dust evolution with Asano model (Asano et al., 2013a, b, 2014; Nozawa et al., 2014)
 */
class AsanoModel {
  private:
    const val  a_val_         = SEDSetting::DustRadiusDustModel();
    const val3 M_SN_val3_     = ReadSNDustFile(a_val_); // Dust mass produced by supernova
    const val4 M_AGB_val4_    = ReadAGBDustFile(); // Dust mass produced by AGB
    const val3 Y_Z_val3_      = ReadMetalYieldFile(); // Metal yield from stars
    const val2 M_rem_val2_    = ReadRemnantFile(); // Stellar remnant mass
    const val4 dest_eta_val4_ = ReadSNDestructionFile(); // Dust destroyed efficiency by SN shocks
    const val  tau_star_val_  = StellarLifeTime();

  public:
    AsanoModel() = default;

    ~AsanoModel() = default;

    /**
     * Calculate dust evolution
     * @param free_params
     */
    inline void Calculate(const FreeParameter& free_params,
                          const std::string& ofn_total_dust_mass,
                          const std::string& ofn_dust_mass_distribution,
                          const std::string& ofn_dust_number_distribution,
                          const int imf_type) const {
        // 年齢の最大値から年齢のビン数を計算
        const auto n_age = static_cast<size_t>(free_params.age_max_ / TIME_BIN_DUST_MODEL) + 1;
        if (n_age > 1e6) {  
            throw std::runtime_error("n_age is too large: " + std::to_string(n_age)); // 年齢が大きすぎる場合はエラーを投げる
        }

        // 計算開始時間を記録
        const auto calc_start = std::chrono::system_clock::now();

        // 指定されたIMFタイプに基づいて正規化されたIMFを取得
        const auto IMF_val = SFH::NormalizedIMF(imf_type);

        // 年齢ごとのガス質量を格納する配列を初期化
        auto M_gas_val = val(n_age);
        if (free_params.is_infall_) M_gas_val[0] = DBL_MIN; // 初期質量を0にすると計算上問題なので、最小値を設定
        else M_gas_val[0] = 1.0; // 初期質量を1に設定

        // 年齢ごとの星の質量、金属量、炭素質量(ダスト以外も含む)、シリケート質量(ダスト以外も含む)、超新星率を格納する配列を初期化
        auto M_star_val    = val(n_age); // 星の質量
        auto M_Z_val       = val(n_age); // 金属量の質量
        auto M_C_ISM_val   = val(n_age); // ISM中の炭素質量
        auto M_sil_ISM_val = val(n_age); // ISM中のシリケート質量
        auto SNR_val       = val(n_age); // 超新星率

        // ダストの生成量を格納するための3次元配列を初期化
        val3 Y_d_val3;
        try {
            // 星によるダストの生成量(yield)
            Y_d_val3 = val3(val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES), n_age); 
        } catch (const std::bad_alloc &e) {
            throw std::runtime_error("Memory allocation failed for Y_d_val3: " + std::string(e.what())); // メモリ確保失敗時のエラー処理
        }

        // 年齢ごとの計算ループ
        for (auto i_age = std::size_t(0); i_age < n_age - 1; ++i_age) {
            const auto i_m_min      = MinimumStellarMassIndex(i_age, tau_star_val_); // 最小星質量のインデックスを取得
            const auto i_birth_valt = StellarBirthIndex(i_age, i_m_min, tau_star_val_); // 星の誕生インデックスを取得 (t - tau)に対応

            const auto i_Z_valt = StellarMetallicityIndex(i_m_min, i_birth_valt,
                                                          M_Z_val / M_gas_val); // 金属量のインデックスを取得
            const auto SFR_val  = SFH::SFRValarray(i_age, free_params.schmidt_index_,
                                                   free_params.tau_SF_, M_gas_val, i_age); // SFRを計算 // TODO:ここのageはinfallのときのageと同じものを入れたい！らしい...

            // Kanoの論文の式(4.10)に基づいて星が放出するISM質量を計算
            const auto R = ReturnMass(i_m_min, i_birth_valt, i_Z_valt, IMF_val, SFR_val,
                                      M_rem_val2_);

            // 加納の論文の式(4.12)に基づいて星がISMに出す金属生成量を計算
            const auto[Y_Z, Y_C, Y_sil] = MetalYield(i_m_min, i_birth_valt, i_Z_valt, IMF_val,
                                                     SFR_val, Y_Z_val3_);

            // 加納の論文の式(4.13)に基づいて星が出すダスト生成量を計算
            Y_d_val3[i_age + 1] = DustYield(i_m_min, i_birth_valt, M_gas_val, M_Z_val, IMF_val,
                                            SFR_val, M_SN_val3_, M_AGB_val4_);
            
            // 超新星発生率を計算
            SNR_val[i_age + 1]  = SNR(i_m_min, i_birth_valt, IMF_val, SFR_val);

            // 加納の論文の式(4.4)に基づいて銀河全体の星の質量を計算
            M_star_val[i_age + 1] = StellarMass(M_star_val[i_age], M_gas_val[i_age], R,
                                                free_params.schmidt_index_, free_params.tau_SF_, i_age);

            // 加納の論文の式(4.5)に基づいて銀河全体のISM質量を計算
            M_gas_val[i_age + 1] = GasMass(free_params.is_infall_, i_age, M_gas_val[i_age], R,
                                           free_params.schmidt_index_, free_params.tau_SF_,
                                           free_params.tau_infall_);

            // 加納の論文の式(4.6)に基づいて銀河全体の金属量を計算
            M_Z_val[i_age + 1]       = MetalMass(M_gas_val[i_age], SFR_val[i_age],
                                                 M_Z_val[i_age],
                                                 Y_Z);
            M_C_ISM_val[i_age + 1]   = MetalMass(M_gas_val[i_age], SFR_val[i_age],
                                                 M_C_ISM_val[i_age],
                                                 Y_C);
            M_sil_ISM_val[i_age + 1] = MetalMass(M_gas_val[i_age],
                                                 SFR_val[i_age], M_sil_ISM_val[i_age], Y_sil);
        }

        // 計算開始からの経過時間を表示
        my_util::MyTime::PrintElapsedTime(calc_start, "sec",
                                          "Asano model (calc before dust evolution): ");

        const auto volume_val = Map(a_val_, Volume); // ボリュームをマッピング

        // 出力ファイルストリームを作成
        auto ofp_total = MakeTotalMassFile(ofn_total_dust_mass); // 総ダスト質量ファイルを作成
        auto ofp_m     = MakeMassDistributionFile(ofn_dust_mass_distribution); // ダスト質量分布ファイルを作成
        auto ofp_n     = MakeNumberDistributionFile(ofn_dust_number_distribution); // ダスト数分布ファイルを作成

        auto M_dust_val2 = val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS); // ダスト質量を格納する2次元配列を初期化

        // 年齢ごとのループ
        for (auto i_age = std::size_t(0); i_age < n_age; ++i_age) {
            // 一定の時間ビンごとにダストサイズ分布をファイルに書き込む
            if (i_age % static_cast<int>(TIME_BIN / TIME_BIN_DUST_MODEL) == 0 && i_age != 0) {
                WriteFiles(i_age, M_Z_val[i_age], M_gas_val[i_age], M_C_ISM_val[i_age],
                           M_sil_ISM_val[i_age], M_star_val[i_age], a_val_, volume_val,
                           M_dust_val2, ofp_total, ofp_n, ofp_m);
            }
            if (i_age == n_age - 1) break; // 最後の年齢の場合はループを終了

            // ダスト質量の進化を計算
            const auto M_evolution_val2 = DustMassEvolution(M_gas_val[i_age],
                                                            {M_C_ISM_val[i_age],
                                                             M_sil_ISM_val[i_age]},
                                                            a_val_, volume_val, M_dust_val2);

            // 超新星によるダストの破壊を計算
            const auto M_SND_val2 = SNDestruction(M_gas_val[i_age], M_Z_val[i_age], SNR_val[i_age],
                                                  a_val_, M_dust_val2, dest_eta_val4_);

            // ダスト質量を計算
            M_dust_val2 = DustMass(M_gas_val[i_age], free_params.schmidt_index_,
                                   free_params.tau_SF_, M_dust_val2, Y_d_val3[i_age], M_SND_val2,
                                   M_evolution_val2, i_age);
        }

        // 3. fclose() の安全な処理
        if (ofp_total) fclose(ofp_total); // 総ダスト質量ファイルを閉じる
        if (ofp_n) fclose(ofp_n); // ダスト数分布ファイルを閉じる
        if (ofp_m) fclose(ofp_m); // ダスト質量分布ファイルを閉じる

        // 計算開始からの経過時間を表示
        my_util::MyTime::PrintElapsedTime(calc_start, "sec", "Asano model (calc): ");
    }
};


}
#endif // ASANO_MODEL_H
