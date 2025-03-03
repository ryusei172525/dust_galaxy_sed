#ifndef READ_H
#define READ_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iterator>
#include <string>
#include <valarray>

#include <constants.h>
#include <file_name.h>

namespace stellar_spectrum
{
    using SEDFile = my_util::FileUtil;

    using vec_str = std::vector<std::string>;
    using vec = std::vector<double>;
    using vec2 = std::vector<vec>;

    using val = std::valarray<double>;
    using val2 = std::valarray<val>;
    using val3 = std::valarray<val2>;

    using valt = std::valarray<std::size_t>;

    double interplinlog(double x1, double x2, double y1, double y2, double x)
    {

        double epsilon = 1.0e-37;

        if (std::fabs(x2 - x1) < epsilon)
            return y1;
        else if (y1 <= epsilon || y2 <= epsilon)
            return 0;
        else
            return y1 * std::pow(y2 / y1, (x - x1) / (x2 - x1));
    }

    double ReadMetallicity(const std::string &fn_track)
    {
        auto ifs_track = SEDFile::IfsOpen(fn_track);
        auto buff = std::string();
        std::getline(ifs_track, buff); // Skip 1 line
        auto Z_SSP = 0.0;
        ifs_track >> buff >> buff >> buff >> Z_SSP; // Skip "Metallicity", "(mass" and "fraction):"
        return Z_SSP;
    }

    vec_str ReadSSPZFileName()
    {
        auto ifs_SSPs = SEDFile::IfsOpen(SSPs_FILE_NAME);
        // read track file names
        auto fn_track_vec = std::vector<std::string>(); // fn is file name
        std::for_each(std::istream_iterator<std::string>(ifs_SSPs), {},
                      [&fn_track_vec](const auto &fn_track)
                      {
                          fn_track_vec.emplace_back(STELLAR_SPECTRUM_DATA_DIR + fn_track);
                      });
        return fn_track_vec;
    }

    vec ReadMetallicityVector(const vec_str &fn_tracks_vector)
    {
        auto Z_SSP_vec = vec(); // Z_SSP is metallicity of track file
        for (const auto &fn_track : fn_tracks_vector)
        {
            Z_SSP_vec.emplace_back(ReadMetallicity(fn_track));
        }
        std::sort(std::begin(Z_SSP_vec), std::end(Z_SSP_vec));
        return Z_SSP_vec;
    }

    double ReadFractionOfCloseBinarySystems(std::ifstream &ifs)
    {
        // Skip "Fraction of close binary systems: "
        for (auto i = std::size_t(0); i < 5; ++i)
        {
            auto buff = std::string();
            ifs >> buff;
        }
        auto fSNIa = 0.0;
        ifs >> fSNIa;
        return fSNIa;
    }

    void ReadStellibLCBcor(val &lambda_angstrom_val, val2 &F_val2)
    {
        // スペクトルデータファイルを開く
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + "stellibLCBcor.dat");

        // スペクトルの金属量の数と波長の数を格納する変数を初期化
        auto n_Z_stellib = std::size_t(0);
        auto n_lambda = std::size_t(0);
        
        // 最初の行から金属量の数と波長の数を読み取る
        SEDFile::GetLineToStringStream(ifs) >> n_Z_stellib >> n_lambda;

        auto n_spectrum = std::size_t(0); // スペクトルの総数を格納する変数を初期化

        // 各金属量に対してスペクトルの総数を読み取り、他のデータはスキップする
        for (auto i = std::size_t(0); i < n_Z_stellib; ++i)
        {
            auto n_spectrum_per_Z = std::size_t(0); // 各金属量に対するスペクトルの数
            auto buff_double = 0.0; // バッファ用の変数
            SEDFile::GetLineToStringStream(ifs) >> n_spectrum_per_Z >> buff_double; // スペクトルの数を読み取る
            for (auto j = std::size_t(0); j < n_spectrum_per_Z; ++j)
            {
                auto buff_string = std::string(); // スペクトルのデータを格納するためのバッファ
                std::getline(ifs, buff_string); // スペクトルデータをスキップ
            }
            n_spectrum += n_spectrum_per_Z; // スペクトルの総数を更新
        }

        // 波長の配列を初期 Read wavelength of continuum
        lambda_angstrom_val.resize(n_lambda);
        for (auto i = std::size_t(0); i < n_lambda; ++i)
        {
            ifs >> lambda_angstrom_val[i]; // 波長を読み取る
        }

        // スペクトルデータの配列を初期化
        F_val2.resize(n_spectrum);
        for (auto i = std::size_t(0); i < n_spectrum; ++i)
        {
            F_val2[i].resize(n_lambda); // 各スペクトルに対して波長の数を設定
            for (auto j = std::size_t(0); j < n_lambda; ++j)
            {
                ifs >> F_val2[i][j]; // スペクトルの値を読み取る
            }
        }
    }

    val2 ReadStellibCM()
    {
        // スペクトルデータファイルを開く
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + "stellibCM.dat");

        // スペクトルの数と波長の数を格納する変数を初期化
        auto nCM = std::size_t(0);
        auto n_lambda = std::size_t(0);
        
        // 最初の行からスペクトルの数と波長の数を読み取る
        SEDFile::GetLineToStringStream(ifs) >> nCM >> n_lambda;

        // nCM行分のデータをスキップする
        for (auto i = std::size_t(0); i < nCM; ++i)
        {
            auto buff_string = std::string(); // バッファを初期化
            std::getline(ifs, buff_string); // スペクトルデータをスキップ
        }

        // 波長のデータをスキップする
        for (auto i = std::size_t(0); i < n_lambda; ++i)
        {
            auto buff = 0.0; // バッファ用の変数
            ifs >> buff; // 波長の値を読み取るが、使用しないためスキップ
        }

        // スペクトルデータを格納するための2次元配列を初期化
        auto F_val2 = val2(val(n_lambda), nCM);
        for (auto i = std::size_t(0); i < nCM; ++i)
        {
            for (auto j = std::size_t(0); j < n_lambda; ++j)
            {
                ifs >> F_val2[i][j]; // スペクトルの値を読み取って配列に格納
            }
        }
        return F_val2; // 読み取ったスペクトルデータを返す
    }

    /**
     * Read stellar library
     * @param lambda_angstrom_val Wavelength [A]
     * @param F_stellar_LCB_val2 Terr < 50,000 K stellar library of Lejeune et al. (1997, 1998)
     * @param F_stellar_CM_val2  Terr > 50,000 K stellar library of Clegg & Middlemass (1987)
     */
    void ReadStellarSpectrumFile(val &lambda_angstrom_val, val2 &F_stellar_LCB_val2,
                                 val2 &F_stellar_CM_val2)
    {
        // "stellibLCBcor.dat": stellar library of Lejeune et al. (1997, 1998) Terr < 50,000 K
        ReadStellibLCBcor(lambda_angstrom_val, F_stellar_LCB_val2);
        
        // "stellibCM.dat": stellar library of Clegg & Middlemass (1987); Terr > 50,000 K
        F_stellar_CM_val2 = ReadStellibCM();
    }

    // C++ には split がなくて書き直すのが面倒なので C style で書いてあるまま
    void ReadNebularSpectrumFile(const val &lambda_angstrom_val, val &F_neb_val, val &lambda_line_val,
                                 val &F_line_val)
    {
        // 各種変数を初期化
        auto g1 = 0.0; // スペクトルの成分1
        auto g2 = 0.0; // スペクトルの成分2
        auto g3 = 0.0; // スペクトルの成分3
        auto g4 = 0.0; // スペクトルの成分4

        // HII.datファイルを開く
        auto ifp = SEDFile::IfpOpen(STELLAR_SPECTRUM_DATA_DIR + "HII.dat");

        // ヘリウムの状態を格納する変数を初期化
        auto HeI = 0.0; // HeIの値
        auto HeII = 0.0; // HeIIの値
        auto alphaTe = 0.0; // 電子温度の値

        // fscanfを使用してHeI, HeII, alphaTeを読み取る
        if (fscanf(ifp, "%lf %lf %lf", &HeI, &HeII, &alphaTe) != 3)
        {
            fprintf(stderr, "Error reading HeI, HeII, and alphaTe\n");
            return; // 読み取りエラーが発生した場合は関数を終了
        }

        // 波長の数を格納する変数を初期化
        auto n_lambda_continuum = std::size_t(0);

        // fscanfを使用して波長の数を読み取る
        if (fscanf(ifp, "%lu", &n_lambda_continuum) != 1)
        {
            fprintf(stderr, "Error reading n_lambda_continuum\n");
            return; // 読み取りエラーが発生した場合は関数を終了
        }

        char hoge[256];
        fgets(hoge, sizeof(hoge), ifp); // ヘッダー行をスキップ

        // 波長の配列を初期化
        auto lambda_continuum = val(n_lambda_continuum); // HII.datの連続波長
        auto gam1 = val(100); // スペクトル成分1
        auto gam2 = val(100); // スペクトル成分2
        auto gam3 = val(100); // スペクトル成分3
        auto gam4 = val(100); // スペクトル成分4

        // 各波長に対してデータを読み取る
        for (auto i = std::size_t(0); i < n_lambda_continuum; ++i)
        {
            char gam[256]; // スペクトル成分を格納するためのバッファ
            fgets(gam, sizeof(gam), ifp); // 行を読み取る
            sscanf(gam, " %lf %lf %lf %lf %lf", &lambda_continuum[i],
                   &gam1[i], &gam2[i], &gam3[i], &gam4[i]); // 各成分を読み取る
        }

        // スペクトルのフラックスを波長に合わせて補間する
        auto j = std::size_t(0);
        for (auto i = std::size_t(0); i < lambda_angstrom_val.size(); i++)
        {
            // lambdacontとlambdaの大小関係が入れ替わるところを探す
            while ((lambda_continuum[j] - lambda_angstrom_val[i]) *
                           (lambda_continuum[j + 1] - lambda_angstrom_val[i]) >
                       0 &&
                   j < n_lambda_continuum)
            {
                j++; // jをインクリメント
            }

            // 補間を行う
            if ((lambda_continuum[j] - lambda_angstrom_val[i]) *
                    (lambda_continuum[j + 1] - lambda_angstrom_val[i]) <=
                0)
            {
                // 各成分の補間を行う
                g1 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam1[j], gam1[j + 1], 1.0 / lambda_angstrom_val[i]);
                g2 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam2[j], gam2[j + 1], 1.0 / lambda_angstrom_val[i]);
                g3 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam3[j], gam3[j + 1], 1.0 / lambda_angstrom_val[i]);
                g4 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam4[j], gam4[j + 1], 1.0 / lambda_angstrom_val[i]);
            }

            // 波長が連続波長の最大値を超えた場合の処理
            if (lambda_angstrom_val[i] >= lambda_continuum[n_lambda_continuum - 1])
            {
                g1 = gam1[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
                g2 = gam2[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
                g3 = gam3[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
                g4 = gam4[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
            }

            // フラックスを計算し、結果をF_neb_valに格納
            F_neb_val[i] = 1e-40 * (g1 + g2 + HeI * g3 + HeII * g4) * cl * 1e8 / alphaTe / (lambda_angstrom_val[i] * lambda_angstrom_val[i]);
        }

        // 星雲線の計算
        auto n_line = std::size_t(0);
        if (fscanf(ifp, "%lu", &n_line) != 1)
        {
            fprintf(stderr, "Error reading n_line\n");
            return; // 読み取りエラーが発生した場合は関数を終了
        }
        fgets(hoge, sizeof(hoge), ifp); // ヘッダー行をスキップ
        auto gline = val(n_line); // 星雲線のフラックスを格納する配列
        lambda_line_val.resize(n_line); // 波長の配列を初期化
        F_line_val.resize(n_line); // 星雲線のフラックスの配列を初期化

        // 各星雲線のデータを読み取る
        for (auto i = std::size_t(0); i < n_line; ++i)
        {
            char gam[256]; // スペクトル成分を格納するためのバッファ
            fgets(gam, sizeof(gam), ifp); // 行を読み取る
            sscanf(gam, "%lf %lf", &lambda_line_val[i], &gline[i]); // 波長とフラックスを読み取る
            F_line_val[i] = 1.244e-25 * gline[i] / alphaTe; // フラックスを計算
        }
        fclose(ifp); // ファイルを閉じる
    }

    // Table 5.4 of Spitzer (1978, p.113)
    // tau_dust_Spitzer_varray is defined as the optical thickness of the dust just shortward of the H ionization limit
    void
    ReadSpitzer(std::size_t &n_spitzer, val &tau_dust_Spitzer_varray, val &y_Spitzer_varray)
    {
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + "Spitzer.dat");
        ifs >> n_spitzer;
        for (auto i = std::size_t(0); i < n_spitzer; ++i)
        {
            ifs >> tau_dust_Spitzer_varray[i] >> y_Spitzer_varray[i];
        }
    }

    /**
     * 星のスペクトルにおける瞬時のバーストから放出されるフラックスを読み取る関数
     * @param fn_track_vec トラックファイル名のリスト
     * @param F_stellar_LCB_val2 50,000 K未満の星のライブラリ
     * @param F_stellar_CM_val2 50,000 K以上の星のライブラリ
     * @param t_SSP_val 時間の配列
     * @param F_bol_SSP_val2 ボルフラックスの配列
     * @param F_SSP_val3 スペクトルフラックスの配列
     * @param m_alive_val2 生存している星の質量の配列
     * @param t_inv_val 時間の逆インデックス
     * @param beta_val 補間係数の配列
     * @param n_Lym_val2 Lymanフラックスの配列
     * @param n_SNII_val2 II型超新星の数の配列
     * @param n_SNIa_val2 Ia型超新星の数の配列
     * @param ejecta_val2 放出物質の質量の配列
     * @param ejecta_Z_val2 放出物質の金属量の配列
     * @param m_BHNS_val2 ブラックホール中性子星の質量の配列
     * @param m_WD_val2 白色矮星の質量の配列
     */
    void ReadStellarFluxFile(const vec_str &fn_track_vec, const val2 &F_stellar_LCB_val2,
                             const val2 &F_stellar_CM_val2, valt &t_SSP_val, val2 &F_bol_SSP_val2,
                             val3 &F_SSP_val3, val2 &m_alive_val2, valt &t_inv_val, val &beta_val,
                             val2 &n_Lym_val2, val2 &n_SNII_val2, val2 &n_SNIa_val2, val2 &ejecta_val2,
                             val2 &ejecta_Z_val2, val2 &m_BHNS_val2, val2 &m_WD_val2)
    {
        auto n_time_SSP = std::size_t(0); // SSPの時間数を格納する変数を初期化

        // トラックファイルを読み込むループ
        for (auto i_Z = std::size_t(0); i_Z < fn_track_vec.size(); ++i_Z)
        {
            // トラックファイルを開く
            auto ifs = SEDFile::IfsOpen(fn_track_vec[i_Z]);

            // ヘッダーを読み飛ばす
            for (auto i = std::size_t(0); i < 4; ++i)
            {
                // getlineは空白のみの行も前の行と一緒に読み込むので見た目は5必要だが4で良い
                auto buff = std::string();
                std::getline(ifs, buff); // ヘッダー行をスキップ
            }

            // 最初の行から時間の数を読み取る
            SEDFile::GetLineToStringStream(ifs) >> n_time_SSP;

            auto n_used = std::size_t(0); // 使用されるスペクトルの数
            auto x = 0.0; // 一時的な変数
            double n_HI_tot_pre; // 中性HIの密度を格納する変数
            // 時間、使用されるスペクトルの数、ボルフラックス、前の中性HIの密度を読み取る
            SEDFile::GetLineToStringStream(ifs) >> t_SSP_val[0] >> n_used >> x >> F_bol_SSP_val2[i_Z][0] >> n_HI_tot_pre;
            
            // 他の行をスキップ
            SEDFile::SkipLine(ifs);
            SEDFile::SkipLine(ifs);

            // 使用されるスペクトルの数だけループ
            for (auto i = std::size_t(0); i < n_used; ++i)
            {
                auto intgr = 0; // スペクトルのインデックス
                auto F_spec = 0.0; // スペクトルのフラックス
                ifs >> intgr >> F_spec; // スペクトルのインデックスとフラックスを読み取る
                
                // インデックスが正の場合
                if (intgr > 0)
                {
                    const auto intgr_s = static_cast<std::size_t>(intgr + 1); // インデックスを調整
                    if (intgr_s > F_stellar_LCB_val2.size() - 1)
                        continue; // インデックスが範囲外の場合はスキップ
                    F_SSP_val3[i_Z][0] += F_spec * F_stellar_LCB_val2[intgr_s]; // フラックスを加算
                }
                // インデックスが負の場合
                else
                {
                    const auto intgr_s = static_cast<std::size_t>(-intgr - 1); // インデックスを調整
                    if (intgr_s > F_stellar_CM_val2.size() - 1)
                        continue; // インデックスが範囲外の場合はスキップ
                    F_SSP_val3[i_Z][0] += F_spec * F_stellar_CM_val2[intgr_s]; // フラックスを加算
                }
            }

            m_alive_val2[i_Z][0] = 1.0; // 生存している星の質量を初期化

            // t_SSPが1でない場合の処理
            for (auto i_t_SSP = std::size_t(1); i_t_SSP < n_time_SSP; ++i_t_SSP)
            {
                auto n_HI_tot = 0.0; // 中性HIの総密度を初期化
                // 時間、使用されるスペクトルの数、ボルフラックス、中性HIの密度を読み取る
                ifs >> t_SSP_val[i_t_SSP] >> n_used >> x >> F_bol_SSP_val2[i_Z][i_t_SSP] >> n_HI_tot;
                SEDFile::SkipLine(ifs); // 他の行をスキップ

                // 各種変数を初期化
                auto n_SNII_IS = 0.0; // II型超新星の数
                auto n_SNII_CB = 0.0; // II型超新星の数（CB）
                auto n_SNIa = 0.0; // Ia型超新星の数
                auto m_BHNS_IS = 0.0; // ブラックホール中性子星の質量（IS）
                auto m_BHNS_CB = 0.0; // ブラックホール中性子星の質量（CB）
                auto m_WD_IS = 0.0; // 白色矮星の質量（IS）
                auto m_WD_CB = 0.0; // 白色矮星の質量（CB）
                auto ejecta_IS = 0.0; // 放出物質の質量（IS）
                auto ejecta_CB = 0.0; // 放出物質の質量（CB）
                auto ejecta_Z_IS = 0.0; // 放出物質の金属量（IS）
                auto ejecta_Z_CB = 0.0; // 放出物質の金属量（CB）

                // 各種データを読み取る
                ifs >> n_SNII_IS >> n_SNII_CB >> n_SNIa >> m_BHNS_IS >> m_BHNS_CB >> m_WD_IS >> m_WD_CB >> ejecta_IS >> ejecta_CB >> ejecta_Z_IS >> ejecta_Z_CB;

                // Ia型超新星の数を調整
                n_SNIa *= F_SNIa;
                // II型超新星の数を計算
                const auto n_SNII = (1.0 - F_SNIa) * n_SNII_IS + F_SNIa * n_SNII_CB;
                // ブラックホール中性子星の質量を計算
                const auto m_BHNS = (1.0 - F_SNIa) * m_BHNS_IS + F_SNIa * m_BHNS_CB;
                // 白色矮星の質量を計算
                const auto m_WD = (1.0 - F_SNIa) * m_WD_IS + F_SNIa * m_WD_CB;
                // 放出物質の質量を計算
                const auto ejecta = (1.0 - F_SNIa) * ejecta_IS + F_SNIa * ejecta_CB;
                // 放出物質の金属量を計算
                const auto ejecta_Z = (1.0 - F_SNIa) * ejecta_Z_IS + F_SNIa * ejecta_Z_CB;

                // 使用されるスペクトルの数だけループ
                for (auto i = std::size_t(0); i < n_used; ++i)
                {
                    auto intgr = 0; // スペクトルのインデックス
                    auto F_spec = 0.0; // スペクトルのフラックス
                    ifs >> intgr >> F_spec; // スペクトルのインデックスとフラックスを読み取る
                    
                    // インデックスが正の場合
                    if (intgr > 0)
                    {
                        const auto intgr_s = static_cast<std::size_t>(intgr); // インデックスを調整
                        if (intgr_s > F_stellar_LCB_val2.size() - 1)
                            continue; // インデックスが範囲外の場合はスキップ
                        F_SSP_val3[i_Z][i_t_SSP] += F_spec * F_stellar_LCB_val2[intgr_s]; // フラックスを加算
                    }
                    // インデックスが負の場合
                    else
                    {
                        const auto intgr_s = static_cast<std::size_t>(-intgr); // インデックスを調整
                        if (intgr_s > F_stellar_CM_val2.size() - 1)
                            continue; // インデックスが範囲外の場合はスキップ
                        F_SSP_val3[i_Z][i_t_SSP] += F_spec * F_stellar_CM_val2[intgr_s]; // フラックスを加算
                    }
                }

                // t_SSP_val[i_t_SSP-1]からt_SSP/val[i_t_SSP]までの間を補間する
                const auto dt = static_cast<double>(t_SSP_val[i_t_SSP] - t_SSP_val[i_t_SSP - 1]);

                // 補間処理
                for (auto t = t_SSP_val[i_t_SSP - 1]; t < t_SSP_val[i_t_SSP]; ++t)
                {
                    // t_SSP_vals 時刻ビンの左側の値をつめる
                    t_inv_val[t] = i_t_SSP - 1; // 時間の逆インデックスを設定
                    beta_val[t] = static_cast<double>(t_SSP_val[i_t_SSP] - t) / dt; // 補間係数を計算
                    
                    // Lymanフラックスの補間
                    n_Lym_val2[i_Z][t] = n_HI_tot_pre + (n_HI_tot - n_HI_tot_pre) *
                                                            static_cast<double>(t - t_SSP_val[i_t_SSP - 1]) / dt;

                    // これらの値は補間しないでそのまま右側の値を使う
                    n_SNII_val2[i_Z][t] = n_SNII; // II型超新星の数
                    n_SNIa_val2[i_Z][t] = n_SNIa; // Ia型超新星の数
                    ejecta_val2[i_Z][t] = ejecta; // 放出物質の質量
                    ejecta_Z_val2[i_Z][t] = ejecta_Z; // 放出物質の金属量
                    m_BHNS_val2[i_Z][t] = m_BHNS; // ブラックホール中性子星の質量
                    m_WD_val2[i_Z][t] = m_WD; // 白色矮星の質量

                    // 生存している星の質量を計算
                    m_alive_val2[i_Z][t + 1] = m_alive_val2[i_Z][t] - ejecta_val2[i_Z][t] - m_BHNS_val2[i_Z][t] - m_WD_val2[i_Z][t];
                }
                n_HI_tot_pre = n_HI_tot; // 前の中性HIの密度を更新
            }
        }
        
        // 後ろの境界条件を設定
        t_inv_val[t_SSP_val[n_time_SSP - 1]] = n_time_SSP - 1 - 1; // 最後の時間ビンの逆インデックスを設定
        beta_val[t_SSP_val[n_time_SSP - 1]] = 0.0; // 最後の時間ビンの補間係数を0に設定
    }

    void ReadSFRFile(int typeSFR, const std::string &fileSFR, std::size_t &n_time_file,
                     val &time_file_varray, val &SFR_file_varray, val &Z_file_varray)
    {
        if (typeSFR == -1)
        {
            auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + fileSFR);
            n_time_file = 0;
            while (ifs >> time_file_varray[n_time_file] >> SFR_file_varray[n_time_file])
            {
                ++n_time_file;
            }
            return;
        }

        // typeSFR == -2
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + fileSFR);
        n_time_file = 0;
        while (ifs >> time_file_varray[n_time_file] >> SFR_file_varray[n_time_file] >> Z_file_varray[n_time_file])
        {
            ++n_time_file;
        }
    }

    void ReadScenariosFile(std::ifstream &ifs, double &Z_gas, bool &is_infall, double &t_infall,
                           double &Z_infall, int &type_SFR, double &SFR_param1, double &SFR_param2,
                           val &SFR_param_varray, std::string &file_name_SFR, bool &is_code_Z,
                           double &Z_SFR, double &f_sub, double &t_wind, std::string &answer_neb,
                           int &code_ext)
    {
        SEDFile::SkipLine(ifs);
        SEDFile::SkipLine(ifs);
        SEDFile::SkipLine(ifs);

        // Read Initial metallicity
        SEDFile::GetLineToStringStream(ifs, ":") >> Z_gas;

        auto answer_infall = SEDFile::GetLineByString(ifs);

        if (answer_infall == "Infall")
        {
            is_infall = true;
            SEDFile::GetLineToStringStream(ifs, ":") >> t_infall;
            SEDFile::GetLineToStringStream(ifs, ":") >> Z_infall;
        }

        SEDFile::GetLineToStringStream(ifs, ":") >> type_SFR;
        if (type_SFR >= 1 && type_SFR <= 3)
        {
            SEDFile::GetLineToStringStream(ifs, ":") >> SFR_param1;
            SEDFile::GetLineToStringStream(ifs, ":") >> SFR_param2;
        }
        if (type_SFR >= 10)
        {
            auto n_param = std::size_t(0);
            SEDFile::GetLineToStringStream(ifs, ":") >> n_param;
            for (auto i = std::size_t(0); i < n_param; ++i)
            {
                SEDFile::GetLineToStringStream(ifs, ":") >> SFR_param_varray[i];
            }
        }

        if (type_SFR <= -1)
        {
            SEDFile::GetLineToStringStream(ifs, ":") >> file_name_SFR;
        }

        if (type_SFR >= -1)
        {
            const auto answer_Z = SEDFile::GetLineByString(ifs);
            if (answer_Z.find("No") != std::string::npos)
            { // When answer_Z contains "No"
                is_code_Z = false;
                SEDFile::GetLineToStringStream(ifs, ":") >> Z_SFR;
            }
            else
            {
                Z_SFR = Z_gas;
            }
        }

        SEDFile::GetLineToStringStream(ifs, ":") >> f_sub;
        auto answer_wind = SEDFile::GetLineByString(ifs);

        // When answer_wind contains "Galactic winds"
        if (answer_wind.find("Galactic winds\n") != std::string::npos)
        {
            SEDFile::GetLineToStringStream(ifs, ":") >> t_wind;
        }

        answer_neb = SEDFile::GetLineByString(ifs);
        auto answer_ext = SEDFile::GetLineByString(ifs);

        if (answer_ext.find("No") != std::string::npos)
            code_ext = 0;
        if (answer_ext.find("spheroidal") != std::string::npos)
            code_ext = 1;
        if (answer_ext.find("averaged") != std::string::npos)
            code_ext = 2;
        if (answer_ext.find("specific") != std::string::npos)
        {
            code_ext = 3;
            SEDFile::SkipLine(ifs);
        }
    }

    val MoveAverage(val &input)
    {
        auto output = input;
        for (auto j = std::size_t(588) + 5; j < N_LAMBDA - 5; ++j)
        {
            output[j] = (input[j - 5] + input[j - 4] + input[j - 3] + input[j - 2] + input[j - 1] + input[j] + input[j + 5] + input[j + 4] + input[j + 3] + input[j + 2] +
                         input[j + 1]) /
                        11;
        }
        return output;
    }
}
#endif
