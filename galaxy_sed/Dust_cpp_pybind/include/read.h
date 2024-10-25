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

    void
    ReadStellibLCBcor(val &lambda_angstrom_val, val2 &F_val2)
    {
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + "stellibLCBcor.dat");

        auto n_Z_stellib = std::size_t(0);
        auto n_lambda = std::size_t(0);
        SEDFile::GetLineToStringStream(ifs) >> n_Z_stellib >> n_lambda;

        auto n_spectrum = std::size_t(0); // The number of spectrum

        // Read only the total number of spectrum, skip other data
        for (auto i = std::size_t(0); i < n_Z_stellib; ++i)
        {
            auto n_spectrum_per_Z = std::size_t(0);
            auto buff_double = 0.0;
            SEDFile::GetLineToStringStream(ifs) >> n_spectrum_per_Z >> buff_double;
            for (auto j = std::size_t(0); j < n_spectrum_per_Z; ++j)
            {
                auto buff_string = std::string();
                std::getline(ifs, buff_string);
            }
            n_spectrum += n_spectrum_per_Z;
        }

        // Read wavelength of continuum
        lambda_angstrom_val.resize(n_lambda);
        for (auto i = std::size_t(0); i < n_lambda; ++i)
        {
            ifs >> lambda_angstrom_val[i];
        }

        F_val2.resize(n_spectrum);
        for (auto i = std::size_t(0); i < n_spectrum; ++i)
        {
            F_val2[i].resize(n_lambda);
            for (auto j = std::size_t(0); j < n_lambda; ++j)
            {
                ifs >> F_val2[i][j];
            }
        }
    }

    val2 ReadStellibCM()
    {
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + "stellibCM.dat");

        auto nCM = std::size_t(0);
        auto n_lambda = std::size_t(0);
        SEDFile::GetLineToStringStream(ifs) >> nCM >> n_lambda;

        for (auto i = std::size_t(0); i < nCM; ++i)
        {
            auto buff_string = std::string();
            std::getline(ifs, buff_string); // Skip nCM lines
        }

        // skip lambdas
        for (auto i = std::size_t(0); i < n_lambda; ++i)
        {
            auto buff = 0.0;
            ifs >> buff;
        }

        auto F_val2 = val2(val(n_lambda), nCM);
        for (auto i = std::size_t(0); i < nCM; ++i)
        {
            for (auto j = std::size_t(0); j < n_lambda; ++j)
            {
                ifs >> F_val2[i][j];
            }
        }
        return F_val2;
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
        // "stellibCM.dat": stellar library of Lejeune et al. (1997, 1998) Terr < 50,000 K
        ReadStellibLCBcor(lambda_angstrom_val, F_stellar_LCB_val2);
        // "stellibLCBcor.dat": stellar library of Clegg & Middlemass (1987); Terr > 50,000 K
        F_stellar_CM_val2 = ReadStellibCM();
    }

    // C++ には split がなくて書き直すのが面倒なので C style で書いてあるまま
    void ReadNebularSpectrumFile(const val &lambda_angstrom_val, val &F_neb_val, val &lambda_line_val,
                                 val &F_line_val)
    {
        auto g1 = 0.0;
        auto g2 = 0.0;
        auto g3 = 0.0;
        auto g4 = 0.0;

        auto ifp = SEDFile::IfpOpen(STELLAR_SPECTRUM_DATA_DIR + "HII.dat");

        auto HeI = 0.0;
        auto HeII = 0.0;
        auto alphaTe = 0.0;

        // fscanfのエラーハンドリング
        if (fscanf(ifp, "%lf %lf %lf", &HeI, &HeII, &alphaTe) != 3)
        {
            fprintf(stderr, "Error reading HeI, HeII, and alphaTe\n");
            return;
        }

        auto n_lambda_continuum = std::size_t(0);

        // fscanfのエラーハンドリング
        if (fscanf(ifp, "%lu", &n_lambda_continuum) != 1)
        {
            fprintf(stderr, "Error reading n_lambda_continuum\n");
            return;
        }

        char hoge[256];
        fgets(hoge, sizeof(hoge), ifp); // Skip line

        auto lambda_continuum = val(n_lambda_continuum); // Wavelength of continuum in HII.dat
        auto gam1 = val(100);
        auto gam2 = val(100);
        auto gam3 = val(100);
        auto gam4 = val(100);

        for (auto i = std::size_t(0); i < n_lambda_continuum; ++i)
        {
            char gam[256]; // gamma
            fgets(gam, sizeof(gam), ifp);
            sscanf(gam, " %lf %lf %lf %lf %lf", &lambda_continuum[i],
                   &gam1[i], &gam2[i], &gam3[i], &gam4[i]);
        }

        // interpolate flax for adjusting to stellar library wavelength step
        auto j = std::size_t(0);
        for (auto i = std::size_t(0); i < lambda_angstrom_val.size(); i++)
        {
            // lambdacont と lambda の大小関係が入れ替わるところを探す
            while ((lambda_continuum[j] - lambda_angstrom_val[i]) *
                           (lambda_continuum[j + 1] - lambda_angstrom_val[i]) >
                       0 &&
                   j < n_lambda_continuum)
            {
                j++;
            }

            if ((lambda_continuum[j] - lambda_angstrom_val[i]) *
                    (lambda_continuum[j + 1] - lambda_angstrom_val[i]) <=
                0)
            {
                g1 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam1[j], gam1[j + 1], 1.0 / lambda_angstrom_val[i]);
                g2 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam2[j], gam2[j + 1], 1.0 / lambda_angstrom_val[i]);
                g3 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam3[j], gam3[j + 1], 1.0 / lambda_angstrom_val[i]);
                g4 = interplinlog(1.0 / lambda_continuum[j], 1.0 / lambda_continuum[j + 1],
                                  gam4[j], gam4[j + 1], 1.0 / lambda_angstrom_val[i]);
            }
            if (lambda_angstrom_val[i] >= lambda_continuum[n_lambda_continuum - 1])
            {
                g1 = gam1[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
                g2 = gam2[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
                g3 = gam3[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
                g4 = gam4[n_lambda_continuum - 1] * pow((lambda_angstrom_val[i] / lambda_continuum[n_lambda_continuum - 1]), 0.1);
            }
            // cl is units in [cm/s] and convert to [um/s]
            F_neb_val[i] = 1e-40 * (g1 + g2 + HeI * g3 + HeII * g4) * cl * 1e8 / alphaTe / (lambda_angstrom_val[i] * lambda_angstrom_val[i]);
        }

        // Calculations for the nebular lines
        auto n_line = std::size_t(0);
        if (fscanf(ifp, "%lu", &n_line) != 1)
        {
            fprintf(stderr, "Error reading n_line\n");
            return;
        }
        fgets(hoge, sizeof(hoge), ifp);
        auto gline = val(n_line);
        lambda_line_val.resize(n_line);
        F_line_val.resize(n_line);

        for (auto i = std::size_t(0); i < n_line; ++i)
        {
            char gam[256]; // gamma
            fgets(gam, sizeof(gam), ifp);
            sscanf(gam, "%lf %lf", &lambda_line_val[i], &gline[i]);
            F_line_val[i] = 1.244e-25 * gline[i] / alphaTe;
        }
        fclose(ifp);
    }

    // 現在は使っていないので中身を把握していない
    void ReadDustFile(std::size_t &n_lambda_ext, val &Z_ext_varray, val2 &coeff_ext_varray2,
                      val &lambda_ext_varray, val2 &tau_ext_varray2, val2 &albedo_ext_varray2,
                      val2 &asym_ext_varray2)
    {
        auto ifs = SEDFile::IfsOpen(STELLAR_SPECTRUM_DATA_DIR + "dust.dat");
        // Skip 2 lines
        auto buff = std::string();
        std::getline(ifs, buff);
        std::getline(ifs, buff);

        for (auto i = std::size_t(1); i < 4; ++i)
        {
            ifs >> Z_ext_varray[i] >> coeff_ext_varray2[i][0] >> coeff_ext_varray2[i][1];
            std::getline(ifs, buff);
        }

        Z_ext_varray[0] = 0;
        coeff_ext_varray2[0][0] = coeff_ext_varray2[1][0];
        coeff_ext_varray2[0][1] = coeff_ext_varray2[1][1];
        Z_ext_varray[4] = 1;
        coeff_ext_varray2[4][0] = coeff_ext_varray2[3][0];
        coeff_ext_varray2[4][1] = coeff_ext_varray2[3][1];

        std::getline(ifs, buff); // Skip line
        ifs >> n_lambda_ext;
        std::getline(ifs, buff); // Skip line

        for (auto i = std::size_t(0); i < n_lambda_ext; ++i)
        {
            ifs >> lambda_ext_varray[i];
            for (auto j = std::size_t(0); j < 2; ++j)
            {
                ifs >> tau_ext_varray2[i][j] >> albedo_ext_varray2[i][j] >> asym_ext_varray2[i][j];
            }
        }
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
     * Reading of the flux emitted for an instantaneous burst in the stellar spectra
     * @param fn_track_vec Track file name list
     * @param F_stellar_LCB_val2 Stellar library
     * @param F_stellar_CM_val2 Stellar library
     * @param t_SSP_val
     * @param F_bol_SSP_val2
     * @param F_SSP_val3
     * @param m_alive_val2
     * @param t_inv_val
     * @param beta_val
     * @param n_Lym_val2
     * @param n_SNII_val2
     * @param n_SNIa_val2
     * @param ejecta_val2
     * @param ejecta_Z_val2
     * @param m_BHNS_val2
     * @param m_WD_val2
     */
    void ReadStellarFluxFile(const vec_str &fn_track_vec, const val2 &F_stellar_LCB_val2,
                             const val2 &F_stellar_CM_val2, valt &t_SSP_val, val2 &F_bol_SSP_val2,
                             val3 &F_SSP_val3, val2 &m_alive_val2, valt &t_inv_val, val &beta_val,
                             val2 &n_Lym_val2, val2 &n_SNII_val2, val2 &n_SNIa_val2, val2 &ejecta_val2,
                             val2 &ejecta_Z_val2, val2 &m_BHNS_val2, val2 &m_WD_val2)
    {
        auto n_time_SSP = std::size_t(0);

        // Read track file
        for (auto i_Z = std::size_t(0); i_Z < fn_track_vec.size(); ++i_Z)
        {
            auto ifs = SEDFile::IfsOpen(fn_track_vec[i_Z]);

            // Read headers
            for (auto i = std::size_t(0); i < 4; ++i)
            {
                // getline は空白のみの行も前の行と一緒に読み込むので見た目は5必要だが4で良い
                auto buff = std::string();
                std::getline(ifs, buff);
            }

            SEDFile::GetLineToStringStream(ifs) >> n_time_SSP;

            auto n_used = std::size_t(0);
            auto x = 0.0;
            double n_HI_tot_pre; // Density of neutral HI
            SEDFile::GetLineToStringStream(ifs) >> t_SSP_val[0] >> n_used >> x >> F_bol_SSP_val2[i_Z][0] >> n_HI_tot_pre;
            // skip others
            SEDFile::SkipLine(ifs);
            SEDFile::SkipLine(ifs);

            for (auto i = std::size_t(0); i < n_used; ++i)
            {
                auto intgr = 0; // 負の値もあるので std::size_t ではなく int
                auto F_spec = 0.0;
                ifs >> intgr >> F_spec;
                if (intgr > 0)
                {
                    const auto intgr_s = static_cast<std::size_t>(intgr + 1);
                    if (intgr_s > F_stellar_LCB_val2.size() - 1)
                        continue;
                    F_SSP_val3[i_Z][0] += F_spec * F_stellar_LCB_val2[intgr_s];
                }
                else
                {
                    const auto intgr_s = static_cast<std::size_t>(-intgr - 1);
                    if (intgr_s > F_stellar_CM_val2.size() - 1)
                        continue;
                    F_SSP_val3[i_Z][0] += F_spec * F_stellar_CM_val2[intgr_s];
                }
            }

            m_alive_val2[i_Z][0] = 1.0;

            // t_SSP != 1
            for (auto i_t_SSP = std::size_t(1); i_t_SSP < n_time_SSP; ++i_t_SSP)
            {
                auto n_HI_tot = 0.0;
                ifs >> t_SSP_val[i_t_SSP] >> n_used >> x >> F_bol_SSP_val2[i_Z][i_t_SSP] >> n_HI_tot;
                SEDFile::SkipLine(ifs);

                auto n_SNII_IS = 0.0;
                auto n_SNII_CB = 0.0;
                auto n_SNIa = 0.0;
                auto m_BHNS_IS = 0.0;
                auto m_BHNS_CB = 0.0;
                auto m_WD_IS = 0.0;
                auto m_WD_CB = 0.0;
                auto ejecta_IS = 0.0;
                auto ejecta_CB = 0.0;
                auto ejecta_Z_IS = 0.0;
                auto ejecta_Z_CB = 0.0;
                ifs >> n_SNII_IS >> n_SNII_CB >> n_SNIa >> m_BHNS_IS >> m_BHNS_CB >> m_WD_IS >> m_WD_CB >> ejecta_IS >> ejecta_CB >> ejecta_Z_IS >> ejecta_Z_CB;

                n_SNIa *= F_SNIa;
                const auto n_SNII = (1.0 - F_SNIa) * n_SNII_IS + F_SNIa * n_SNII_CB;
                const auto m_BHNS = (1.0 - F_SNIa) * m_BHNS_IS + F_SNIa * m_BHNS_CB;
                const auto m_WD = (1.0 - F_SNIa) * m_WD_IS + F_SNIa * m_WD_CB;
                const auto ejecta = (1.0 - F_SNIa) * ejecta_IS + F_SNIa * ejecta_CB;
                const auto ejecta_Z = (1.0 - F_SNIa) * ejecta_Z_IS + F_SNIa * ejecta_Z_CB;

                for (auto i = std::size_t(0); i < n_used; ++i)
                {
                    auto intgr = 0;
                    auto F_spec = 0.0;
                    ifs >> intgr >> F_spec;
                    if (intgr > 0)
                    {
                        const auto intgr_s = static_cast<std::size_t>(intgr);
                        if (intgr_s > F_stellar_LCB_val2.size() - 1)
                            continue;
                        F_SSP_val3[i_Z][i_t_SSP] += F_spec * F_stellar_LCB_val2[intgr_s];
                    }
                    else
                    {
                        const auto intgr_s = static_cast<std::size_t>(-intgr);
                        if (intgr_s > F_stellar_CM_val2.size() - 1)
                            continue;
                        F_SSP_val3[i_Z][i_t_SSP] += F_spec * F_stellar_CM_val2[intgr_s];
                    }
                }

                // t_SSP_val[i_t_SSP-1] から t_SSP/val[i_t_SSP] までの間を補間する
                const auto dt = static_cast<double>(t_SSP_val[i_t_SSP] - t_SSP_val[i_t_SSP - 1]);

                for (auto t = t_SSP_val[i_t_SSP - 1]; t < t_SSP_val[i_t_SSP]; ++t)
                {
                    // t_SSP_vals 時刻ビンの左側の値をつめる
                    t_inv_val[t] = i_t_SSP - 1;
                    beta_val[t] = static_cast<double>(t_SSP_val[i_t_SSP] - t) / dt;
                    // 実質的に補間しているのはこれだけ
                    n_Lym_val2[i_Z][t] = n_HI_tot_pre + (n_HI_tot - n_HI_tot_pre) *
                                                            static_cast<double>(t - t_SSP_val[i_t_SSP - 1]) / dt;

                    // これらの値は補間しないでそのまま右側の値を使う
                    n_SNII_val2[i_Z][t] = n_SNII;
                    n_SNIa_val2[i_Z][t] = n_SNIa;
                    ejecta_val2[i_Z][t] = ejecta;
                    ejecta_Z_val2[i_Z][t] = ejecta_Z;
                    m_BHNS_val2[i_Z][t] = m_BHNS;
                    m_WD_val2[i_Z][t] = m_WD;

                    // Mass in stars still alive.
                    m_alive_val2[i_Z][t + 1] = m_alive_val2[i_Z][t] - ejecta_val2[i_Z][t] - m_BHNS_val2[i_Z][t] - m_WD_val2[i_Z][t];
                }
                n_HI_tot_pre = n_HI_tot;
            }
        }
        // 後ろの境界条件
        t_inv_val[t_SSP_val[n_time_SSP - 1]] = n_time_SSP - 1 - 1;
        beta_val[t_SSP_val[n_time_SSP - 1]] = 0.0;
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
