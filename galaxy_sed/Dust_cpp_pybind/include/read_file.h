#ifndef ASANO_MODEL_READ_FILE_H
#define ASANO_MODEL_READ_FILE_H

#include <string>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>

#include <constants.h>
#include <file_name.h>
#include <file_util.h>
#include <dust_util.h>

#include "function.h"

namespace asano_model
{

    using val = std::valarray<double>;
    using val2 = std::valarray<val>;
    using val3 = std::valarray<val2>;
    using val4 = std::valarray<val3>;

    using SEDFile = my_util::SED_util::SEDFileUtil;

    inline val2 ReadSNDustSizeDistributionFile(const std::string &m_star_string)
    {
        auto f_SN_val2 = val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES - 1);
        if (SN_DATA_TYPE == 0)
        { // data of Nozawa et al.(2007)
            auto ifs = SEDFile::IfsOpen(SNDustSizeDistributionFileName(m_star_string));

            for (auto i_a = std::size_t(0); i_a < N_MAX_DUST_RADIUS; ++i_a)
            {
                SEDFile::ReadDouble(ifs); // Skip dust radius [cm]
                for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
                {
                    ifs >> f_SN_val2[i_s][i_a];
                    //                if (i_a < N_RADIUS_MIN) f_SN_val2[i_s][i_a] = 0;
                }
                SEDFile::ReadDouble(ifs); // Skip total dust mass
            }
        }
        else
        { // Bianchi & Schneider (2007)
            auto ifp = SEDFile::IfpOpen(SNdustSizeDistributionBsFileName());
            for (auto i_a = std::size_t(0); i_a < N_MAX_DUST_RADIUS; ++i_a)
            {
                auto buff = 0.0;
                if (fscanf(ifp, "%lf", &buff) != 1)
                {
                    fprintf(stderr, "Error reading value for buff at index %zu\n", i_a);
                    break; // 読み取りエラーが発生した場合、ループを終了
                }

                if (fscanf(ifp, "%lf", &f_SN_val2[0][i_a]) != 1)
                {
                    fprintf(stderr, "Error reading f_SN_val2[0][%zu] at index %zu\n", i_a, i_a);
                    break; // 読み取りエラーが発生した場合、ループを終了
                }
                for (auto i_s = std::size_t(1); i_s < N_DUST_SPECIES; ++i_s)
                {
                    f_SN_val2[i_s][i_a] = f_SN_val2[0][i_a];
                }
                if (feof(ifp))
                    break;
            }
        }
        return f_SN_val2;
    }

    /**
     * Read dust size distribution produced by supernova
     * @param m_star_string progenitor stellar mass in std::string type [Msun]
     * @param a_val dust radius units [cm]
     * @return val2
     */
    inline val2
    ReadAndNormalizeSNDustSizeDistributionFile(const std::string &m_star_string, const val &a_val)
    {
        auto f_SN_val2 = ReadSNDustSizeDistributionFile(m_star_string);

        const auto da_val = DeltaRadius(a_val);
        const auto volume_val = Map(a_val, Volume);

        for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
        {
            const auto total = (f_SN_val2[i_s] * volume_val * da_val).sum();
            for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
            {
                if (total < DBL_MIN)
                    f_SN_val2[i_s][i_a] = 0;
                else
                    f_SN_val2[i_s][i_a] *= volume_val[i_a] * da_val[i_a] / total;
            }
        }

        return f_SN_val2;
    }

    /**
     * Read dust mass produced by SN from file
     * The file is written the result Figure.13 in Nozawa et al. (2003) and Table.3 in Nozawa et al. (2007)
     * @param progenitor_stellar_mass units in solar mass
     * @return val [dust_species] C Si Fe FeS Al2O3 MgSiO3 Mg2SiO4 SiO2 MgO
     */
    inline val ReadStellarMassFile(const std::string &stellar_mass)
    {
        auto ifs = SEDFile::IfsOpen(SNDustMassFileName(stellar_mass));

        auto i_n_H = std::size_t(0);
        if (n_H <= 0.03)
        {
            i_n_H = 0;
        }
        else if (n_H <= 0.1)
        {
            i_n_H = 1;
        }
        else if (n_H <= 0.3)
        {
            i_n_H = 2;
        }
        else if (n_H <= 1.0)
        {
            i_n_H = 3;
        }
        else if (n_H <= 3.0)
        {
            i_n_H = 4;
        }
        else if (n_H <= 10.0)
        {
            i_n_H = 5;
        }
        else
        {
            i_n_H = 6;
        }
        SEDFile::SkipLine(ifs); // SKip header
        for (auto i = std::size_t(0); i < i_n_H; ++i)
            SEDFile::SkipLine(ifs); // SKip to i_n_H

        auto m_SN_per_species = val(N_DUST_SPECIES - 1);
        for (auto i_species = std::size_t(0); i_species < N_DUST_SPECIES - 1; ++i_species)
        {
            ifs >> m_SN_per_species[i_species];
        }

        if (isSN_ADJUST == 1)
        {
            std::for_each(std::begin(m_SN_per_species), std::end(m_SN_per_species),
                          [](auto x)
                          {
                              return x * 0.1;
                          });
        }
        return m_SN_per_species;
    }

    /**
     * Read dust size distribution produced by AGB stars from file
     * Log-normal distribution
     * @param ifn In file name. The file is from Yasuda & Kozasa (2012).
     * @param a_val Dust radius [cm]
     * @param n_AGB_per_radius Dust size distribution per radius [cm^-4]
     */
    inline void
    ReadAGBDustSizeDistributionFile(const std::string &ifn, val &a_val, val &n_AGB_per_radius)
    {
        auto ifp = SEDFile::IfpOpen(ifn);
        for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
        {
            if (fscanf(ifp, "%lf", &a_val[i]) != 1)
            {
                fprintf(stderr, "Error reading a_val[%lu]\n", i);
                break; // 読み取りエラーが発生した場合、ループを終了
            }

            if (fscanf(ifp, "%lf", &n_AGB_per_radius[i]) != 1)
            {
                fprintf(stderr, "Error reading n_AGB_per_radius[%lu]\n", i);
                break; // 読み取りエラーが発生した場合、ループを終了
            }
            if (feof(ifp))
                break;
        }
        fclose(ifp);
    }

    /**
     * Read file and normalize dust size distribution produced by AGB stars.
     * @return val Normalized dust size distribution produced by AGB stars. f_AGB_per_radius_val.sum() == 1
     */
    inline val ReadAndNormalizeAGBDustSizeDistributionFile()
    {
        auto n_AGB_per_a = val(N_MAX_DUST_RADIUS);
        auto a_val = val(N_MAX_DUST_RADIUS);
        ReadAGBDustSizeDistributionFile(AGBDustSizeDistributionFileName(), a_val, n_AGB_per_a);
        const auto da_val = DeltaRadius(a_val);
        const auto volume_val = Map(a_val, Volume);
        const auto M_total = (volume_val * n_AGB_per_a * da_val).sum();
        const auto f_AGB_per_a_val = n_AGB_per_a * volume_val * da_val / M_total;
        return f_AGB_per_a_val;
    }

    /**
     * Read dust mass distribution each dust species produced by AGB stars from files.
     * @return val3 [metallicity][progenitor mass][dust species].
     * Dust species: stellar_mass, forsterite, fayalite, enstatite, ferrosilite, quartz(Mstar), iron(Mstar), quartz(Sstar), iron(Sstar), SiC, carbon, iron(Cstar), and total.
     */
    inline val3 ReadAGBMassFile()
    {
        if (AGB_DATA_TYPE == 0)
        { // data of Zhukovska
            const auto AGB_Z_array = std::array<std::string, N_METALICITY_AGB>(
                {"0.001", "0.002", "0.004", "0.008", "0.015", "0.020"});

            auto m_AGB_per_species = val3(val2(val(13), N_M_AGB_STAR), N_METALICITY_AGB);

            for (auto i_Z = std::size_t(0); i_Z < N_METALICITY_AGB; ++i_Z)
            {
                auto ifs = SEDFile::IfsOpen(AGBDustMassFileName(AGB_Z_array[i_Z]));
                SEDFile::SkipLine(ifs); // Skip header
                for (auto i_m = std::size_t(0); i_m < N_M_AGB_STAR; ++i_m)
                {
                    SEDFile::ReadDouble(ifs); // Skip progenitor stellar mass
                    for (auto i_s = std::size_t(0); i_s < 11; ++i_s)
                    {
                        ifs >> m_AGB_per_species[i_Z][i_m][i_s];
                    }
                    SEDFile::ReadDouble(ifs); // Skip total dust mass
                }
            }
            return m_AGB_per_species;
        }
        else
        { // data of Ventura
            const auto AGB_metalicity_array = std::array<std::string, 2>({"0.001", "0.008"});
            auto m_AGB_per_species = val3(val2(val(8), 14), 2);
            for (auto i_m_star = std::size_t(0); i_m_star < 2; ++i_m_star)
            {
                auto filename = AGBDustMassVenturaFileName(AGB_metalicity_array[i_m_star]);
                auto ifp = SEDFile::IfpOpen(filename);

                for (auto i = std::size_t(0); i < 14; ++i)
                {
                    for (auto k = std::size_t(0); k < 8; ++k)
                    {
                        if (fscanf(ifp, "%lf", &m_AGB_per_species[i_m_star][i][k]) != 1)
                        {
                            fprintf(stderr, "Error reading m_AGB_per_species[%lu][%lu][%lu]\n", i_m_star, i, k);
                            break;
                        }
                    }
                    if (feof(ifp))
                        break;
                }
            }
            return m_AGB_per_species;
        }
    }

    inline val2 NormalizeSNDustMass(const val &m_SN_per_species, const val2 &f_SN_val2) noexcept
    {
        auto M_dust_SN_val2 = val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES - 1);

        for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
        {
            for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
            {
                M_dust_SN_val2[i_s][i_a] = f_SN_val2[i_s][i_a] * m_SN_per_species[i_s];
            }
        }
        return M_dust_SN_val2;
    }

    /**
     * Weighting AGB dust mass by dust size distribution
     * @param f_AGB_val Fraction of dust mass each dust radius bin
     * @param M_AGB_total_val3 Total dust mass each dust species
     * @return val4 [metallicity][progenitor stellar mass][dust species][dust radius]
     */
    inline val4 WeightingAGBDustMass(const val &f_AGB_val, const val3 &M_AGB_total_val3) noexcept
    {
        auto M_dust_AGB_val4 = val4(
            val3(val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES - 1), N_M_AGB_STAR), N_METALICITY_AGB);

        if (AGB_DATA_TYPE == 0)
        {
            for (auto i_Z = std::size_t(0); i_Z < N_METALICITY_AGB; ++i_Z)
            {
                for (auto i_m = std::size_t(0); i_m < N_M_AGB_STAR; ++i_m)
                {
                    for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
                    {
                        M_dust_AGB_val4[i_Z][i_m][0][i_a] =
                            M_AGB_total_val3[i_Z][i_m][9] * f_AGB_val[i_a]; // carbon

                        M_dust_AGB_val4[i_Z][i_m][2][i_a] =
                            (M_AGB_total_val3[i_Z][i_m][5] + M_AGB_total_val3[i_Z][i_m][7] +
                             M_AGB_total_val3[i_Z][i_m][10]) *
                            f_AGB_val[i_a]; // iron

                        M_dust_AGB_val4[i_Z][i_m][5][i_a] =
                            (M_AGB_total_val3[i_Z][i_m][2] + M_AGB_total_val3[i_Z][i_m][3]) *
                            f_AGB_val[i_a]; // enstatite and ferrosilite to MgSiO3

                        M_dust_AGB_val4[i_Z][i_m][6][i_a] =
                            (M_AGB_total_val3[i_Z][i_m][0] + M_AGB_total_val3[i_Z][i_m][1]) *
                            f_AGB_val[i_a]; // forsterite and fayalite to Mg2SiO4

                        M_dust_AGB_val4[i_Z][i_m][7][i_a] =
                            (M_AGB_total_val3[i_Z][i_m][4] + M_AGB_total_val3[i_Z][i_m][6] +
                             M_AGB_total_val3[i_Z][i_m][8]) *
                            f_AGB_val[i_a]; // quartz and SiC to SiO2
                    }
                }
            }
        }
        else
        {
            for (auto i_m_star = std::size_t(0); i_m_star < 2; ++i_m_star)
            {
                for (auto i = std::size_t(0); i < 14; ++i)
                {
                    //                m_normalized_AGB[i_m_star][i][0] =
                    //                m_AGB_per_species[i_m_star][i][6] / m / RHO_GRAIN[0];  //carbon dust
                    //                m_normalized_AGB[i_m_star][i][2] =
                    //                m_AGB_per_species[i_m_star][i][5] / m / RHO_GRAIN[2];  //iron dust
                    //                m_normalized_AGB[i_m_star][i][5] =
                    //                m_AGB_per_species[i_m_star][i][2] / m / RHO_GRAIN[5];  //olivine
                    //                m_normalized_AGB[i_m_star][i][6] =
                    //                m_AGB_per_species[i_m_star][i][3] / m / RHO_GRAIN[6];  //pyroxene
                    //                m_normalized_AGB[i_m_star][i][7] =
                    //                (m_AGB_per_species[i_m_star][i][4] + m_AGB_per_species[i_m_star][i][7]) / m
                    //                / RHO_GRAIN[7];  //quartz and silicon carbide
                }
            }
        }
        return M_dust_AGB_val4;
    }

    inline val2 ReadTotalMetalMassFile(const std::string &metalicity)
    {
        std::string filename_str = TotalMetalMassFileName(metalicity);
        const char *filename = filename_str.c_str();
        auto ifp = SEDFile::IfpOpen(filename);

        auto M_Z_per_star = val2(val(3), N_M_STAR);

        for (auto i = std::size_t(0); i < N_M_STAR; ++i)
        {
            auto m_star = 0.0;
            if (fscanf(ifp, "%lf", &m_star) != 1)
            {
                fprintf(stderr, "Error reading m_star from file: %s\n", filename); // ファイル名を表示

                break; // 読み取りエラーが発生した場合、ループを終了
            }
            for (auto j = std::size_t(0); j < 3; ++j)
            {
                if (fscanf(ifp, "%lf", &M_Z_per_star[i][j]) != 1)
                {
                    fprintf(stderr, "Error reading M_Z_per_star[%lu][%lu]\n", i, j);
                    break; // 読み取りエラーが発生した場合、ループを終了
                }
                if (feof(ifp))
                    break; // ファイルの末尾まで読み込み
            }
        }
        fclose(ifp);
        return M_Z_per_star;
    }

    inline val ReadTotalRemnantMassFile(const std::string &metalicity)
    {
        auto ifp = SEDFile::IfpOpen(TotalRemnantMassFileName(metalicity));
        auto M_rem_val = val(N_M_STAR);

        for (auto i_m_star = std::size_t(0); i_m_star < N_M_STAR; ++i_m_star)
        {
            auto m_star = 0.0;
            if (fscanf(ifp, "%lf", &m_star) != 1)
            {
                fprintf(stderr, "Error reading m_star\n");
                break; // 読み取りエラーが発生した場合、ループを終了
            }

            if (fscanf(ifp, "%lf", &M_rem_val[i_m_star]) != 1)
            {
                fprintf(stderr, "Error reading M_rem_val[%lu]\n", i_m_star);
                break; // 読み取りエラーが発生した場合、ループを終了
            }
            if (feof(ifp))
                break; // ファイルの末尾まで読み込み
        }
        fclose(ifp);
        return M_rem_val;
    }

    inline val3 ReadDestEtaFile(const std::string &metalicity)
    {
        auto ifp = SEDFile::IfpOpen(DestEtaFileName(metalicity));
        auto dest_eta_val3 = val3(val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS);

        for (auto i_a = std::size_t(0); i_a < N_MAX_DUST_RADIUS; ++i_a)
        {
            auto a = 0.0;
            if (fscanf(ifp, "%lf", &a) != 1)
            {
                fprintf(stderr, "Error reading value for a\n");
                break; // 読み取りエラーが発生した場合、ループを終了
            }

            for (auto j_a = std::size_t(0); j_a < N_MAX_DUST_RADIUS; ++j_a)
            {
                for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES; ++i_s)
                {
                    if (fscanf(ifp, "%lf", &dest_eta_val3[i_a][j_a][i_s]) != 1)
                    {
                        fprintf(stderr, "Error reading dest_eta_val3[%zu][%zu][%zu]\n", i_a, j_a, i_s);
                        break; // 読み取りエラーが発生した場合、内側のループを終了
                    }
                }
            }
            if (feof(ifp))
                break;
        }
        fclose(ifp);
        return dest_eta_val3;
    }

    /**
     * Read dust mass produced by supernova
     * @param a_val dust radius
     * @return val3 [progenitor stellar mass][dust species][dust radius]
     */
    inline val3 ReadSNDustFile(const val &a_val)
    {
        auto M_dust_SN_val3 = val3(val2(val(N_MAX_DUST_RADIUS), N_DUST_SPECIES), N_M_SN_STAR);

        const auto m_star = std::array<std::string, 4>({"13", "20", "25", "30"});
        for (auto i_m_star = std::size_t(0); i_m_star < 4; ++i_m_star)
        {
            const auto f_SN_val2 = ReadAndNormalizeSNDustSizeDistributionFile(m_star[i_m_star],
                                                                              a_val);
            const auto M_SN_total_val2 = ReadStellarMassFile(m_star[i_m_star]);
            M_dust_SN_val3[i_m_star] = NormalizeSNDustMass(M_SN_total_val2, f_SN_val2);
        }
        return M_dust_SN_val3;
    }

    /**
     * Read dust mass distribution produced by AGB stars
     * @return val4 [metallicity][progenitor mass][dust species][dust radius]
     */
    inline val4 ReadAGBDustFile()
    {
        const auto f_AGB_per_a_val = ReadAndNormalizeAGBDustSizeDistributionFile();
        const auto M_AGB_total_val3 = ReadAGBMassFile();
        return WeightingAGBDustMass(f_AGB_per_a_val, M_AGB_total_val3);
    }

    inline val3 ReadMetalYieldFile()
    {
        auto M_Z_per_star_val3 = val3(val2(val(3), N_M_STAR), N_METALLICITY);
        for (auto i_Z = std::size_t(5); i_Z < N_METALLICITY; ++i_Z)
        {
            auto Z = std::stringstream();
            Z << std::fixed << std::setprecision(3) << static_cast<double>(i_Z) / 1000.0;
            M_Z_per_star_val3[i_Z] = ReadTotalMetalMassFile(Z.str());
        }
        return M_Z_per_star_val3;
    }

    inline val2 ReadRemnantFile()
    {
        auto M_rem_val2 = val2(val(N_M_STAR), N_METALLICITY);

        for (auto i_Z = std::size_t(5); i_Z < N_METALLICITY; ++i_Z)
        {
            auto Z = std::stringstream();
            Z << std::fixed << std::setprecision(3) << static_cast<double>(i_Z) / 1000.0;
            M_rem_val2[i_Z] = ReadTotalRemnantMassFile(Z.str());
        }
        return M_rem_val2;
    }

    inline val4 ReadSNDestructionFile()
    {
        auto dest_eta_val4 = val4(val3(val2(val(N_DUST_SPECIES), N_MAX_DUST_RADIUS), N_MAX_DUST_RADIUS),
                                  4);

        const auto SND_Z_array = std::array<std::string, 4>({"-4", "-2", "-1", "0"});

        for (auto i_Z = std::size_t(0); i_Z < 4; ++i_Z)
        {
            dest_eta_val4[i_Z] = ReadDestEtaFile(SND_Z_array[i_Z]);
        }
        return dest_eta_val4;
    }

    inline val2 ReadGrainVelocityFile(const std::string &ISM_type)
    {
        auto grain_velocity = val2(val(2), N_MAX_DUST_RADIUS);
        auto ifs_carbon = SEDFile::IfpOpen(GrainVelocityFileName(ISM_type, "car"));
        for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
        {
            auto radius = 0.0;
            // fscanf(ifs_carbon, "%lf   %lf", &radius, &grain_velocity[i][0]);
            if (fscanf(ifs_carbon, "%lf   %lf", &radius, &grain_velocity[i][0]) != 2)
            {
                // 予期しないエラー処理（エラーログを出力する、異常終了するなど）
                fprintf(stderr, "Error reading carbon grain velocity data at index %lu\n", i);
            }
            if (feof(ifs_carbon))
                break;
        }
        fclose(ifs_carbon);
        auto ifs_silicate = SEDFile::IfpOpen(GrainVelocityFileName(ISM_type, "sil"));
        for (auto i = std::size_t(0); i < N_MAX_DUST_RADIUS; ++i)
        {
            auto radius = 0.0;
            if (fscanf(ifs_silicate, "%lf   %lf", &radius, &grain_velocity[i][1]) != 2)
            {
                // 予期しないエラー処理（エラーログを出力する、異常終了するなど）
                fprintf(stderr, "Error reading carbon grain velocity data at index %lu\n", i);
            }
            if (feof(ifs_silicate))
                break;
        }
        fclose(ifs_silicate);
        return grain_velocity;
    }

    inline FILE *MakeTotalMassFile(const std::string &ofn)
    {
        auto ofp_m_total = SEDFile::OfpOpen(ofn);
        fprintf(ofp_m_total,
                "age[Myr] M_galaxy[Msun] M_gas[Msun] M_star[Msun] m_gas_carbon[Msun] m_gas_silicate[Msun] M_Z_val[Msun] m_dust_C[Msun] m_dust_Si[Msun], m_dust_Fe[Msun], m_dust[Msun] metallicity[/Zsun] m_dust/M_gas m_dust/M_star\n");
        return ofp_m_total;
    }

    /**
     * Make file and write header. n_Sil represents total silicate number (= n_dust - n_C)
     * @param ofn output file name
     * @return FILE*
     */
    inline FILE *MakeNumberDistributionFile(const std::string &ofn)
    {
        auto ofp_n = SEDFile::OfpOpen(ofn);
        fprintf(ofp_n,
                "radius[cm] n_C n_Sil n_Si n_Fe n_FeS n_Al2O3 n_MgSiO3 n_Mg2SiO4 n_SiO2 n_MgO n_dust\n");
        return ofp_n;
    }

    /**
     * Make file and write header. m_Sil represents total silicate mass (= m_dust - m_C)
     * @param ofn output file name
     * @return FILE*
     */
    inline FILE *MakeMassDistributionFile(const std::string &ofn)
    {
        auto ofp_m = SEDFile::OfpOpen(ofn);
        fprintf(ofp_m,
                "radius[cm] m_C[Msun] m_Sil[Msun] m_Si[Msun] m_Fe[Msun] m_FeS[Msun] m_Al2O3[Msun] m_MgSiO3[Msun] m_Mg2SiO4[Msun] m_SiO2[Msun] m_MgO[Msun] m_dust[Msun]\n");
        return ofp_m;
    }

    /**
     * Write dust size distribution
     * @param i_age Index of age
     * @param m_metal Metal mass in gas
     * @param m_gas Gas mass
     * @param m_carbon Carbon mass in gas
     * @param m_silicate Silicate mass in gas
     * @param m_star Stellar mass
     * @param a_val dust radius [cm]
     * @param V_val volume [cm^3]
     * @param m_val2 Dust mass each dust species and radius
     * @param ofp_m_total File pointer of total dust mass
     * @param ofp_n File pointer of number distribution
     * @param ofp_m File pointer of mass distribution
     */
    inline void
    WriteFiles(std::size_t i_age, double m_metal, double m_gas, double m_carbon, double m_silicate,
               double m_star, const val &a_val, const val &V_val, const val2 &m_val2,
               FILE *ofp_m_total, FILE *ofp_n, FILE *ofp_m)
    {
        auto m_total_per_species_val = val(N_DUST_SPECIES - 1);
        auto m_total_per_radius_val = val(N_MAX_DUST_RADIUS);
        auto m_total = 0.0;

        auto n_total_per_radius = val(N_MAX_DUST_RADIUS);
        auto n_val2 = val2(val(N_DUST_SPECIES - 1), N_MAX_DUST_RADIUS);

        for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
        {
            m_total += m_val2[i_a].sum();
            for (auto i_s = std::size_t(0); i_s < N_DUST_SPECIES - 1; ++i_s)
            {
                n_val2[i_a][i_s] = m_val2[i_a][i_s] * Msun / (V_val[i_a] * RHO_GRAIN[i_s]);
                n_total_per_radius[i_a] += n_val2[i_a][i_s];

                m_total_per_species_val[i_s] += m_val2[i_a][i_s];
                m_total_per_radius_val[i_a] += m_val2[i_a][i_s];
            }
        }

        const auto age = static_cast<std::size_t>(static_cast<double>(i_age) * TIME_BIN_DUST_MODEL /
                                                  1e6);
        fprintf(ofp_m_total,
                "%lu %5e %5e %5e %5e %5e %5e %5e %5e %5e %5e %5e %5e %5e\n",
                age, m_gas + m_star, m_gas, m_star, m_carbon, m_silicate, m_metal,
                m_total_per_species_val[0], m_total_per_species_val[1], m_total_per_species_val[2],
                m_total, m_metal / m_gas / Zsun, m_total / m_gas, m_total / m_star);

        for (auto i_a = N_RADIUS_MIN; i_a < N_MAX_DUST_RADIUS; ++i_a)
        {
            fprintf(ofp_n,
                    "%5e %5e %5e %5e %5e %5e %5e %5e %5e %5e %5e %5e\n",
                    a_val[i_a], n_val2[i_a][0],
                    std::accumulate(std::begin(n_val2[i_a]) + 1, std::end(n_val2[i_a]),
                                    0.0), // Silicate
                    n_val2[i_a][1], n_val2[i_a][2], n_val2[i_a][3], n_val2[i_a][4],
                    n_val2[i_a][5], n_val2[i_a][6], n_val2[i_a][7], n_val2[i_a][8],
                    n_total_per_radius[i_a]);

            fprintf(ofp_m,
                    "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
                    a_val[i_a], m_val2[i_a][0],
                    std::accumulate(std::begin(m_val2[i_a]) + 1, std::end(m_val2[i_a]),
                                    0.0), // Silicate
                    m_val2[i_a][1], m_val2[i_a][2], m_val2[i_a][3], m_val2[i_a][4],
                    m_val2[i_a][5], m_val2[i_a][6], m_val2[i_a][7], m_val2[i_a][8],
                    m_total_per_radius_val[i_a]);
        }
    }
}
#endif // ASANO_MODEL_READ_FILE_H
