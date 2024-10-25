/*
 * file_util.h
 *
 *  Created on: 2017/08/04
 *      Author: kazuki
 */

#ifndef FILE_UTIL_H_
#define FILE_UTIL_H_

#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>

#include <container_util.h>

namespace my_util {
class FileUtil {
  private:
    template<class T>
    static void Open(const std::string& file_name, const std::ios::openmode& mode, T& fs) {
        fs.open(file_name, mode);
        if (fs.fail()) {  // オープン失敗
            std::cout << "Error: \"" << file_name << "\" can not be opened !!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    static std::FILE* Open(const std::string& file_name, const std::string& mode) {
        auto fp = std::fopen(file_name.c_str(), mode.c_str());
        if (!fp) {  // オープン失敗
            std::cout << "Error: \"" << file_name << "\" can not be opened !!" << std::endl;
            exit(EXIT_FAILURE);
        }
        return fp;
    }

    static bool CheckFileExistence(const std::string& file_name) noexcept {
        const auto ifs = std::ifstream(file_name, std::ios::in);
        return ifs.is_open();
    }

    // すでにあるファイル名を指定していて上書きしないよう場合に使う
    // file name の後ろに_2 を追加する。_2がある場合もさらに_2が追加される
    template<class T>
    static void AddNumberFileNameEnd(T& path) noexcept {
        // 例：../../Calculation/file_name_1.dat の場合をコメントに記す
        auto path_vector      = Split(path, "/"); // [.. .. Calculation file_name_1.dat]
        auto file_name_vector = Split(*(std::end(path_vector) - 1), "."); // [file_name_1 dat]
        auto base_end_vector  = Split(*(std::end(file_name_vector) - 2), "_"); // [file name 1]
        auto base_end         = *(std::end(base_end_vector) - 1); // [1]

        if (CheckInt(base_end)) {
            // 末尾が数字のときはその数字に1を足したものをファイル名にする
            // base_end_vector = [file name 2]
            *(std::end(base_end_vector) - 1) = std::to_string(std::stoi(base_end) + 1);
        } else {
            // 末尾が数字出ない場合は "_1" を追加
            *(std::end(base_end_vector) - 1) += "_1";
        }

        // file_name_end_vector を再び "_" で繋いで一つの文字列にする
        const auto base = Combine(base_end_vector, "_"); // "file_name_2"s
        *(std::end(file_name_vector) - 2) = base; // [file_name_2 dat]

        const auto file_name = Combine(file_name_vector, "."); // "file_name_2.dat"s
        *(std::end(path_vector) - 1) = file_name; // [.. .. Calculation file_name_2.dat]

        path = Combine(path_vector, "/"); // "../../Calculation/file_name_2.dat"s
    }

    template<class T>
    static bool CheckInt(const T& input) {
        try {
            std::stoi(input);
        } catch (...) {
            return false;
        }
        return true;
    }

  public:
    FileUtil() = default;

    virtual ~FileUtil() = default;

    static std::ifstream IfsOpen(const std::string& filename) {
        auto ifs = std::ifstream();
        Open(filename, std::ios::in, ifs);
        return ifs;
    }

    /**
     * Open file with read only mode
     * @param filename
     * @return FILE pointer
     */
    static std::FILE* IfpOpen(const std::string& filename) {
        return Open(filename, "r");
    }


    static std::ofstream OfsOpen(const std::string& file_name, const bool print_message = false,
                                 const bool is_over_write = true) {
        auto ofs = std::ofstream();
        if (is_over_write) {
            Open(file_name, std::ios::out, ofs);
            if (print_message)
                std::cout << "Open file: \"" << file_name << "\"" << std::endl;
            return ofs;
        }
        auto new_file_name = file_name;
        while (CheckFileExistence(new_file_name)) {
            // 同名ファイルが存在する場合
            AddNumberFileNameEnd(new_file_name);
        }
        Open(new_file_name, std::ios::out, ofs);
        if (print_message)
            std::cout << "Open file: \"" << new_file_name << "\"" << std::endl;
        return ofs;
    }

    static std::FILE* OfpOpen(const std::string& file_name, const bool print_message = false,
                              const bool is_over_write = true) {
        if (is_over_write) {
            auto ofp = Open(file_name, "w");
            if (print_message)
                std::cout << "Open file: \"" << file_name << "\"" << std::endl;
            return ofp;
        }

        auto new_file_name = file_name;
        while (CheckFileExistence(new_file_name)) {
            // 同名ファイルが存在する場合
            AddNumberFileNameEnd(new_file_name);
        }
        auto ofp = Open(new_file_name, "w");
        if (print_message)
            std::cout << "Open file: \"" << new_file_name << "\"" << std::endl;
        return ofp;
    }

    // 入出力両方を行う場合の open
    static std::fstream AfsOpen(const std::string& file_name, const bool is_over_write = false) {
        auto afs = std::fstream();
        if (is_over_write) {
            Open(file_name, std::ios::app, afs);
            // Debug
            std::cout << "Open file: \"" << file_name << "\"" << std::endl;
            return afs;
        }
        auto new_file_name = file_name;
        while (CheckFileExistence(new_file_name)) {
            // 同名ファイルが存在する場合
            AddNumberFileNameEnd(new_file_name);
        }
        Open(new_file_name, std::ios::app, afs);
        std::cout << "Open file: \"" << new_file_name << "\"" << std::endl;
        return afs;
    }

    // 入出力両方を行う場合の open
    static std::FILE* AfpOpen(const std::string& file_name, const bool is_over_write = false) {
        if (is_over_write) {
            auto afp = Open(file_name, "a");
            // Debug
            std::cout << "Open file: \"" << file_name << "\"" << std::endl;
            return afp;
        }
        auto new_file_name = file_name;
        while (CheckFileExistence(new_file_name)) {
            // 同名ファイルが存在する場合
            AddNumberFileNameEnd(new_file_name);
        }
        auto afp = Open(new_file_name, "a");
        std::cout << "Open file: \"" << new_file_name << "\"" << std::endl;
        return afp;
    }

    // Filename を与えた場合と、ファイルオブジェクトを与えた場合でオーバーロードしてある
    template<class MapType>
    static void WriteMapToFile(const MapType& map, const std::string& filename) {
        auto ofs = OfsOpen(filename);
        for (const auto&& p : map) {
            ofs << std::scientific << std::setprecision(3) << p.first << " " << p.second
                << std::endl;
        }
        ofs.close();
    }

    template<class MapType, class FileType>
    static void WriteMapToFile(const MapType& map, FileType& ofs) {
        for (const auto& p : map) {
            ofs << std::scientific << std::setprecision(10) << p.first << " " << p.second
                << std::endl;
        }
    }

    template<class MapType>
    static void WriteMapToFile(const MapType& map, FILE* ofp) {
        for (const auto& p : map) {
            std::fprintf(ofp, "%1.10e %1.10e\n", p.first, p.second);
        }
    }

    // ファイルからデータを読み込んで、first, second のじゅんで詰めていく
    template<typename FirstType, typename SecondType, class FileType>
    static bool ReadAndInsertMap(FileType& ifs, std::map<FirstType, SecondType>& map) {
        FirstType  first;
        SecondType second;
        if (ifs >> first) {
            ifs >> second;
            map.insert(std::make_pair(first, second));
            return true;
        } else {
            return false;
        }
    }

// StreamType は fstream 系か stringstream を入れられる
    template<class StreamType, class ContainerType>
    static bool
    ReadAndPushBack(StreamType& stream, ContainerType& contianer, double multiplicand = 1) {
        double value = 0.0;
        if (stream >> value) {
            contianer.push_back(value * multiplicand);
            return true;
        }
        return false;
    }

    //! getline without create unused variable
    inline static auto SkipLine(FILE* ifp) {
        char* buff = nullptr;
        std::fgets(buff, 256, ifp);
        return buff;
    }

    //! Read std::string without create unused variable
    template<typename CharT, typename Traits>
    inline static std::string ReadString(std::basic_istream<CharT, Traits>& is) {
        auto str = std::string();
        is >> str;
        return str;
    }

    //! Read double without create unused variable
    template<typename CharT, typename Traits>
    inline static double ReadDouble(std::basic_istream<CharT, Traits>& is) {
        auto v = 0.0;
        is >> v;
        return v;
    }

    //! getline without create unused variable
    template<typename CharT, typename Traits>
    inline static std::string SkipLine(std::basic_istream<CharT, Traits>& is) {
        auto str = std::string();
        std::getline(is, str);
        return str;
    }

    template<typename CharT, typename Traits>
    inline static std::string GetLineByString(std::basic_istream<CharT, Traits>& is) {
        //inline static std::stringstream GetLineToStringStream(std::ifstream& is) {
        auto line_string = std::string();
        std::getline(is, line_string);
        return line_string;
    }

    template<typename CharT, typename Traits>
    inline static std::stringstream GetLineToStringStream(std::basic_istream<CharT, Traits>& is) {
        //inline static std::stringstream GetLineToStringStream(std::ifstream& is) {
        auto line_string = std::string();
        std::getline(is, line_string);
        return std::stringstream(line_string);
    }

    template<typename CharT, typename Traits>
    inline static std::stringstream
    GetLineToStringStream(std::basic_istream<CharT, Traits>& is, const std::string& delim) {
        auto line_ss = GetLineToStringStream(is);
        auto buff    = std::string();
        std::getline(line_ss, buff, *delim.c_str()); // Split by delim
        return line_ss;
    }

    template<typename CharT, typename Traits>
    inline static std::string SkipOneWord(std::basic_istream<CharT, Traits>& is) {
        auto s = std::string();
        is >> s;
        return s;
    }

};
} /* namespace my_util */
#endif /* FILE_UTIL_H_ */
