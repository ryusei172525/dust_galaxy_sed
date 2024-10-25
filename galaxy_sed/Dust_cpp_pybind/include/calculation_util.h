//====================================================================
// Name        : calculation_util.h
// Author      : kazuki
// Version     : 1.0
// Copyright   : Your copyright notice
// Date        : 2018/06/22
// Input file  :
// Output file :
// Description :
// Usage       :
//====================================================================
#ifndef CALCULATION_UTIL_H_
#define CALCULATION_UTIL_H_

#include <algorithm>
#include <ctime>
#include <iostream>
#include <string>
#include <chrono>

namespace my_util {
class MyTime {
  public:
    MyTime() = default;

    ~MyTime() = default;

    static time_t PrintDateTime(const std::string& message) noexcept {
        const auto now   = std::time(nullptr);
        const auto p_now = std::localtime(&now);
        // Transform time format for people. (e.g. 2018/6/12 17:42:39)
        std::cout << p_now->tm_year + 1900 << "/" << p_now->tm_mon + 1 << "/" << p_now->tm_mday
                  << " "
                  << p_now->tm_hour << ":" << p_now->tm_min << ":" << p_now->tm_sec << " "
                  << message
                  << std::endl;
        return now;
    }

    static auto GetTime() noexcept {
        return std::chrono::system_clock::now();
    }

    static void
    PrintElapsedTime(const std::chrono::system_clock::time_point& start,
                     const std::string& unit = "micro", const std::string& comment = "") noexcept {
        const auto end = std::chrono::system_clock::now();
        std::cout << comment;
        if (unit == "sec") {
            std::cout << std::chrono::duration_cast<std::chrono::seconds>(
            end - start).count() << " sec\n";
        } else if (unit == "milli") {
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            end - start).count() << " msec\n";
        } else {
            std::cout << std::chrono::duration_cast<std::chrono::microseconds>(
            end - start).count() << " usec\n";
        }
    }


    static void PrintElapsedTime(const time_t& start_time,
                                 const time_t& end_time) noexcept {
        if (auto elapsed_time_sec = static_cast<double>(end_time - start_time);
        elapsed_time_sec < 60) {
            std::cout << "所要時間は" << elapsed_time_sec << "秒です。" << std::endl;
        } else if (const auto elapsed_time_min = elapsed_time_sec / 60.0; elapsed_time_min < 60) {
            std::cout << "所要時間は" << elapsed_time_min << "分です。" << std::endl;
        } else {
            std::cout << "所要時間は" << elapsed_time_min / 60.0 << "時間です。" << std::endl;
        }
    }
};

// Construct 時に container を2つ与えてそれをメンバ変数としてもつ。
// 線形補間を行う。log 補完をしたい場合は container の要素をすべて log にしておけば良い。
template<class XContainer, class YContainer, class XValueType>
class Interpolate {
  private:
    XContainer x_container_;
    YContainer y_container_;
    size_t     max_iteration_;
    bool       flag_valid_;

  public:
    // Check arguments are valid or not.
    Interpolate(XContainer x_container, YContainer y_container)
    : x_container_(x_container),
      y_container_(y_container),
      max_iteration_(x_container_.size() - 1),
      flag_valid_(true) {

        // 要素数が2個以上かどうかを調べる
        if (x_container_.size() < 2 || y_container_.size() < 2) {
            std::cout << "Error! Interpolate: The size of X or Y must be greater than 2!"
                      << std::endl;
            flag_valid_ = false;
            return;
        }
        // x_container と y_container の要素数が一致しているかどうかを調べる
        if (!(x_container_.size() == y_container_.size())) {
            std::cout << "Error! Interpolate: The size of X and Y must be same!" << std::endl;
            flag_valid_ = false;
            return;
        }

        // x が昇順にソートされているかどうかを調べる
        if (!std::is_sorted(std::begin(x_container_), std::end(x_container_))) {
            std::cout << "Error! Interpolate: X must be sorted in ascending!" << std::endl;
            flag_valid_ = false;
            return;
        }
        // y_container が昇順か降順にソートされているかどうかを調べる
        if (!std::is_sorted(std::begin(y_container_), std::end(y_container_)) &&
            !std::is_sorted(std::begin(y_container_), std::end(y_container_),
                            [](const auto& x, const auto& y) {
                                return x >= y;
                            })) {
            // y_container がソートされていなければエラー
            std::cout << "Error! Interpolate: Y must be sorted in ascending (decending)!"
                      << std::endl;
            flag_valid_ = false;
            return;
        }
    }

    // Calculate interpolate.
    double operator()(XValueType x) {
        // When given Invalid container in constructor.
        if (!flag_valid_) {
            std::cout << "Error! Interpolate: This object is invalid!" << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
        // When x is smaller than x_container range.
        if (x < x_container_[0]) {
            return y_container_[0] +
                   (y_container_[1] - y_container_[0]) / (x_container_[1] - x_container_[0]) *
                   (x - x_container_[0]);
        }
        // When x is bigger than x_container range.
        if (x > x_container_[max_iteration_]) {
            return y_container_[max_iteration_] +
                   (y_container_[max_iteration_] - y_container_[max_iteration_ - 1]) /
                   (x_container_[max_iteration_] - x_container_[max_iteration_ - 1]) *
                   (x - x_container_[max_iteration_]);
        }
        // search index which is most near x in container.
        const auto lower_iterator = x_container_.lower_bound(x);
        const auto lower_index    = std::distance(std::begin(x_container_), lower_iterator);
        return std::move(y_container_[lower_index] +
                         (y_container_[lower_index + 1] - y_container_[lower_index]) /
                         (x_container_[lower_index + 1] - x_container_[lower_index]) *
                         (x - x_container_[lower_index]));
    };
};
} /* namespace my_util */
#endif /* CALCULATION_UTIL_H_ */
