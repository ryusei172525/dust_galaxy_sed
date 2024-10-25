/*
 * container_util.h
 *
 *  Created on: 2018/06/22
 *      Author: NISHIDA kazuki
 *
 *      Container に関する便利な関数が詰まっている。
 *      テンプレートが多いのでヘッダのみ。
 */

#ifndef ARRAY_UTIL_H_
#define ARRAY_UTIL_H_

#include <algorithm>
#include <string>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

namespace my_util {
template<class InputType, class DelimiterType>
std::vector<InputType> Split(const InputType& input, const DelimiterType& delimiter) noexcept {
    // 今後、delimiter が見つかるたびに pos を見つかった位置に
    // 更新していき、pos 以降を探すようにする
    auto current_position = std::string::size_type(0);
    auto result           = std::vector<InputType>();

    while (current_position != std::string::npos) {
        const auto finded_position = input.find(delimiter, current_position);
        // 見つからなかった場合
        if (finded_position == InputType::npos) {
            result.emplace_back(input.substr(current_position));
            break;
        }
        // 見つかった場合
        result.emplace_back(input.substr(current_position, finded_position - current_position));
        current_position = finded_position + std::size(delimiter) - 1;
    }
    return result;
}

// vector などの成分を std::string に変換して delimiter でつなぐ
// 例： input = [.. .. Calculation file_name], delimiter = "/"
//     result = "../../Calculation/file_name"s
template<class Container, class DelimiterType>
std::string Combine(const Container& input, const DelimiterType& delimiter) noexcept {
    auto result = std::string();
    for (const auto& buffer : input) {
        result += static_cast<std::string>(buffer);
        if (buffer != *(std::end(input) - 1)) {
            result += delimiter;
        }
    }
    return result;
}

template<typename Container>
double SumIndexRange(size_t first_index, size_t last_index, const Container& container) {
    // コンテナの範囲外のインデックスを参照していないか確認
    const auto container_size = container.size();
    if (container_size < first_index || container_size < last_index) {
        std::cout << "Error: Container_size is " << container_size << ". But first_index is "
                  << first_index << " and last_index is " << last_index << std::endl;
        exit(EXIT_FAILURE);
    }
    const auto first_iterator = std::cbegin(container) + first_index;
    const auto last_iterator  = std::cbegin(container) + last_index;

    return std::move(std::accumulate(first_iterator, last_iterator, 0.0));
}

// 線形空間でもっとも近い値のイテレータを戻す
//     例: container = {1, 5, 50, 100}, value = 70 のときは、50 の成分を示すイテレータを戻す
// ２つの成分の中心の value のときは大きい方の成分の iterator を戻す
//     例: container = {1, 5, 50, 100}, value = 75 のときは、50 の成分を示すイテレータを戻す
// container が降順の場合は begin_iterator に rbegin、end_iterator に rend を与えれば正しく計算可能
template<class IteratorType, class ValueType>
IteratorType SearchContainerMostNearIterator(const IteratorType& begin_iterator,
                                             const IteratorType& end_iterator,
                                             const ValueType value) {
    // std::lower_bound は value が配列の最大値より大きいとき未定義になるのでここで例外処理が必要。
    // end_iterator は最後尾の一つ後ろを指しているので注意が必要
    if (value >= *(end_iterator - 1)) return end_iterator - 1;

    // value 以上になる最初のイテレータを戻す
    const auto iterator = std::lower_bound(begin_iterator, end_iterator, value);

    // これが成り立つときは次の iterator - 1 が未定義になる
    if (iterator == begin_iterator) return iterator;
    // ２つのイテレータの値のうち, value がより近い方を返す
    const auto median = (*(iterator - 1) + *iterator) * 0.5;
    return (value < median) ? iterator - 1 : iterator;
}

// Container が二重以上の構造になっているとエラーがでる。
template<class T>
void PrintContainer(const T& container, const std::string& message = "") {
    std::cout << message << "[";
    for (const auto& value : container) std::cout << value << " ";
    std::cout << "], size is " << container.size() << std::endl;
}

// コンテナと関数だけの引数で計算できるように rap した transform
template<typename Container, typename Func>
auto Transform(Container& container, const Func& function) {
    return std::transform(std::begin(container), std::end(container), std::begin(container),
                          function);
}

// web サイト [http://kou-yeung.hatenablog.com/entry/2013/12/19/233000] から取ってきた
// 今の所まったくしようしていないのでよくわからない
class Range {
  private:
    int index_;
    int max_;

  public:
    Range(int start, int count) : index_(start), max_(start + count) {
        if (count < 0) throw std::out_of_range("count is less than 0.");
    }

    // [0...count)
    Range(int count) : Range(0, count) {}

    // イテレータの条件を満たす
    int operator*() { return index_; }

    void operator++() { ++index_; }

    bool operator!=(Range&) { return index_ < max_; }

    // 範囲 for 文の条件を満たす
    Range begin() const { return *this; }

    Range end() const { return *this; }
};

} /* namespace my_util */

#endif /* ARRAY_UTIL_H_ */
