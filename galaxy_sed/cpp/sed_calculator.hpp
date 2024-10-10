#ifndef SED_CALCULATOR_HPP
#define SED_CALCULATOR_HPP

#include <vector>

// クラス全体に対してシンボルの可視性を設定
class __attribute__((visibility("default"))) SEDCalculator {
public:
    SEDCalculator();
    std::vector<double> calculate(const std::vector<double>& params);
};

#endif
