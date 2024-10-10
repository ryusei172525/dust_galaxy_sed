cmake_minimum_required(VERSION 3.10)
project(pybind_wrapper)

# pybind11のパッケージを探す
find_package(pybind11 REQUIRED)

# C++のサブディレクトリを追加
add_subdirectory(../cpp sed_calculator)

# ヘッダーファイルがあるディレクトリをインクルードパスに追加
include_directories(../cpp)

# pybind11のモジュールをビルド
pybind11_add_module(pybind_wrapper pybind_wrapper.cpp)

# SEDCalculatorライブラリをリンクする
target_link_libraries(pybind_wrapper PRIVATE sed_calculator)
