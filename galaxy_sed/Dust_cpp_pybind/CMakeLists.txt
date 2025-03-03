# cmake_minimum_required(VERSION 3.11.4)
cmake_minimum_required(VERSION 3.13)  # CMakeのバージョンを指定（3.13以降を推奨）
cmake_policy(SET CMP0079 NEW)  # CMP0079ポリシーを有効化

project(test)

# CMakeのポリシーを設定 (CMP0074を使用する)
cmake_policy(SET CMP0074 NEW)

# コンパイラを指定
set(CMAKE_CXX_COMPILER g++)

# C++20を使用
set(CMAKE_CXX_STANDARD 20)

# コンパイルオプションを指定
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

# デフォルトのモデルディレクトリを設定
set(SED_MODEL_ROOT_DIR ${CMAKE_SOURCE_DIR}/)

# ヘッダーファイルにパスを渡す
add_definitions(-DSED_MODEL_ROOT_DIR="${SED_MODEL_ROOT_DIR}")

# インクルードディレクトリを指定
include_directories(${CMAKE_SOURCE_DIR}/include)

# 出力ディレクトリを設定
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# pybind11を見つける
find_package(pybind11 REQUIRED)

# サブディレクトリを追加
add_subdirectory(main)
add_subdirectory(src)

find_package(OpenMP REQUIRED)
# find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SED_calculator PUBLIC OpenMP::OpenMP_CXX)
endif()

# デバック用
set(CMAKE_BUILD_TYPE Debug) # デバッグビルドを有効化
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

