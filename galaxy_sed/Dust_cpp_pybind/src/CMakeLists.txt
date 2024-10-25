# ソースファイルをワイルドカードで指定
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# pybind11を見つける
find_package(pybind11 REQUIRED)

# C++フラグに-fPICを追加
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# ライブラリを作成
add_library(SED_calculator ${SOURCES})

# pybind11のインクルードディレクトリを追加
target_include_directories(SED_calculator PRIVATE ${pybind11_INCLUDE_DIRS})

# OpenMP を SED_calculator にリンク
if (OpenMP_CXX_FOUND)
    target_link_libraries(SED_calculator PUBLIC OpenMP::OpenMP_CXX)
endif()
