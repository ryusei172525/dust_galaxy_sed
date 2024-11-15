# ソースファイルをワイルドカードで指定
file(GLOB MAIN_SOURCES "${CMAKE_SOURCE_DIR}/main/*.cpp")

# Python 3.11を指定
# find_package(PythonInterp 3.11 REQUIRED)
# find_package(PythonLibs 3.11 REQUIRED)
find_package(Python 3.11 REQUIRED COMPONENTS Development)

# pybind11を見つける
find_package(pybind11 REQUIRED)

# SED_calculatorライブラリのパスを指定
set(SED_CALCULATOR_LIB_PATH "${CMAKE_SOURCE_DIR}/build/src")
link_directories(${SED_CALCULATOR_LIB_PATH})

# pybind11モジュールを作成
pybind11_add_module(sed_module ${MAIN_SOURCES})

# 出力ディレクトリを設定
set_target_properties(sed_module PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# 実行ファイルを作成 (必要に応じて)
# add_executable(out ${MAIN_SOURCES})

# 出力ディレクトリを設定
# set_target_properties(out PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# 出力ディレクトリを設定
set_target_properties(sed_module PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# ライブラリへのリンクを指定
# target_link_libraries(out PRIVATE test_lib)
target_link_libraries(sed_module PRIVATE SED_calculator)

