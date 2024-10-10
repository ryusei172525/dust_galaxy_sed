# プロジェクト構成
galaxy_sed/
├── README.md
├── setup.py
├── galaxy_sed/
│   ├── __init__.py
│   ├── sed_model.py       # PythonからSED計算を呼び出す
│   ├── mcmc.py            # MCMCの実装
│   ├── utils.py           # 便利関数
│   ├── cpp/
│   │   ├── sed_calculator.cpp  # C++でのSED計算
│   │   ├── sed_calculator.hpp  # ヘッダファイル
│   │   ├── CMakeLists.txt  # C++コードのビルド設定
│   └── bindings/
│       ├── pybind_wrapper.cpp  # pybind11を使ったラッパー
│       └── CMakeLists.txt  # pybind11のビルド設定

# 必要なライブラリ
・pybind11: PythonとC++を連携するために使用。
・NumPy: 数値計算をサポート。
・scipy: MCMCのサンプル生成のために使用。
・cmake: C++部分のビルドのため。

# 開始方法
1. 最初にbuildする必要がある。
　`mkdir galaxy_sed/bindings/build`
　`cd galaxy_sed/bindings/build`
  `cmake ..`
  `make`

2. `make`コマンドが失敗せずにビルドが完了したなら`example.ipynb`ファイルが使えるようになる

