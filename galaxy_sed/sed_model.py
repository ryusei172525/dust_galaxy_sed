from galaxy_sed.Dust_cpp_pybind.build.bin.sed_module import SED_calculator

class GalaxySED:
    def __init__(self, galaxy_age):
        self.galaxy_age = galaxy_age
        # SED_calculatorの初期化など、必要な処理をここに記述できます

    def calculate_sed(self):
        """SED計算のためのパラメータを設定し、C++で計算を行う"""
        sed_calculator = SED_calculator(self.galaxy_age)  # インスタンス変数を使用
        print("sed_calculator", sed_calculator)
        return sed_calculator