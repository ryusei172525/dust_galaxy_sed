from .bindings.build.pybind_wrapper import SEDCalculator

class GalaxySED:
    def __init__(self, age_range, mass_range):
        self.age_range = age_range
        self.mass_range = mass_range
        self.sed_calculator = SEDCalculator()

    def calculate_sed(self, params):
        """SED計算のためのパラメータを設定し、C++で計算を行う"""
        return self.sed_calculator.calculate(params)