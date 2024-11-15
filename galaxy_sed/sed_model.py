# from galaxy_sed.Dust_cpp_pybind.build.bin.sed_module import SED_calculator
from galaxy_sed.Dust_cpp_pybind.build.bin.sed_module import SED_calculator, Hyperparameters as CppHyperparameters
from galaxy_sed.hyperparams import Hyperparameters

class GalaxySED:
    def __init__(self, galaxy_age, galaxy_mass, starformation_timescale, gas_infall_timescale, dust_scale_height, galaxy_radius, n0_cnm, is_closedbox, imf_type):
        self.galaxy_age = galaxy_age
        self.galaxy_mass = galaxy_mass
        self.starformation_timescale = starformation_timescale
        self.gas_infall_timescale = gas_infall_timescale
        self.dust_scale_height = dust_scale_height
        self.galaxy_radius = galaxy_radius
        self.n0_cnm = n0_cnm
        self.is_closedbox = is_closedbox
        self.imf_type = imf_type
        # SED_calculatorの初期化など、必要な処理をここに記述できます

    def calculate_sed(self):
        # SED_calculatorを呼び出す際に、Hyperparametersオブジェクトを渡す
        # params = Hyperparameters(self.galaxy_age, self.galaxy_mass, self.starformation_timescale)
        # print("Created Hyperparameters:", params.galaxy_age, params.galaxy_mass, params.starformation_timescale)
        
        # C++のHyperparametersオブジェクトを作成
        params = CppHyperparameters()
        params.galaxy_age = int(self.galaxy_age)
        params.galaxy_mass = self.galaxy_mass
        params.starformation_timescale = self.starformation_timescale
        params.gas_infall_timescale = self.gas_infall_timescale
        params.dust_scale_height = self.dust_scale_height
        params.galaxy_radius = self.galaxy_radius
        params.n0_cnm = self.n0_cnm
        params.is_closedbox = self.is_closedbox
        params.imf_type = self.imf_type

        sed_calculator = SED_calculator(params)
        print("sed_calculator", sed_calculator)
        return sed_calculator