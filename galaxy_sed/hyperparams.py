# hyperparams.py

# デフォルト値の設定
default_galaxy_age = 129  # 例: 130Myr
default_galaxy_mass = 100  # 1.0e11
default_starformation_timescale = 3.9 # 3.9Gyr
default_gas_infall_timescale = 15 # 15Gyr
default_dust_scale_height = 150 # 150pc
default_galaxy_radius = 10000 # 10kpc
default_n0_cnm = 1000
default_is_infall = False
default_imf_type = 2 # 0:Salpeter IMF, 1:Larson IMF, 2:Chabrier IMF, 3:Chabrier_1+x=2, 4:Chabrier_1+x=1.5, 5:Chabrier_1+x=1.1 7:Chabrier_1+x=0.9

class Hyperparameters:
    def __init__(self, galaxy_age=None, galaxy_mass=None, starformation_timescale=None, gas_infall_timescale=None, dust_scale_height=None, galaxy_radius=None, n0_cnm=None, is_infall=None, imf_type=None):
        self.galaxy_age = galaxy_age if galaxy_age is not None else default_galaxy_age
        self.galaxy_mass = galaxy_mass if galaxy_mass is not None else default_galaxy_mass
        self.starformation_timescale = starformation_timescale if starformation_timescale is not None else default_starformation_timescale
        self.gas_infall_timescale = gas_infall_timescale if gas_infall_timescale is not None else default_gas_infall_timescale
        self.dust_scale_height = dust_scale_height if dust_scale_height is not None else default_dust_scale_height
        self.galaxy_radius = galaxy_radius if galaxy_radius is not None else default_galaxy_radius
        self.n0_cnm = n0_cnm if n0_cnm is not None else default_n0_cnm
        self.is_infall = is_infall if is_infall is not None else default_is_infall
        self.imf_type = imf_type if imf_type is not None else default_imf_type

# ハイパーパラメータを辞書として作成する関数
def get_hyperparams(galaxy_age=None, galaxy_mass=None, starformation_timescale=None, gas_infall_timescale=None, dust_scale_height=None, galaxy_radius=None, n0_cnm=None, is_infall=None, imf_type=None):
    return Hyperparameters(galaxy_age, galaxy_mass, starformation_timescale, gas_infall_timescale, dust_scale_height, galaxy_radius, n0_cnm, is_infall, imf_type)