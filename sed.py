from galaxy_sed.sed import SEDModel
from galaxy_sed.diagram import delete_past_data
from galaxy_sed.hyperparams import get_hyperparams

galaxy_age = 5  # 4: 500Myr, 5: 600Myr, 129: 130Myr=13Gyr
galaxy_mass = 100  # 1: 1e9, 30: 3e10, 100:1e11
starformation_timescale_Gyr = 0.1 # 1: 1Gyr
gas_infall_timescale = 0.5
is_infall = True
n0_cnm = 3e5

delete_past_data()

# ハイパーパラメータを取得
hyperparams = get_hyperparams(galaxy_age=galaxy_age, galaxy_mass=galaxy_mass, starformation_timescale=starformation_timescale_Gyr, gas_infall_timescale=gas_infall_timescale, is_infall=is_infall, n0_cnm=n0_cnm)  # galaxy_age のみを指定
SEDModel.calculate_sed(hyperparams) 
# メモリの解放?

# galaxy_age = 3  # 4: 500Myr, 5: 600Myr, 129: 130Myr=13Gyr
# galaxy_mass = 200  # 1: 1e9, 100:1e11
# starformation_timescale_Gyr = 1 # 1: 1Gyr

# hyperparams = get_hyperparams(galaxy_age=galaxy_age, galaxy_mass=galaxy_mass, starformation_timescale=starformation_timescale_Gyr)  # galaxy_age のみを指定
# SEDModel.calculate_sed(hyperparams) 