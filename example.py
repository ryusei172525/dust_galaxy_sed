import nest_asyncio

import sys
import asyncio
from galaxy_sed.mcmc import MCMCOptimizer
from galaxy_sed.diagram import delete_past_data
from galaxy_sed.hyperparams import get_hyperparams

sys.path.append('/home/ryuseikano/Study/dust_galaxy_sed')

# ネストされたイベントループを許可
nest_asyncio.apply()

galaxy_age = 4  # 4: 500Myr, 5: 600Myr, 129: 130Myr=13Gyr
galaxy_mass = 100  # 1: 1e9, 100:1e11
starformation_timescale_Gyr = 3.9 # 1: 1Gyr

# 初期値の設定
n_walkers = 50
n_dimensions = 6

# 各パラメータの下限と上限を設定
param_limits = {
    "galaxy_mass_min": 1, "galaxy_mass_max": 700,
    "galaxy_age_min": 0, "galaxy_age_max": 10,
    "starformation_timescale_min": 0.01, "starformation_timescale_max": 5,
    "gas_infall_timescale_min": 0.01, "gas_infall_timescale_max": 15,
    "n0_cnm_min": 1, "n0_cnm_max": 1000,
    "sigma_min": 0.01, "sigma_max": 1.0
}

async def main():
    # delete_past_data()
    # ハイパーパラメータを取得
    hyperparams = get_hyperparams(galaxy_age=galaxy_age, galaxy_mass=galaxy_mass, starformation_timescale=starformation_timescale_Gyr)  # galaxy_age のみを指定
    observational_data = 'observation/MACS0416Y1.dat'

    result = await MCMCOptimizer.run_mcmc(hyperparams, observational_data, n_walkers, n_dimensions, param_limits)  # 辞書を渡す

    print("計算終了")

# イベントループを実行
if __name__ == "__main__":
    asyncio.run(main())  # main関数を非同期に実行
