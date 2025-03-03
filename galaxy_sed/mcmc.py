import numpy as np
from scipy.optimize import minimize
from galaxy_sed.sed_model import GalaxySED
from galaxy_sed.hyperparams import Hyperparameters
import emcee
import glob
from galaxy_sed.diagram import delete_past_data
import os
import matplotlib.pyplot as plt
import time
import multiprocessing


# 尤度関数の定義（上限データと非対称エラーバーに対応）
def log_likelihood(theta, observed_wavelength, y, error_plus, error_minus, upper_limit, params):
    galaxy_mass, galaxy_age, starformation_timescale, gas_infall_timescale, n0_cnm, sigma = theta
    
    delete_past_data()

    print("_______________________________________________")    
    print("galaxy mass: "+ str(galaxy_mass) +"e9 Msolar")
    print("galaxy age: "+ str((int(galaxy_age)+1)) + "00Myr")
    print("starformation timescale: " + str(starformation_timescale) + "Gyr")
    print("gas infall timescale: " + str(gas_infall_timescale) + "Gyr")
    print("n0_cnm: " + str(n0_cnm))
    
    # GalaxySED インスタンスを作成
    # ユーザが直接指定したパラメータを取りたいときはこれ.paramsの中にデータが入っている
    # sed_model = GalaxySED(galaxy_age, galaxy_mass, params.starformation_timescale, params.gas_infall_timescale, params.dust_scale_height, params.galaxy_radius, params.n0_cnm, params.is_infall, params.imf_type)
    
    # TODO: ユーザーから見てMCMCでもハイパーパラメータでも変数を指定できるようにしたい。現状はプログラム内で最初から決めている。
    # MCMCでパラメータをふったとき
    sed_model = GalaxySED(galaxy_age, galaxy_mass, starformation_timescale, gas_infall_timescale, params.dust_scale_height, params.galaxy_radius, n0_cnm, params.is_infall, params.imf_type)


    # SEDを計算する。計算結果はcalculation_resultに入る
    print("C++呼び出し")
    sed_model.calculate_sed()
    print("Python処理再開")
    
    formatted_galaxy_mass = f"{galaxy_mass:.2e}"  # 科学的記法にフォーマット
    # '+0' を削除して返す
    formatted_galaxy_mass = formatted_galaxy_mass.replace("e+", "e")
    # 'e' を基準に分割して基数と指数を取得
    base, exponent = formatted_galaxy_mass.split('e')
    base = float(base)
    exponent = int(exponent)

    # シフト処理
    new_exponent = exponent + 9

    # 再構成して返す
    formatted_galaxy_mass = f"{base:.2f}e{new_exponent}"
        
    # SEDファイルが生成されるのを待つ
    sed_file = None
    retry_count = 300  # 最大リトライ回数
    while retry_count > 0:
        print("Current working directory:", os.getcwd())
        search_string = 'galaxy_sed/Dust_cpp_pybind/calculation_result/makeSED/SED_M=' + str(formatted_galaxy_mass) + 'Msun_*.dat'
        sed_files = glob.glob(search_string)
        if sed_files:
            sed_file = sed_files[0]
            print("SEDファイルできました!: ", sed_file)
            break
        print("見つけたいsed_file: ", search_string)
        # print("sed_files見つからない: ", sed_files)
        time.sleep(1)  # 1秒待機してリトライ
        retry_count -= 1

    # ファイルが見つからなかった場合はエラーを出力
    if sed_file is None:
        raise FileNotFoundError("SEDファイルが見つかりませんでした。C++の計算が完了しているか確認してください。")
    
    # makeSEDフォルダ内の'SED'から始まるファイル名のdatファイルを読み込む
    # sed_file = glob.glob('calculation_result/makeSED/SED*.dat')[0]

    print("SED file: ", sed_file)
    calculation_result = np.loadtxt(sed_file, skiprows=1)

    
    # observed_wavelength(um)に最も近い値をcalculation_resultのWavelength(cm)列から持ってくる
    y_model = []
    for obs_wl in observed_wavelength:
        # Wavelength列（1列目）をcm単位からum単位に変換し、観測波長に最も近いインデックスを探す
        wl_column_um = calculation_result[:, 0] * 1e4  # cm to μm変換
        print("observed wavelength: ", obs_wl)
        closest_idx = np.abs(wl_column_um - obs_wl).argmin()
        print("closest_idx: ", closest_idx)
        
        # 対応するフラックス値（2列目）を取得してy_modelに追加
        y_model.append(calculation_result[closest_idx, 1])

    y_model = np.array(y_model)
    print("y_model: ", y_model)
    
    # 残差の計算
    residuals = y - y_model
    print("残差: ", residuals)
    
    # 非対称エラーバーを使用
    errors = np.where(residuals >= 0, error_minus, error_plus)

    # 通常データと上限データの尤度を計算
    log_likelihood_values = np.where(
        upper_limit,
        -0.5 * ((y_model - y) / error_plus) ** 2,  # 上限データ（観測値以下に制限）
        -0.5 * (residuals / errors) ** 2 - np.log(np.sqrt(2 * np.pi) * errors)  # 通常データ
    )
    print("尤度: ", log_likelihood_values)
    
    return np.sum(log_likelihood_values)

# プライオリティの定義
def log_prior(theta, param_limits):
    galaxy_mass, galaxy_age, starformation_timescale, gas_infall_timescale, n0_cnm, sigma = theta

    # param_limits から下限と上限を取得
    galaxy_mass_min = param_limits["galaxy_mass_min"]
    galaxy_mass_max = param_limits["galaxy_mass_max"]
    galaxy_age_min = param_limits["galaxy_age_min"]
    galaxy_age_max = param_limits["galaxy_age_max"]
    starformation_timescale_min = param_limits["starformation_timescale_min"]
    starformation_timescale_max = param_limits["starformation_timescale_max"]
    gas_infall_timescale_min = param_limits["gas_infall_timescale_min"]
    gas_infall_timescale_max = param_limits["gas_infall_timescale_max"]
    n0_cnm_min = param_limits["n0_cnm_min"]
    n0_cnm_max = param_limits["n0_cnm_max"]
    sigma_min = param_limits["sigma_min"]
    sigma_max = param_limits["sigma_max"]
    if galaxy_mass_min < galaxy_mass < galaxy_mass_max and galaxy_age_min < galaxy_age < galaxy_age_max and starformation_timescale_min < starformation_timescale < starformation_timescale_max and gas_infall_timescale_min < gas_infall_timescale < gas_infall_timescale_max and n0_cnm_min < n0_cnm < n0_cnm_max and sigma_min < sigma < sigma_max:
        return 0.0
    return -np.inf

# 統合関数の定義
def log_probability(theta, observed_wavelength, y, error_plus, error_minus, upper_limit, params, param_limits):
    lp = log_prior(theta, param_limits)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, observed_wavelength, y, error_plus, error_minus, upper_limit, params)


class MCMCOptimizer:
    @staticmethod
    async def run_mcmc(params: Hyperparameters, observationnal_data, n_walkers, n_dimensions, param_limits):
        # ユーザーがinputしたハイパーパラメータを読み込む
        # galaxy_age = params.galaxy_age
        # galaxy_mass = params.galaxy_mass
        # starformation_timescale = params.starformation_timescale
        # gas_infall_timescale = params.gas_infall_timescale
        # dust_scale_height = params.dust_scale_height
        # galaxy_radius = params.galaxy_radius
        # n0_cnm = params.n0_cnm
        # is_infall = params.is_infall
        # imf_type = params.imf_type

        if params.galaxy_age is None or params.galaxy_mass is None:
            raise ValueError("galaxy_age, galaxy_mass, and starformation_timescale must be provided in params.")

        # 観測データの読み込み
        observed_data = np.loadtxt(observationnal_data, skiprows=1)
        observed_wavelength = observed_data[:, 0]
        y_data = observed_data[:, 1]
        error_plus = observed_data[:, 2]
        error_minus = observed_data[:, 3]
        upper_limit = observed_data[:, 4].astype(bool)

        # 各パラメータの下限と上限を設定
        galaxy_mass_min, galaxy_mass_max = param_limits["galaxy_mass_min"], param_limits["galaxy_mass_max"]
        galaxy_age_min, galaxy_age_max = param_limits["galaxy_age_min"], param_limits["galaxy_age_max"]
        starformation_timescale_min, starformation_timescale_max = param_limits["starformation_timescale_min"], param_limits["starformation_timescale_max"]
        gas_infall_timescale_min, gas_infall_timescale_max = param_limits["gas_infall_timescale_min"], param_limits["gas_infall_timescale_max"]
        n0_cnm_min, n0_cnm_max = param_limits["n0_cnm_min"], param_limits["n0_cnm_max"]        
        sigma_min, sigma_max = param_limits["sigma_min"], param_limits["sigma_max"]

        # 初期位置を設定（下限と上限を考慮）
        initial_positions = np.empty((n_walkers, n_dimensions))
        initial_positions[:, 0] = np.random.rand(n_walkers) * (galaxy_mass_max - galaxy_mass_min) + galaxy_mass_min
        initial_positions[:, 1] = np.random.rand(n_walkers) * (galaxy_age_max - galaxy_age_min) + galaxy_age_min  # galaxy_ageを整数に丸める
        initial_positions[:, 2] = np.random.rand(n_walkers) * (starformation_timescale_max - starformation_timescale_min) + starformation_timescale_min
        initial_positions[:, 3] = np.random.rand(n_walkers) * (gas_infall_timescale_max - gas_infall_timescale_min) + gas_infall_timescale_min
        initial_positions[:, 4] = np.random.rand(n_walkers) * (n0_cnm_max - n0_cnm_min) + n0_cnm_min
        initial_positions[:, 5] = np.random.rand(n_walkers) * (sigma_max - sigma_min) + sigma_min
        
        # initial_positions[:, 0] = np.random.uniform(param_limits["galaxy_mass_min"], param_limits["galaxy_mass_max"], n_walkers)
        # initial_positions[:, 1] = np.random.randint(param_limits["galaxy_age_min"], param_limits["galaxy_age_max"] + 1, n_walkers)
        # initial_positions[:, 2] = np.random.uniform(param_limits["sigma_min"], param_limits["sigma_max"], n_walkers)
        # print("initial_positions",initial_positions)

        # MCMCの実行（並列)
        print("MCMCの実行")
        # プロセスプールを作成
        # with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        #     sampler = emcee.EnsembleSampler(
        #         n_walkers, 
        #         n_dimensions, 
        #         log_probability, 
        #         args=(observed_wavelength, y_data, error_plus, error_minus, upper_limit, params, param_limits), 
        #         pool=pool
        #     )
        #     # MCMCの実行
        #     # sampler.run_mcmc(initial_positions, 5000, progress=True)
        #     n_steps = 5000
        #     for step, result in enumerate(sampler.sample(initial_positions, iterations=n_steps, progress=True)):
        #         if step % 100 == 0:  # 100ステップごとに収束状況を表示
        #             samples = sampler.get_chain(flat=True, discard=step)
        #             if len(samples) > 0:
        #                 mean_values = np.mean(samples, axis=0)
        #                 std_values = np.std(samples, axis=0)
        #                 print(f"Step {step}: mean={mean_values}, std={std_values}")
        
        
        # 非並列処理
        sampler = emcee.EnsembleSampler(n_walkers, n_dimensions, log_probability, args=(observed_wavelength, y_data, error_plus, error_minus, upper_limit, params, param_limits))
        # sampler.run_mcmc(initial_positions, 5000, progress=True)
        
        n_steps = 5000
        for step, result in enumerate(sampler.sample(initial_positions, iterations=n_steps, progress=True)):
            if step % 2 == 0:  # 100ステップごとに収束状況を表示
                samples = sampler.get_chain(flat=True)
                if len(samples) > 0:
                    mean_values = np.mean(samples, axis=0)
                    std_values = np.std(samples, axis=0)
                    print(f"Step {step}: mean={mean_values}, std={std_values}")

        print("MCMCの終了")

        # 結果のプロット
        samples = sampler.get_chain(flat=True)
        
        # k, y0, sigmaのヒストグラム
        fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
        axes[0].hist(samples[:, 0], bins=50, color='blue', alpha=0.7, label='galaxy_mass')
        axes[0].set_ylabel('Frequency')
        axes[0].legend()

        axes[1].hist(samples[:, 1], bins=50, color='orange', alpha=0.7, label='galaxy_age')
        axes[1].set_ylabel('Frequency')
        axes[1].legend()

        axes[2].hist(samples[:, 2], bins=50, color='green', alpha=0.7, label='sigma')
        axes[2].set_ylabel('Frequency')
        axes[2].legend()
        axes[2].set_xlabel('Parameter Value')

        plt.tight_layout()
        plt.show()

        # 推定結果の表示
        k_estimate = np.mean(samples[:, 0])
        y0_estimate = np.mean(samples[:, 1])
        sigma_estimate = np.mean(samples[:, 2])
        print(f"Estimated k: {k_estimate}")
        print(f"Estimated y0: {y0_estimate}")
        print(f"Estimated sigma: {sigma_estimate}")
        
        # # GalaxySED インスタンスを作成
        # sed_model = GalaxySED(galaxy_age, galaxy_mass, starformation_timescale, gas_infall_timescale, dust_scale_height, galaxy_radius, n0_cnm, is_infall, imf_type)
        # sed_model.calculate_sed()
        
        """MCMCの代わりに簡単な最適化を実行（例えば、パラメータの平均を返す）"""
        return None


    # MCMC用の関数定義
    # def log_posterior(self, params):
    #     """事後分布の対数を計算する"""
    #     sed = self.sed_model.calculate_sed(params)
    #     # 評価する目的関数を定義（例：観測データとのフィッティング）
    #     return -np.sum((sed - observed_data)**2)
