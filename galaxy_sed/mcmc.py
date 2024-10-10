import numpy as np
from scipy.optimize import minimize
from .sed_model import GalaxySED

class MCMCOptimizer:
    def __init__(self, age_range, mass_range):
        self.sed_model = GalaxySED(age_range, mass_range)

    def run_mcmc(self, initial_params):
        # """MCMCで最適化を行う"""
        # result = minimize(self.log_posterior, initial_params)
        # return result.x

        """MCMCの代わりに簡単な最適化を実行（例えば、パラメータの平均を返す）"""
        best_params = np.mean(initial_params)  # 適当に平均を取る
        return self.sed_model.calculate_sed([best_params])

    # MCMC用の関数定義
    # def log_posterior(self, params):
    #     """事後分布の対数を計算する"""
    #     sed = self.sed_model.calculate_sed(params)
    #     # 評価する目的関数を定義（例：観測データとのフィッティング）
    #     return -np.sum((sed - observed_data)**2)
