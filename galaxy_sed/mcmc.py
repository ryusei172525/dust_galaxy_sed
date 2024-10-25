import numpy as np
from scipy.optimize import minimize
from galaxy_sed.sed_model import GalaxySED

class MCMCOptimizer:
    async def run_mcmc(galaxy_age):
        # """MCMCで最適化を行う"""
        # result = minimize(self.log_posterior, initial_params)
        # return result.x
        sed_model = GalaxySED(galaxy_age)
        print("sed_model", sed_model)

        """MCMCの代わりに簡単な最適化を実行（例えば、パラメータの平均を返す）"""
        # best_params = np.mean(galaxy_age)  # 適当に平均を取る
        return sed_model.calculate_sed()

    # MCMC用の関数定義
    # def log_posterior(self, params):
    #     """事後分布の対数を計算する"""
    #     sed = self.sed_model.calculate_sed(params)
    #     # 評価する目的関数を定義（例：観測データとのフィッティング）
    #     return -np.sum((sed - observed_data)**2)
