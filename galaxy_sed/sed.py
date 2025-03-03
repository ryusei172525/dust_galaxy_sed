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

class SEDModel:
    @staticmethod
    def calculate_sed(params: Hyperparameters):
        # ユーザーがinputしたハイパーパラメータを読み込む
        galaxy_age = params.galaxy_age
        galaxy_mass = params.galaxy_mass
        starformation_timescale = params.starformation_timescale
        gas_infall_timescale = params.gas_infall_timescale
        dust_scale_height = params.dust_scale_height
        galaxy_radius = params.galaxy_radius
        n0_cnm = params.n0_cnm
        is_infall = params.is_infall
        imf_type = params.imf_type

        if params.galaxy_age is None or params.galaxy_mass is None:
            raise ValueError("galaxy_age, galaxy_mass, and starformation_timescale must be provided in params.")

        # GalaxySED インスタンスを作成
        sed_model = GalaxySED(galaxy_age, galaxy_mass, starformation_timescale, gas_infall_timescale, dust_scale_height, galaxy_radius, n0_cnm, is_infall, imf_type)
        sed_model.calculate_sed()
        
        return None