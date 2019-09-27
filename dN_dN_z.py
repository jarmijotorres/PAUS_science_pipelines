import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import time,sys ### for py3
from astropy.io import fits
from astropy.cosmology import WMAP7 as cosmo
from scipy.interpolate import interp1d
from scipy.optimize import newton,brentq#,ridder,bisect,brentq
sys.path.append('./dev_PAUS_science_pipelines/')
from binning_data import binning_function
from glob import glob

L_M = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LMi_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)

L_Z = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LZ_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)

L_Vgal = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LVgal_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)