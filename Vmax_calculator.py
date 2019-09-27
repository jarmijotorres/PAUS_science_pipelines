import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import time,sys ### for py3
from astropy.io import fits
from astropy.cosmology import WMAP7 as cosmo
from scipy.interpolate import interp1d
from scipy.optimize import newton#,ridder,bisect,brentq
sys.path.append('./dev_PAUS_science_pipelines/')
from binning_data import binning_function
from glob import glob
# useful mathematicals definitions
h = cosmo.h
def dcomv(z):
    return cosmo.comoving_distance(z).value * h #in Mpc/h
def dLum(z):
    return cosmo.luminosity_distance(z).value * h #in Mpc/h
def DM(z):
    return 5*np.log10((dLum(z)))
def kz(z):
    return 2.5*np.log10(1+z)
# 
Kz_per_color = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Kz_per_color_binning_8cbins_20zbins.npy')
Kz_functions = []
for i,Kc in enumerate(Kz_per_color):
    ki = interp1d(Kc[:,0],Kc[:,1],kind='linear',fill_value='extrapolate') #K(z) per color interpolations
    Kz_functions.append(ki)
K_median = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Kz_median_1cbin_20zbins.npy',allow_pickle=True)
# ================================================
Kz = interp1d(K_median[:,0],K_median[:,1],kind='linear',fill_value='extrapolate')#interpolation
###
f2 = fits.open('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/catalogs/LCmock_radecz_MiCFHTLS_SDSScolors_Kz_ColorID_MB_SFR_23cut.fits')
data1 = f2[1].data
dL = dLum(data1['Z'])
# =============================================================================
# Absolute maginute I-band Luminosity function
Omega_deg = (data1['DEC'].max() - data1['DEC'].min())* (data1['RA'].max() - data1['RA'].min()) # cos(alpha) where alpha is an angle (?) 
Omega_rad = Omega_deg * (np.pi/180.)**2. # to sq. rad
b = 32 #N of bins LF
M_min = -16.
M_max = -24.
Mbins= np.linspace(M_max,M_min,b+1,endpoint=True)
dM = abs(M_max-M_min)/(b) # Nb = 32 between (-16,-24) --> dM = 0.25 (Nb = 40 --> dM = 0.2 for the same range)
#########################################################
M_new = data1['CFHTLS_i'] - 25. - 5*np.log10(dL) - data1['Kz_gal'] #+ 2.5*np.log10(1+data1['Z']) # New I-band 
bb = Mbins[:-1] + np.diff(Mbins)[0]/2.
######################## to get  Vmax #######################
micut=23
def get_zmax(Ms,Zi,color_id,zini,zend):
    X0 = micut - Ms - 25.
    fz = lambda x: DM(x) + Kz_functions[color_id](x) - X0
    #fz = lambda x: DM(x) - X0
#    root = ridder(fz,-0.001,zend)# Choose wisely
    root = newton(fz,Zi,maxiter=int(1e3),tol=1.48e-05)
    if root >= zend:
        return zend
    else:
        return root
#
get_zmax_v = np.vectorize(get_zmax)

zmin = 0.00
zmax = 1.2
Nz = 10 # Number of redshift slices to compute the luminosity function
zbins = np.linspace(zmin,zmax,Nz+1,endpoint=True)#### redshift bins
#zb = 0 #from z to 12
#zi = zbins[z]
#zf = zbins[z+1]
print('')
Zmax = get_zmax_v(M_new,data1['Z'],data1['color_id'],zmin,zmax)
Vmax = Omega_rad/3. * (dcomv(Zmax)**3. - dcomv(zmin)**3)
Vgal = Omega_rad/3. * (dcomv(data1['Z'])**3. - dcomv(zmin)**3)

Vm_fits = fits.Column(name='Vmax',format='D',array=Vmax)
Vg_fits = fits.Column(name='Vgal',format='D',array=Vgal)
orig_cols = data1.columns
new_cols = fits.ColDefs((Vm_fits,Vg_fits))
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/catalogs/LCmock_radecz_MiCFHTLS_SDSScolors_Kz_ColorID_MB_SFR_VmaxVgal_23cut.fits')

print('New catalogue created in /cosma5/data/dp004/dc-armi2/PAU/PAU_test/catalogs/LCmock_radecz_MiCFHTLS_SDSScolors_Kz_ColorID_MB_SFR_VmaxVgal_23cut.fits.\n')

print('end of program.\n')