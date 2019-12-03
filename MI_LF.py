# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

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
# =============================================================================
# def Schecter_mag_log(M,ps,Ms,alpha):
#     return np.log10(0.4*np.log(10)*ps*(10**(0.4*(Ms-M)))**(alpha+1)*np.exp(-10**(0.4*(Ms-M))))
# 
# =============================================================================
# read data and define the K - correction
f1 = fits.open('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/catalogs/LCmock_radecz_MiCFHTLS_MB_SFR_23cut.fits')
#f1 = fits.open('/Users/Joaquin/Documents/Catalogs/mocks/mocks_radecz_MIMB_SFRHaplha_23cut_fix_nfrf.fits')
data1 = f1[1].data
dL = dLum(data1['Z'])
#Kfc = data1['SDSS_i'] - data1['SDSS_I'] -25. - 5*np.log10(dL) # aparent_magnitude -  absolute_magnitude
Kfc = data1['CFHTLS_i'] - data1['CFHTLS_I'] - 25. - 5*np.log10(dL)#CFHTLS photometry
# =============================================================================
# Colors based in u-g sdss magnitudes bins with the same number of points
# =============================================================================
cl = 0.25
cu = 1.75
Nc = 8
cbins = np.linspace(cl,cu,Nc+1,endpoint=True)
u_g = data1['SDSS_U'] - data1['SDSS_G']
Ncolor,edges_color = np.histogram(u_g,bins=40,range=(0.25,1.75),density=True)
Ncolor *= np.diff(edges_color)[0] # in this way the distribution sums 1
# cumulative distribution
cum_col = np.cumsum(Ncolor)
cb_color = edges_color[:-1] + np.diff(edges_color)[0]/2.
f_cumcolor = interp1d(cb_color,cum_col,kind='linear',fill_value='extrapolate') 
# new function domain
color_range = np.arange(cbins[0],cbins[-1],0.00001)
n_percentiles = np.linspace(0,1,Nc+1) # in octiles
ec_edges = np.zeros_like(n_percentiles)
ec_edges[0] = cl
K_color_bins = []
for i in range(Nc):
    R = np.isclose(f_cumcolor(color_range),n_percentiles[i+1],rtol=1e-4)
    ec_edges[i+1] = np.average(color_range[R])
    Ki = Kfc[(u_g > ec_edges[i])&(u_g < ec_edges[i+1])]
    K_color_bins.append(Ki)
#
# give to each K(z) function 20 points.
Ns = 20
zmin = 0.00 #range valid only for I-band magnitude
zmax = 1.2
zbins = np.linspace(zmin,zmax,Ns+1,endpoint=True)#### redshift bins
Kz_per_color = []
for i,Ki in enumerate(K_color_bins):
    Z_binned = data1['Z'][(u_g > ec_edges[i])&(u_g < ec_edges[i+1])]
    K_bin, _ = binning_function(Z_binned,Ki,zmin,zmax,Ns,percentile=16)
    Kz_per_color.append(K_bin)
# =============================================================================
Kz_per_color = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Kz_per_color_binning_8cbins_20zbins.npy')
# =============================================================================
Kz_functions = []
Kz_gen = []
for i,Kc in enumerate(Kz_per_color):
    xx = np.append(np.array([0]),Kc[:,0])
    yy = np.append(np.array([0]),Kc[:,3])
    ki = interp1d(xx,yy,kind='linear',fill_value='extrapolate') #K(z) per color interpolations
    Kz_functions.append(ki)
    Kz_gen.append(np.array([xx,yy]).T)
np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Kz_per_color_binning_8cbins_20zbins.npy',Kz_gen,allow_pickle=True)
K_binned, _ = binning_function(data1['Z'],Kfc,zmin,zmax,Ns,percentile=16)
# ====================================================
K_binned = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Kz_median_1cbin_20zbins.npy',allow_pickle=True)
# ================================================
xx = np.append(np.array([0]),K_binned[:,0])
yy = np.append(np.array([0]),K_binned[:,3]) # median all bins
Kz = interp1d(xx,yy,kind='linear',fill_value='extrapolate')#interpolation
Kz_median_gen = np.array([xx,yy]).T
np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Kz_median_1cbin_20zbins.npy',Kz_median_gen,allow_pickle=True)
  # =============================================================================
color_id = np.zeros_like(u_g,dtype=int)
for i,color in enumerate(u_g):
    for j in range(Nc):
        bc = (color > ec_edges[j]) and (color< ec_edges[j+1])
        if bc: idc = j
    color_id[i] = idc
Kz_gals = np.zeros_like(u_g,dtype=float)
for i,idc in enumerate(color_id):
    Kz_gals[i] = Kz_functions[idc](data1['Z'][i])
# Apply this K-correction to all abosultes magnitudes now on
#========== save new table with Kz corrections per each galaxy =============#
Kz_fits = fits.Column(name='Kz_gal',format='D',array=Kz_gals)
Cid_fits = fits.Column(name='Color_id',format='K',array=color_id)
orig_cols = data1.columns
new_cols = fits.ColDefs((Kz_fits,Cid_fits))
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/catalogs/LCmock_radecz_MiCFHTLS_SDSScolors_Kz_ColorID_MB_SFR_23cut.fits')
#
#
#
#
#=============== load from her if the K-correction was created already ====#
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
#M_new = data1['CFHTLS_i'] - 25. - 5*np.log10(dL) - data1['Kz_gal'] #+ 2.5*np.log10(1+data1['Z']) # New I-band 
bb = Mbins[:-1] + np.diff(Mbins)[0]/2.
######################## to get  Vmax #######################
micut=23
def get_zmax(Ms,Zs,color_id,zini,zend):
    X0 = micut - Ms - 25.
    fz = lambda x: DM(x) + Kz_functions[color_id](x) - X0
    #fz = lambda x: DM(x) - X0
#    root = ridder(fz,-0.001,zend)# Choose wisely
    #root = newton(fz,(zini+zend)*0.3,maxiter=int(1e3),tol=1.48e-05)
    root = brentq(fz,Zs,2.,maxiter=int(1e2),xtol=1.01e-05)
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
L_LF = []
L_M = []
#L_B = [] # this is the blue magnitude
L_Vmax = []
L_Vgal = []
Vi_Mbin = []
L_Z = []
t = time.process_time()
for z in range(Nz):
    zi = zbins[z]
    zf = zbins[z+1]
    N = M_new[(data1['Z']>zi)&(data1['Z']<zf)]
    #K = Kz_gals[(data1['Z']>zi)&(data1['Z']<zf)]
    C = data1['color_id'][(data1['Z']>zi)&(data1['Z']<zf)]
    #B = data1['MB'][(data1['Z']>zi)&(data1['Z']<zf)] # only available for the z = 0.11-0.9 catalog version
    Z = data1['Z'][(data1['Z']>zi)&(data1['Z']<zf)]
    Zmax = get_zmax_v(N,Z,C,zi,zf)
    Vmax = Omega_rad/3. * (dcomv(Zmax)**3. - dcomv(zi)**3)
    Vgal = Omega_rad/3. * (dcomv(Z)**3. - dcomv(zi)**3)
    # In each bin of the histogram find the minimum and maximum redshift 
# in order to compute the comoving volume surveyed by that bin.
    LF = np.zeros(b)
    for i in range(b):
        Vi = Vmax[(N > Mbins[i])&(N < Mbins[i+1])]
        M_in_bin = N[(N > Mbins[i])&(N < Mbins[i+1])]
        Vi_Mbin.append(Vi)
        LF[i] = np.sum(1./Vi)
#        if len(M_in_bin) > 0:
#            if ~np.isclose(M_in_bin.max(),Mbins[i+1],rtol=1e-3):
 #               LF[i] = 0.0
    L_LF.append(LF)
    L_M.append(N)
    L_Z.append(Z)
    L_Vmax.append(Vmax)
    L_Vgal.append(Vgal)
    #L_B.append(B)
e_t = time.process_time() - t 
#
L_LF = np.array(L_LF)
L_M = np.array(L_M)
L_Z = np.array(L_Z)
#L_B = np.array(L_B)
L_Vmax = np.array(L_Vmax)
L_Vgal = np.array(L_Vgal)
np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LLF_Mim'+str(abs(M_max))+'_m'+str(abs(M_min))+'_dM'+str(dM)+'_z'+str(zmin)+'_'+str(zmax)+'_'+str(Nz)+'zbins_Kcorr8cbins_mi23cut.npy',L_LF)

np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LZ_Mim'+str(abs(M_max))+'_m'+str(abs(M_min))+'_dM'+str(dM)+'_z'+str(zmin)+'_'+str(zmax)+'_'+str(Nz)+'zbins_Kcorr8cbins_mi23cut.npy',L_M,allow_pickle=True)

np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LMi_Mim'+str(abs(M_max))+'_m'+str(abs(M_min))+'_dM'+str(dM)+'_z'+str(zmin)+'_'+str(zmax)+'_'+str(Nz)+'zbins_Kcorr8cbins_mi23cut.npy',L_M,allow_pickle=True)

#np.save('/Users/jarmijo/Documents/Mocks/MB_8bins_z0.11_z.9_mi23cut.npy',L_B,allow_pickle=True)

np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LVmax_Mim'+str(abs(M_max))+'_m'+str(abs(M_min))+'_dM'+str(dM)+'_z'+str(zmin)+'_'+str(zmax)+'_'+str(Nz)+'zbins_Kcorr8cbins_mi23cut.npy',L_Vmax,allow_pickle=True)

np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LVgal_Mim'+str(abs(M_max))+'_m'+str(abs(M_min))+'_dM'+str(dM)+'_z'+str(zmin)+'_'+str(zmax)+'_'+str(Nz)+'zbins_Kcorr8cbins_mi23cut.npy',L_Vgal,allow_pickle=True)
################################################
