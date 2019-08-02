# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import time ### for py3
from astropy.io import fits
from astropy.cosmology import WMAP7 as cosmo
from scipy.interpolate import interp1d
from scipy.optimize import ridder
from binning_data import binning_function
from glob import glob
h = cosmo.h


plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='inout')
plt.rc('ytick',direction='inout')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='sans-serif')
plt.rc('font',size=18)


def dcomv(z):
    return cosmo.comoving_distance(z).value * h #in Mpc/h
# =============================================================================
# def Schecter_mag_log(M,ps,Ms,alpha):
#     return np.log10(0.4*np.log(10)*ps*(10**(0.4*(Ms-M)))**(alpha+1)*np.exp(-10**(0.4*(Ms-M))))
# 
# =============================================================================
f1 = fits.open('/home/jarmijo/Documents/mocks/mocks_radecz_SDSSphotometry_SFR_23cut_z0.00_1.2.fits')
#f1 = fits.open('/Users/Joaquin/Documents/Catalogs/mocks/mocks_radecz_MIMB_SFRHaplha_23cut_fix_nfrf.fits')
data1 = f1[1].data
dL = cosmo.luminosity_distance(data1['Z'])
dL = dL.value
MI_noK = data1['SDSS_i'] - 25. - 5.*np.log10(dL) - 5*np.log10(h) # 5logh units if mock
Kfc = data1['SDSS_i'] - data1['SDSS_I'] -25. - 5*np.log10(dL) - 5*np.log10(h)
# =============================================================================
# Colors based in u-g sdss magnitudes bins with the same number of points
# =============================================================================
cl = 0.25
cu = 1.75
Nc = 8
cbins = np.linspace(cl,cu,Nc+1,endpoint=True)
u_g = data1['SDSS_U'] - data1['SDSS_G']
Ncolor,edges_color = np.histogram(u_g,bins=40,range=(0.25,1.75),density=True)
Ncolor *= np.diff(edges_color)[0]
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
    Ki = u_g[(u_g > ec_edges[i])&(u_g < ec_edges[i+1])]
    K_color_bins.append(Ki)
#
    
Ns = 16
zmin = 0.11
zmax = 0.9
zbins = np.linspace(zmin,zmax,Ns+1,endpoint=True)#### redshift bins
Kz_per_color = []
for i,Ki in enumerate(K_color_bins):
    Z_binned = data1['Z'][(u_g > ec_edges[i])&(u_g < ec_edges[i+1])]
    K_bin, _ = binning_function(Z_binned,Ki,zmin,zmax,Ns,percentile=16)
    Kz_per_color.append(K_bin)
# =============================================================================
Kz_functions = []
for i,Kc in enumerate(Kz_per_color):
    ki = interp1d(Kc[:,0],Kc[:,3],kind='linear',fill_value='extrapolate') #K(z) per color interpolations
    Kz_functions.append(ki)
     
K_binned, _ = binning_function(data1['Z'],Kfc,zmin,zmax,Ns,percentile=16)
Kmedian = K_binned[:,3] # median all bins
Kz = interp1d(K_binned[:,0],Kmedian,kind='linear',fill_value='extrapolate')#interpolation
 
  # =============================================================================
#   f,ax = plt.subplots(1,1,figsize=(7,6))
#   ax.scatter(data1['Z'],Kfc,s=4,c='k')
#   for i in range(len(Kz_per_color)):
#       ax.plot(Kz_per_color[i][:,0],Kz_per_color[i][:,3],'-')
#   plt.show()
#  # 
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
# =============================================================================
# 
# f,ax =plt.subplots(1,1,figsize=(7,6))
# ax.scatter(data1['Z'],MI_noK - data1['M_I'],s=5,c='b')
# plt.show()
# 
# =============================================================================

 # %load 92 107
Omega_deg = (data1['DEC'].max() - data1['DEC'].min())* (data1['RA'].max() - data1['RA'].min())
Omega_rad = Omega_deg * (np.pi/180.)**2.
Ns = 8
zmin = 0.0
zmax = 1.2
zbins = np.linspace(zmin,zmax,Ns+1,endpoint=True)#### redshift bins
b = 40 #N of bins LF
M_min = -16.
M_max = -23.5
Mbins= np.linspace(M_max,M_min,b+1,endpoint=True)
dM = abs(M_max-M_min)/(b)
#########################################################
Mmax = np.max(data1['SDSS_I'])# m_i = 23 i-band cut
zMmax = data1['Z'][np.where(data1['SDSS_I'] == Mmax ) ]
mMmax = Mmax + 25 + 5*np.log10(cosmo.luminosity_distance(zMmax).value) + 5*np.log10(h) #should be m_i = 23
zrange = np.arange(0.01,1.,0.01) # complete range of redshift
dlzrange = cosmo.luminosity_distance(zrange).value
Mrange = mMmax - 25 - 5*np.log10(dlzrange) - 5*np.log10(h)
fz = interp1d(zrange,Mrange,kind='cubic',fill_value='extrapolate')

M_new = data1['SDSS_i'] - 25. - 5*np.log10(dL) - 5*np.log10(h) - Kz_gals # K-corrected
######################## to get  Vmax #######################
def get_zmax2(M,zend):
    Mend = fz(zend)
    if M < Mend:
        r1 = zend
    if M > Mend:
        interp_fn2 = lambda x: fz(x) - M
        r1 = ridder(interp_fn2,-0.1,zend)
    return r1
########################################################
t = time.process_time()
get_zmax_v = np.vectorize(get_zmax2)
#Zmax = get_zmax_v(M_new)
elapsed_time = time.process_time() - t
#S = np.load('/Users/jarmijo/Documents/Mocks/mocks_zmax_Vmax.npy')
#Zmax = S[:,0]
#Vmax = S[:,1]
##############################################################
# %load 96
L_LF = []
L_bLF = []
L_M = []
L_B = []
L_Vmax = []
L_Vgal = []
for z in range(Ns-1):
    zi = zbins[z]
    zf = zbins[z+1]
    N = M_new[(data1['Z']>zi)&(data1['Z']<zf)]
    B = data1['MB'][(data1['Z']>zi)&(data1['Z']<zf)]
    Z = data1['Z'][(data1['Z']>zi)&(data1['Z']<zf)]
    Zmax = get_zmax_v(N,zf)
    V = Omega_rad/3. * (dcomv(Zmax)**3. - dcomv(zi)**3)
    Vgals = Omega_rad/3. * (dcomv(Z)**3. - dcomv(zi)**3)
    # In each bin of the histogram find the minimum and maximum redshift 
# in order to compute the comoving volume surveyed by that bin.
    LF = np.zeros(b)
    for i in range(b):
        Vi = V[(N > Mbins[i])&(N < Mbins[i+1])]
        LF[i] = np.sum(1./Vi)
    bb = Mbins[:-1] + np.diff(Mbins)[0]/2.
    L_LF.append(LF)
    L_bLF.append(bb)
    L_M.append(N)
    L_Vmax.append(V)
    L_Vgal.append(Vgals)
    L_B.append(B)
#
L_LF = np.array(L_LF)
L_bLF = np.array(L_bLF)
L_M = np.array(L_M)
L_B = np.array(L_B)
L_Vmax = np.array(L_Vmax)
L_Vgal = np.array(L_Vgal)
np.save('/Users/jarmijo/Documents/Mocks/LF_8bins_z0.11_z.9_mi23cut.npy',L_LF)

np.save('/Users/jarmijo/Documents/Mocks/Mi_8bins_z0.11_z.9_mi23cut.npy',L_M,allow_pickle=True)

np.save('/Users/jarmijo/Documents/Mocks/MB_8bins_z0.11_z.9_mi23cut.npy',L_B,allow_pickle=True)

np.save('/Users/jarmijo/Documents/Mocks/Vmax_8bins_z0.11_z.9_mi23cut.npy',L_Vmax,allow_pickle=True)

np.save('/Users/jarmijo/Documents/Mocks/Vgals_8bins_z0.11_z.9_mi23cut.npy',L_Vgal,allow_pickle=True)
################################################
nf = 2
nc = 4
f,ax = plt.subplots(nf,nc,figsize=(4*nf,4.),sharex=True,
                        sharey=True,gridspec_kw={'width_ratios': [1]*nc, 'height_ratios': [1]*nf})
for f in range(nf):
    for c in range(nc):
        Vr = L_Vgal[nc*f+c]/L_Vmax[nc*f+c]
        zi = zbins[nc*f+c]
        zf = zbins[nc*f + (c+1)]
        ax[f,c].hist(Vr,bins=20,range=(0,1),histtype='step',label = "%.2f < z < %.2f"%(zi,zf))
        ax[f,c].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
        ax[f,c].legend(prop={'size':10})
plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0.0)
plt.show()
############################
f,ax = plt.subplots(1,1,figsize=(7,6))
for c in range(len(L_M)):
    zi = zbins[c]
    zf = zbins[c+1]
    ax.hist(L_M[c],bins=30,range=(M_max,M_min),histtype='step',label = "%.2f < z < %.2f"%(zi,zf),linewidth=2.5)
ax.tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
ax.set_xlabel(r'$M_{i} - 5\log_{10}h$')
ax.set_ylabel('Counts')
ax.legend(prop = {'size':10},loc=2)
plt.tight_layout()
plt.show()
#==============================================
f,ax = plt.subplots(1,1,figsize=(7,6))
for c in range(len(L_LF)):
    zi = zbins[c]
    zf = zbins[c+1]
    ax.plot(L_bLF[c],np.log10(L_LF[c]),'o-',ms=5.,label = "%.2f < z < %.2f"%(zi,zf))
ax.legend(prop={'size':12})
# =============================================================================
ax.legend(prop = {'size':12})
ax.set_xlabel('$M_{i} - 5\log_{10}h$')
ax.set_ylabel('$\log\ [dN/V_{max}$ Mpc$^3$/$h^{-3}$ (0.25 mag)$^{-1}]$')
#plt.savefig('./Dropbox/PhD/Durham/Projects/PAU/figures/July/Mblue_LF_mocks_mi23cut_8zbins_maxcut_2.png',bbox_inches='tight')
plt.tight_layout()
plt.show()
# =============================================================================
lfs = glob('/')