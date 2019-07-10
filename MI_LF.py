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
f1 = fits.open('/Users/jarmijo/Documents/Mocks/mocks_radecz_MIMB_SFRHaplha_23cut_fix_nofrf.fits')
f1 = fits.open('/Users/Joaquin/Documents/Catalogs/mocks/mocks_radecz_MIMB_SFRHaplha_23cut_fix_nfrf.fits')
data1 = f1[1].data
dL = cosmo.luminosity_distance(data1['Z'])
dL = dL.value
MI_noK = data1['m_i'] - 25. - 5.*np.log10(dL) - 5*np.log10(h) # 5logh units if mock
Kfc = data1['m_i'] - data1['M_I'] -25. - 5*np.log10(dL) - 5*np.log10(h)
K_binned = binning_function(data1['Z'],Kfc,percentile=16)
Kmedian = K_binned[:,3]
Kz = interp1d(K_binned[:,0],Kmedian,kind='linear',fill_value='extrapolate')
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
Ns = 9
zmin = 0.11
zmax = 0.9
zbins = np.linspace(zmin,zmax,Ns,endpoint=True)#### redshift bins
b = 31 #N of bins LF
M_min = -16.
M_max = -23.5
Mbins= np.linspace(M_max,M_min,b,endpoint=True)
dM = abs(M_max-M_min)/(b-1)
#########################################################
Mmax = np.max(data1['M_I'])# m_i = 23 i-band cut
zMmax = data1['Z'][np.where(data1['M_I'] == Mmax ) ]
mMmax = Mmax + 25 + 5*np.log10(cosmo.luminosity_distance(zMmax).value) + 5*np.log10(h) #should be m_i = 23
zrange = np.arange(0.01,1.,0.01) # complete range of redshift
dlzrange = cosmo.luminosity_distance(zrange).value
Mrange = mMmax - 25 - 5*np.log10(dlzrange) - 5*np.log10(h)
fz = interp1d(zrange,Mrange,kind='cubic',fill_value='extrapolate')

M_new = data1['m_i'] - 25. - 5*np.log10(dL) - 5*np.log10(h) - Kz(data1['Z']) # K-corrected
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
get_zmax_v = np.vectorize(get_zmax2)
#t = time.process_time()
#Zmax = get_zmax_v(M_new)
#elapsed_time = time.process_time() - t
#S = np.load('/Users/jarmijo/Documents/Mocks/mocks_zmax_Vmax.npy')
#Zmax = S[:,0]
#Vmax = S[:,1]
##############################################################
Ns = 5
zmin = 0.11
zmax = 0.9
zbins = np.linspace(zmin,zmax,Ns,endpoint=True)#### redshift bins
# %load 96
L_LF = []
L_bLF = []
L_M = []
for z in range(Ns-1):
    zi = zbins[z]
    zf = zbins[z+1]
    N = M_new[(data1['Z']>zi)&(data1['Z']<zf)]
    Zmax = get_zmax_v(N,zf)
    V = Omega_rad/3. * ()
    # In each bin of the histogram find the minimum and maximum redshift 
# in order to compute the comoving volume surveyed by that bin.
    LF = np.zeros(b-1)
    for i in range(b-1):
        Vi = V[(N > Mbins[i])&(N < Mbins[i+1])]
        LF[i] = np.sum(1./Vi)
    bb = Mbins[:-1] + np.diff(Mbins)[0]/2.
    L_LF.append(LF)
    L_bLF.append(bb)
    L_M.append(N)
    
################################################
nf = 2
nc = 2
f,ax = plt.subplots(nf,nc,figsize=(4*nf,4.),sharex=True,
                        sharey=True,gridspec_kw={'width_ratios': [1]*nc, 'height_ratios': [1]*nf})
for f in range(nf):
    for c in range(nc):
        zi = zbins[nf*f+c]
        zf = zbins[nf*f + (c+1)]
        ax[f,c].hist(L_M[nf*f+c],bins=30,range=(M_max,M_min),histtype='step',label = "%.2f < z < %.2f"%(zi,zf))
        ax[f,c].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
        ax[f,c].legend()
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
ax.legend(prop = {'size':12})
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
