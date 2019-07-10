# coding: utf-8
# %load 1 2 6 12 13
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import interp1d
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='in')
plt.rc('ytick',direction='in')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='serif')
plt.rc('font',size=18)
fdata = fits.open('/cosma/home/durham/jarmijo/PAU_test/catalogs/W1_vipers.fits')
W1 = fdata[1].data
z = W1['zspec']
uband = W1['u_T07'][z>0.]; gband = W1['g_T07'][z>0.]
typ = W1['CE_PT'][z>0.] ; z = z[z>0]
c_ug = uband - gband#U - G
bl = 0.4
bu = 1.
Nb = 7
BINS = np.linspace(bl,bu,Nb)#np.linspace(0.01,0.95,5)
b = BINS[:-1] + np.diff(BINS)[0]/2.#bin centre
cmds = []
Ns = []
for i in range(4):
    xx = z[(typ == i+1)&(abs(c_ug)<2.5)]
    yy = c_ug[(typ == i+1)&(abs(c_ug)<2.5)]
    bmedian = np.zeros_like(b)
    N = np.zeros_like(b)
    for i in range(Nb-1):
        bmedian[i] = np.median(yy[(xx>BINS[i])&(xx<BINS[i+1])])
        N[i] = len(yy[(xx>BINS[i])&(xx<BINS[i+1])])
    cmds.append(bmedian)
    Ns.append(N)
limits = np.zeros((len(cmds)-1,len(cmds[0])),dtype=object)
for i in range(3):
    limits[i] = (cmds[i]+cmds[i+1])/2.
f12 = interp1d(b,limits[0],kind='linear',fill_value='extrapolate') 
f23 = interp1d(b,limits[1],kind='linear',fill_value='extrapolate')
f34 = interp1d(b,limits[2],kind='linear',fill_value='extrapolate')
def K_type(color,z):
    #color = u - g# by definition
    opts = np.array([4,3,2,1])
    bol = color<f34(z),f34(z) < color < f23(z),f23(z) < color < f12(z),color > f12(z)
    return opts[np.array(bol)]
