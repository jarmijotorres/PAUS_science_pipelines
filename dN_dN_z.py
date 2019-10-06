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

L_LF = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LLF_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy')

L_M = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LMi_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)

L_Z = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LZ_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)

L_Vgal = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LVgal_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)

L_Vmax = np.load('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/LVmax_Mim24.0_m16.0_dM0.25_z0.0_1.2_10zbins_Kcorr8cbins_mi23cut.npy',allow_pickle=True)

f1 = fits.open('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/catalogs/LCmock_radecz_MiCFHTLS_SDSScolors_Kz_ColorID_MB_SFR_VmaxVgal_23cut.fits')

zmin = 0.00
zmax = 1.2
Nz = 10 # Number of redshift slices to compute the luminosity function
zbins = np.linspace(zmin,zmax,Nz+1,endpoint=True)#### redshift bins

data1 = f1[1].data
ra_dec = np.array([data1['RA'],data1['DEC']]).T
sv_id = np.zeros_like(ra_dec[:,0],dtype = int)
pid = np.zeros_like(ra_dec,dtype = int)
ramin = ra_dec[:,0].min()
ramax = ra_dec[:,0].max() + 0.01
decmin = ra_dec[:,1].min()
decmax = ra_dec[:,1].max() + 0.01
dRA2 = (ramax - ramin)/2.
dDec2 = (decmax - decmin)/2.
for i in range(len(ra_dec)):
    pid[i] = [int((ra_dec[i,0]-ramin)/dRA2),int((ra_dec[i,1]-decmin)/dDec2)]
    sv_id[i] = pid[i][0]*2 + pid[i][1]
#    
L_IDs = []
for z in range(Nz):
    zi = zbins[z]
    zf = zbins[z+1]
    ID = sv_id[(data1['Z']>zi)&(data1['Z']<zf)]
    L_IDs.append(ID)
    
#ra_rad = np.radians(ra_dec[:,0])
#f = plt.figure(figsize=(10,6))
#ax = f.add_subplot(111,projection='polar')
#ax.scatter(ra_rad,data1['Z'],c='k',s=0.01)
#ax.set_rlim(zbins[0],zbins[1])
#ax.set_thetalim(np.radians(ramin),np.radians(ramax-1))
#ax.grid(False)
#plt.savefig('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/figures/LC_zbin0_allradec.png')
#plt.show()

Omega_deg = (data1['DEC'].max() - data1['DEC'].min())* (data1['RA'].max() - data1['RA'].min()) # cos(alpha) where alpha is an angle (?) 
Omega_rad = Omega_deg * (np.pi/180.)**2. # to sq. rad
Nb = 32 #N of bins LF
M_min = -16.
M_max = -24.
Mbins= np.linspace(M_max,M_min,Nb+1,endpoint=True)
dM = abs(M_max-M_min)/(Nb) # Nb = 32 between (-16,-24) --> dM = 0.25 (Nb = 40 --> dM = 0.2 for the same range)
bb = Mbins[:-1] + np.diff(Mbins)[0]/2.
### LF in 4 chunks for 0.00 < z < 0.12
z = 1
N = L_M[z]
Z = L_Z[z]
ID = L_IDs[z]
Vmax = L_Vmax[z]
LF_chunks=[]
for c in range(4):
    LF = np.zeros(Nb)
    for i in range(Nb):
        Vi = Vmax[ID == c][(N[ID == c] > Mbins[i])&(N[ID == c] < Mbins[i+1])]
        LF[i] = np.sum(1./Vi)
    LF_chunks.append(LF)
    
f,ax = plt.subplots(1,1,figsize=(7,6))
for c in range(len(LF_chunks)):
    zi = zbins[c]
    zf = zbins[c+1]
    ax.plot(bb,np.log10(4*LF_chunks[c]),'o-',c=np.random.random(3),ms=5.,label = "subvol %d"%(c+1))
ax.plot(bb,np.log10(L_LF[1]),'o--',c='k',ms=5.,label = r'%.2f < z < %.2f'%(zbins[1],zbins[2]))
ax.legend(prop={'size':12})
# =============================================================================
ax.set_xticks(np.arange(-24,-15,1))
ax.set_xlim(-23,-16)
ax.set_yticks(np.arange(-7,-1,1))
ax.set_ylim(-6.,-1.8)
ax.legend(prop = {'size':12})
ax.set_xlabel('$M_{i} - 5\log_{10}h$')
ax.set_ylabel('$\log$ [$\Phi$ Mpc$^3$/$h^{-3}$ (0.25 mag)$^{-1}]$')
plt.savefig('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/figures/MiCFHTLS_LF_mocks_mi23cut_zbin1_maxcut_23_subvols.png',bbox_inches='tight')
plt.tight_layout()
plt.show()
    


Nz,ez= np.histogram(Z,bins=50,range=(zbins[z],zbins[z+1]))
zb = ez[:-1] + np.diff(ez)[0]/2.
dV_dz = Omega_rad/3. * (dcomv(zb)**3. - dcomv(zi)**3)
f,ax = plt.subplots(1,1,figsize=(7,6))
ax.plot(zb,nz,'r--')
for c in range(4):
    Nzc,_ = np.histogram(Z[ID==c],bins=50,range=(zbins[z],zbins[z+1]))
    nzc = Nzc / dV_dz / 4.
    ax.plot(zb,nzc,'o-',label='subvol %d'%(c+1))
ax.legend(loc=1)
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$n(z)$')
plt.savefig('',bbox_inches='tight')