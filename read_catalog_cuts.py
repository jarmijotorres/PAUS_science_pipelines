# coding: utf-8
# %load codes/read_catalog_cuts.py
# %load 1-15
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
from pandas import read_csv
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
# %load PAU_test/binning_data.py
import numpy as np    
plt.rc('font',family='sans-serif')
plt.rc('font',size=18)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
NPAU = np.loadtxt('/cosma/home/durham/jarmijo/PAU_test/catalogs/gals_imag_866.dat')
i_AB = NPAU[:,7]
NPAU = NPAU[i_AB<22.5]
iauto = NPAU[:,7]
IDs = NPAU[:,0]
zs = NPAU[:,3]
zp = NPAU[:,4]
zpw = NPAU[:,5]
zc = NPAU[:,6]
uband = NPAU[:,8]
gband = NPAU[:,9]
conf = NPAU[:,10]
odss = NPAU[:,11]
qz = NPAU[:,12]
chi2 = NPAU[:,13]
Nb = NPAU[:,14]
N_qz,e_qz = np.histogram(np.log10(qz),bins=40,density=True,range=(-1.5,2.5))
#N_qzspec, e_qzspec = np.histogram(np.log10(qzspec),bins=40,density=True,range=(-1.5,2.5))
Sqz = np.cumsum(N_qz)#; Sqzspec = np.cumsum(N_qzspec)
b1 = e_qz[:-1] + np.diff(e_qz)[0]/2.#; b2 = e_qzspec[:-1] + np.diff(e_qzspec)[0]/2.
B1=10**b1# linear space
#B2=10**b2 #
Sqz*=np.diff(e_qz)[0]#renormalization
#Sqzspec *= np.diff(e_qzspec)[0] #renorm.
S_int = interp1d(B1,Sqz,kind='linear',fill_value='extrapolate')### iterpolated function
X = np.linspace(B1[0],B1[-1],int(1e5))## new x axis
t = 0.0001
t2 = 1e-5
w50 = np.isclose(S_int(X),0.5,atol=t)
Nodds,eodds = np.histogram(np.log10(odss[odss!=0]),bins=40,density=True,range=(-1.,0.0))
Sodds= np.cumsum(Nodds[::-1])
b = eodds[:-1][::-1] + np.diff(e_qz)[0]/2.
B=10**b# linear space
Sodds*=np.diff(eodds)[0]#renormalization
Sodds_int = interp1d(B,Sodds,kind='linear',fill_value='extrapolate')### iterpolated function
#ax.plot(XO,Sodds_int(XO),'b.')
cflag1 = conf == 2.5
cflag2 = (conf >= 3.0)&(conf<5.0)
cflag3 = conf == 9.5
zp_f = zp[(cflag1|cflag2|cflag3)&(Nb>=39)]
zs_f = zs[(cflag1|cflag2|cflag3)&(Nb>=39)]
imag_f = iauto[(cflag1|cflag2|cflag3)&(Nb>=39)]
odd_f = odss[(cflag1|cflag2|cflag3)&(Nb>=39)]
Qz_f = qz[(cflag1|cflag2|cflag3)&(Nb>=39)]
chi2_f = chi2[(cflag1|cflag2|cflag3)&(Nb>=39)]
dL_p = cosmo.luminosity_distance(zp_f).value
dL_s = cosmo.luminosity_distance(zs_f).value
DDL= -5.*np.log10(dL_s/dL_p)
dz = (zp_f-zs_f)/(1.+zp_f)
zp_f_no = zp_f[(DDL>-2.5)&(DDL<2.5)]
zs_f_no = zs_f[(DDL>-2.5)&(DDL<2.5)]
dz_no = dz[(DDL>-2.5)&(DDL<2.5)]
odd_f_no = odd_f[(DDL>-2.5)&(DDL<2.5)]
Qz_f_no = Qz_f[(DDL>-2.5)&(DDL<2.5)]
chi2_f_no = chi2_f[(DDL>-2.5)&(DDL<2.5)]
