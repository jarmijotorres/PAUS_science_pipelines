
# coding: utf-8

# In[21]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.special import erf
from scipy.optimize import curve_fit
import sys,os,glob
from astropy.cosmology import WMAP9 as cosmo
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='in')
plt.rc('ytick',direction='in')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='sans-serif')
plt.rc('font',size=18)
fs = 20


# In[2]:


data = np.loadtxt('/cosma/home/durham/jarmijo/PAU_test/catalogs/mags+kcorr_866.dat')


# In[5]:


iauto = data[:,1]
uband = data[:,2]
gband = data[:,3]
zp = data[:,4]
zs = data[:,5]
types = data[:,6].astype(int)
Kcorr = data[:,7]


# In[8]:


bands = np.array(['NB455', 'NB465', 'NB475', 'NB485', 'NB495', 'NB505', 'NB515',
       'NB525', 'NB535', 'NB545', 'NB555', 'NB565', 'NB575', 'NB585',
       'NB595', 'NB605', 'NB615', 'NB625', 'NB635', 'NB645', 'NB655',
       'NB665', 'NB675', 'NB685', 'NB695', 'NB705', 'NB715', 'NB725',
       'NB735', 'NB745', 'NB755', 'NB765', 'NB775', 'NB785', 'NB795',
       'NB805', 'NB815', 'NB825', 'NB835', 'NB845'], dtype=object)
bdic = []
c = 4550
b = 455
for l in bands:
    bdic.append(['NB%d'%b, c])
    c += 100
    b +=10
bdc = dict(bdic)


# In[9]:


with open('/cosma/home/durham/jarmijo/PAU_test/fluxes/ID_narrow_fluxes.dat','r') as f:
    fluxes = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
fluxes = [x.strip() for x in fluxes]
with open('/cosma/home/durham/jarmijo/PAU_test/fluxes/ID_narrow_bands.dat','r') as f:
    band = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
band = [x.strip() for x in band]
with open('/cosma/home/durham/jarmijo/PAU_test/fluxes/ID_narrow_fluxerror.dat','r') as f:
    fluxerr = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
fluxerr = [x.strip() for x in fluxerr]


# In[10]:


spcs = []
for i in range(len(fluxes)):
    spcs.append(np.array(fluxes[i].split(' '),dtype=float))
#
spcs_err = []
for i in range(len(fluxerr)):
    spcs_err.append(np.array(fluxerr[i].split(' '),dtype=float))
#
filters_obj = []
for i in range(len(band)):
    filters_obj.append(np.array(band[i].split(' '),dtype=str))


# In[11]:


xs = [] #center band de-shifted
xl = [] #low limit
spcs2=[]
xh = [] # high limit
for i in range(len(filters_obj)):
    #if zs[i] != 0:
    x = np.zeros(len(filters_obj[i]))
    x1,x2 = np.zeros(len(filters_obj[i])),np.zeros(len(filters_obj[i]))
    for c in range(len(x)):
        x[c] = bdc[filters_obj[i][c]]/(1+zp[i])
        x1[c]= (bdc[filters_obj[i][c]] - 5.)/(1+zp[i])
        x2[c]= (bdc[filters_obj[i][c]] + 5.)/(1+zp[i])
    #spcs2.append(spcs[i])
    xs.append(x)
    xl.append(x1)
    xh.append(x2)


# In[12]:


P_B = []
#dxp = 10
#dxp2 = 0.5*dxp
xi = 4050
xf = 4450
for i in range(len(xs)):
    dx = np.diff(xs[i])[0]
    dx2 = 0.5*dx
    Ix = xs[i][(xs[i]>xi - dx2)&(xs[i]<xf + dx2)]
    IF = spcs[i][(xs[i]>xi - dx2)&(xs[i]<xf + dx2)]
    F_in = np.zeros_like(IF)
    for c in range(len(Ix)):
        F_in[c] = IF[c]*(1.98e-10/(Ix[c])#flux in cgs units considering the energy of every photon emitted.
        if c == 0:
            if (0<(Ix[c]- xi)<dx2):
                F_in[c] *= ((dx2 - abs(Ix[c]- xi))/dx)
            elif (0<(xi - Ix[c])<dx2):
                F_in[c] *= ((dx2 + abs(Ix[c]- xi))/dx)
        if c == len(Ix) - 1:
            if (0<(xf - Ix[c])<dx2):
                F_in[c] *= ((dx2 - abs(Ix[c]- xf))/dx)
            elif (0<(Ix[c] - xf) < dx2):
                F_in[c] *= ((dx2 + abs(Ix[c]- xf))/dx)
    P_B.append(np.sum(F_in*dx))


# In[13]:


PB = np.array(P_B)


# In[34]:


PAUS_BLUE = PB[(PB>0)&(zp<0.9)&(zp>0.11)]
zB = zp[(PB>0)&(zp<0.9)&(zp>0.11)]
zS = zs[(PB>0)&(zp<0.9)&(zp>0.11)]
imag = iauto[(PB>0)&(zp<0.9)&(zp>0.11)]
gmag = gband[(PB>0)&(zp<0.9)&(zp>0.11)]
umag = uband[(PB>0)&(zp<0.9)&(zp>0.11)]
kcorr = Kcorr[(PB>0)&(zp<0.9)&(zp>0.11)]
tp = types[(PB>0)&(zp<0.9)&(zp>0.11)]


# In[35]:


#PAUS_BLUE = 26. -2.5*np.log10(PAUS_BLUE)# 2 AB magnitude
mab_blue = -2.5*np.log10(PAUS_BLUE) - 48 #flux (in cgs) to AB mags.

# In[36]:


dL = cosmo.luminosity_distance(zB)



# In[37]:


Imag = imag - 25. - 5.*np.log10(dL.value) - kcorr# 2 Absolute Magnitude
Imag_nk = imag - 25. - 5.*np.log10(dL.value)# 2 Absolute Magnitude
MB =  mab_blue- 25. - 5.*np.log10(dL.value)# 2 Absolute Magnitude



S = np.array([Imag,MB,zB,zS]).T


# In[56]:

