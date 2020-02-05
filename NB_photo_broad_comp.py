import numpy as np
from mpi4py import MPI
import h5py,time
from glob import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import newton,brentq#,ridder,bisect,brentq
from scipy.integrate import quad,trapz,simps
from astroML.stats import binned_statistic_2d

### open snapshot

subvols_dir = np.sort(glob('/data/dega1/ddmn39/galform_out/r504/lightcone_input/Gonzalez13.PAU.MillGas/ivol_*'))

S_iband = np.zeros(0,dtype=object)
#S_rband = np.zeros(0,dtype=object)

for sn,sv in enumerate(subvols_dir):
    
    spName = sv + '/galaxies.hdf5'
    sp_ivol = h5py.File(spName,'r')
    #O_x = np.array(sp_ivol.get("/Output00"+str(i+8)+"/xgal"))
    #O_y = np.array(sp_ivol.get("/Output00"+str(i+8)+"/ygal"))
    #O_z = np.array(sp_ivol.get("/Output00"+str(i+8)+"/zgal"))
    ##### test i-band at z=0 snapshot
    Mi_res = np.array(sp_ivol.get("/Output001/magMir_tot_ext"))
    #Mr_res = np.array(sp_ivol.get("/Output001/magrSr_tot_ext"))
    if Mi_res.any() == None: 
        print('M_i not available on subvol %d\n'%sn)
        #break
    elif Mi_res.any() != None: 
    #
        for nb in range(10,50): #i-band range defined by the broad i-band
            NB_i = np.array(sp_ivol.get("/Output001/mag%dr_tot_ext"%nb))
            Mi_res = np.vstack([Mi_res,NB_i])
        if sn == 0: 
            S_iband = Mi_res.T
        else:
            S_iband = np.append(S_iband,Mi_res.T,axis=0)

Mi = S_iband[:,0]
NBs = S_iband[:,1:]

#save data per band (broad band definition)
Ms_table.write('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Mibroad_allNBs_z_0.hdf5',format='hdf5',path='data',compression=True)
#load data if saved 
Sh5 = h5py.File('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/Mibroad_allNBs_z_0.hdf5')
c = np.array(Sh5['data'])
Mi_broad = c['Mi']

NB_mags = []
for i in range(22,41):
    NB_mags.append(c['NB_%d'%i])
NB_mags = np.array(NB_mags).T

#load PAUS narrow band responses
l = np.sort(glob('/cosma/home/dp004/dc-armi2/PAU_test/filters/PAUS_*_band_nt.dat'))
NB_response = []
for i in l:
    NB_r = np.loadtxt(i)
    NB_r[:,1] #normalized between 0 and 1
    NB_response.append(NB_r)
#load broad band responses
CFHT_Mi_response = np.loadtxt('/cosma/home/dp004/dc-armi2/PAU_test/filters/CFHT_MegaCam-iband_response.dat')
CFHT_Mr_response = np.loadtxt('/cosma/home/dp004/dc-armi2/PAU_test/filters/CFHT_MegaCam-rband_response.dat')

#filter response function interpolation
NB_functions = []
for i in range(40):
    f = interp1d(NB_response[i][:,0],NB_response[i][:,1],kind='linear')
    NB_functions.append(f)
Mi_broad_function = interp1d(CFHT_Mi_response[:,0],CFHT_Mi_response[:,1],kind='linear')
#narrow band integration
#weight calculation (ratio between NB peak and BB value)
c_nm = np.zeros(40)
for i in range(40):
    c_nm[i] = np.trunc(np.round(np.average(NB_response[i][:,0],weights=NB_response[i][:,1])/10.)/10.)
c_nm*=100
c_nm+=50 #locate in the band centre
#for the i-band the bands goes from
f_NBs_Mi = NB_response[22:41]
wfi = np.zeros(len(f_NBs_Mi))
for i in range(len(wfi)):
    wfi[i] = Mi_broad_function(c_nm[i+22])/NB_functions[i+22](c_nm[i+22])
#


#find intersections between narrowbands

L = f_NBs_Mi
roots = []
for i in range(len(L)-1):
    f0 = L[i]
    f1 = L[i+1]
    xcom = np.sort(list(set(f0[:,0][f0[:,1]>1e-2])&set(f1[:,0][f1[:,1]>1e-2])))
    x_xcom,x_xcom_1 = [],[]
    for nx in xcom:
        xs = np.where(f0[:,0]==nx)[0]
        x_xcom.append(xs[0])
        xs = np.where(f1[:,0]==nx)[0]
        x_xcom_1.append(xs[0])
    x_xcom = np.array(x_xcom)
    x_xcom_1 = np.array(x_xcom_1)
    h = np.zeros_like(f0[x_xcom])
    h[:,0] = f0[:,0][x_xcom]
    h[:,1] = f0[:,1][x_xcom] - f1[:,1][x_xcom_1]
    h_x = interp1d(h[:,0],h[:,1],kind='linear')
    roots.append(newton(h_x,x0=np.mean(h[:,0]),tol=1e-6))
rs = np.array(roots)

##### NB integration considering weights and overlap
wfi2 = np.zeros(len(f_NBs_Mi)) #including NB range in the i-band window
NB_area = np.zeros(len(f_NBs_Mi))
for i in range(len(wfi2)):
    nb = f_NBs_Mi[i]
    dl = np.diff(nb[:,0])[0]
    NB_area[i] = trapz(nb[:,1],dx=dl)
    if i == 0:
        Fa = trapz(nb[:,1],dx=dl)
        Fb = trapz(nb[:,1][nb[:,0]<rs[i]],dx=dl)
        c = Fb/Fa
        wfi2[i] = c
    elif i == len(f_NBs_Mi)-1:
        Fa = trapz(nb[:,1],dx=dl)
        Fb = trapz(nb[:,1][nb[:,0]>rs[i-1]],dx=dl)
        c = Fb/Fa
        wfi2[i] = c
    else:
        Fa = trapz(nb[:,1],dx=dl)
        Fb = trapz(nb[:,1][(nb[:,0]<rs[i])&(nb[:,0]>rs[i-1])],dx=dl)
        c = Fb/Fa
        wfi2[i] = c


Mi_narrow = np.zeros(len(NB_mags))
for i in range(len(Mi_narrow)):
    flux = 10**((NB_mags[i]+48.6)/-2.5)
    flux *= (NB_area/Mi_broad_area)
    F = np.dot(flux,(wfi*wfi2))#including 2 weigths
    Mi_narrow[i] =  -2.5*np.log10(F) - 48.6



##### Plot #######
cb =np.random.random((len(f_NBs_Mi),3))
f,ax = plt.subplots(1,1,figsize=(12,6))
for i in range(len(f_NBs_Mi)):
    nb = f_NBs_Mi[i]
    ax.plot(nb[:,0],nb[:,1],c=cb[i])
    if i == 0:
        ax.fill_between(nb[:,0][nb[:,0]<rs[i]],0,nb[:,1][nb[:,0]<rs[i]],facecolor=cb[i],alpha=0.3)
    elif i == len(f_NBs_Mi)-1:
        ax.fill_between(nb[:,0][nb[:,0]>rs[i-1]],0,nb[:,1][nb[:,0]>rs[i-1]],facecolor=cb[i],alpha=0.3)
    else:
        ax.fill_between(nb[:,0][(nb[:,0]<rs[i])&(nb[:,0]>rs[i-1])],0,nb[:,1][(nb[:,0]<rs[i])&(nb[:,0]>rs[i-1])],facecolor=cb[i],alpha=0.3)
ax.plot(CFHT_Mi_response[:,0],CFHT_Mi_response[:,1],c='darkred')
    
ax.set_xlabel(r'$\lambda$ [\AA]')
ax.set_ylabel(r'$R(\lambda)$')
ax.set_ylim(0.0,1.01)
ax.set_xlim(6600,8700)
#plt.savefig('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/figures/January2020/NBMi_CFHTMi_integration_nooverlap_wa.png')
plt.show()
#

f,ax = plt.subplots(1,1,figsize=(12,6))
ax.plot(CFHT_Mi_response[:,0],CFHT_Mi_response[:,1],c='darkred')
for i in range(len(f_NBs_Mi)):
    nb = f_NBs_Mi[i]
    ax.plot(nb[:,0],nb[:,1]*wfi[i])
ax.set_xlim(6600,8700)
ax.set_xlabel(r'$\lambda$ [\AA]')
ax.set_ylabel(r'$R(\lambda)$')
ax.set_title('Same effective area')
plt.show()  


##########################

###### in r-band l=555nm and l=695nm in OF

for sn,sv in enumerate(subvols_dir):
    spName = sv + '/galaxies.hdf5'
    sp_ivol = h5py.File(spName,'r')
    #O_x = np.array(sp_ivol.get("/Output00"+str(i+8)+"/xgal"))
    #O_y = np.array(sp_ivol.get("/Output00"+str(i+8)+"/ygal"))
    #O_z = np.array(sp_ivol.get("/Output00"+str(i+8)+"/zgal"))
    ##### test i-band at z=0 snapshot
    Mi_res = np.array(sp_ivol.get("/Output001/magMir_tot_ext"))
    Mr_res = np.array(sp_ivol.get("/Output001/magrSr_tot_ext"))
    
    if Mi_res.any() == None: 
        print('M_i not available on subvol %d\n'%sn)
        #break
    elif Mi_res.any() != None: 
    #
        for nb in range(35,50):
            NB_i = np.array(sp_ivol.get("/Output001/mag%dr_tot_ext"%nb))
            Mi_res = np.vstack([Mi_res,NB_i])
        if sn == 0: 
            S_iband = Mi_res.T
        else:
            S_iband = np.append(S_iband,Mi_res.T,axis=0)

Mi = S_iband[:,0]
NBs = S_iband[:,1:]

rlow = 555
rhigh = 695
