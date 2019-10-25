import numpy as np
from mpi4py import MPI
import h5py,time
from glob import glob
import matplotlib.pyplot as plt

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
        for nb in range(35,50):
            NB_i = np.array(sp_ivol.get("/Output001/mag%dr_tot_ext"%nb))
            Mi_res = np.vstack([Mi_res,NB_i])
        if sn == 0: 
            S_iband = Mi_res.T
        else:
            S_iband = np.append(S_iband,Mi_res.T,axis=0)

Mi = S_iband[:,0]
NBs = S_iband[:,1:]

##### i band integration over galaxies at z = 0 snapshot.
PAUS_ifrd = np.sum(NBs,axis=1) / 15.

N1 = np.histogram(Mi,bins=30,range=(-25.,-4.5))
N2 = np.histogram(PAUS_ifrd,bins=30,range=(-25.,-4.5))
##### 
bb = N1[1][:-1] + np.diff(N1[1])[0]/2.
ts = np.array([bb,N1[1][:-1],N1[1][1:],N1[0],N2[0]]).T
print('file saved in /cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/N_Mi_PAUSi_z0.0_all.txt\n' )
np.savetxt('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/N_Mi_PAUSi_z0.0_all.txt',ts,header='center_bin ei ef N_Mi N_PAUSi',fmt='%.6lf %.6lf %.6lf %d %d')
#f,ax = plt.subplots(1,1,figsize=(7,6))
#ax.hist(Mi,bins=30,range=(Mi.min(),Mi.max()),histtype='step',linewidth=2.,color='lightcoral',label='i-band MegaCam')
#ax.hist(PAUS_ifrd,bins=30,range=(Mi.min(),Mi.max()),histtype='step',linewidth=2.,color='maroon',label='PAUS i (15 NBs)')
#ax.legend(prop={'size':14})
#ax.set_xlabel('Restframe magnitude')
#ax.set_ylabel('Counts')
#ax.set_yscale('log')
#plt.tight_layout()
#plt.savefig('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/figures/October/Mi_PAUSi_z0.0.png',bbox_inches='tight')

##### read response function for NBs



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
