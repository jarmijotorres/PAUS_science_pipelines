import numpy as np
from mpi4py import MPI
import h5py
from glob import glob
from astropy.table import Table

subvols_dir = np.sort(glob('/data/dega1/ddmn39/galform_out/r504/lightcone_input/Gonzalez13.PAU.MillGas/ivol_*'))

S_iband = np.zeros(0,dtype=object)
#S_rband = np.zeros(0,dtype=object)

print('reading snapshots... \n')
for sn,sv in enumerate(subvols_dir):
    
    spName = sv + '/galaxies.hdf5'
    sp_ivol = h5py.File(spName,'r')
    #O_x = np.array(sp_ivol.get("/Output00"+str(i+8)+"/xgal"))
    #O_y = np.array(sp_ivol.get("/Output00"+str(i+8)+"/ygal"))
    #O_z = np.array(sp_ivol.get("/Output00"+str(i+8)+"/zgal"))
    ##### test i-band at z=0 snapshot
    Mi_res = np.array(sp_ivol.get("/Output007/magrSr_tot_ext"))
    #Mr_res = np.array(sp_ivol.get("/Output001/magrSr_tot_ext"))
    if Mi_res.any() == None: 
        print('M_i not available on subvol %d\n'%sn)
        #break
    elif Mi_res.any() != None: 
    #
        for nb in range(10,50): #i-band range defined by the broad i-band
            NB_i = np.array(sp_ivol.get("/Output007/mag%dr_tot_ext"%nb))
            Mi_res = np.vstack([Mi_res,NB_i])
        if sn == 0: 
            S_iband = Mi_res.T
        else:
            S_iband = np.append(S_iband,Mi_res.T,axis=0)
print('done. \n Creating taable...\n')

head = []
head.append('Mi_broad')
for i in range(10,50):
    head.append('NB_%d'%i)
Ms_table = Table(data=S_iband,names=head,)

print('done.\n Table written on /cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/MrbroadSDSS_NBs_10_50_z_0.144_s007.hdf5 \n')
Ms_table.write('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/outputs_pipelines/MrbroadSDSS_NBs_10_50_z_0.144_s007.hdf5',format='hdf5',path='data',compression=True)
print('\end of program.\n')