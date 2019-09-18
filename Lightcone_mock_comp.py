import numpy as np
import h5py

spName = '/data/dega1/ddmn39/galform_out/r504/lightcone_input/Gonzalez13.PAU.MillGas/ivol_0/galaxies.hdf5'

lcChunkName = '/data/dega1/ddmn39/lightcone_out/LC148/GAL504/PAU.MillGas/Gonzalez13.PAU.MillGas.field1.500.hdf5'
lcCoreName = '/cosma5/data/dp004/dc-armi2/LC_PAUS/Gonzalez13.PAU.MillGas.field1.core.0.hdf5'

LC_chunk = h5py.File(lcChunkName,'r')
LC_Core = h5py.File(lcCoreName,'r')
sp_ivol = h5py.File(spName,'r')

sp_O1 = np.array(sp_ivol.get('/Output001'))
O1_x = np.array(sp_ivol.get("/Output001/xgal"))
O1_y = np.array(sp_ivol.get("/Output001/ygal"))
O1_z = np.array(sp_ivol.get("/Output001/zgal"))
## open Output 5 to 8 for 0.1 < z < 0.2 bin

subvols_dir = glob('/data/dega1/ddmn39/galform_out/r504/lightcone_input/Gonzalez13.PAU.MillGas/ivol_*')

sp_Pos = np.zeros(4,dtype=object)
for sv in subvols_dir:
    spName = sv + '/galaxies.hdf5'
    sp_ivol = h5py.File(spName,'r')
    for i in range(5,9):
        O_x = np.array(sp_ivol.get("/Output00"+str(i)+"/xgal"))
        O_y = np.array(sp_ivol.get("/Output00"+str(i)+"/ygal"))
        O_z = np.array(sp_ivol.get("/Output00"+str(i)+"/zgal"))
