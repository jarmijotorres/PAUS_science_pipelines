import numpy as np
from mpi4py import MPI
import h5py,time
from glob import glob
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='inout')
plt.rc('ytick',direction='inout')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='sans-serif')
plt.rc('font',size=16)

### read snapshot files (sorted in 512 subvolumes, each one with 34 redshift outputs)
spName = '/data/dega1/ddmn39/galform_out/r504/lightcone_input/Gonzalez13.PAU.MillGas/ivol_0/galaxies.hdf5'

### read lightcone fragment (several parts)
lcChunkName = '/data/dega1/ddmn39/lightcone_out/LC148/GAL504/PAU.MillGas/Gonzalez13.PAU.MillGas.field1.500.hdf5'

### read full lightcone
lcCoreName = '/cosma5/data/dp004/dc-armi2/LC_PAUS/Gonzalez13.PAU.MillGas.field1.core.0.hdf5'

LC_chunk = h5py.File(lcChunkName,'r')
LC_Core = h5py.File(lcCoreName,'r')
sp_ivol = h5py.File(spName,'r')

#sp_O1 = np.array(sp_ivol.get('/Output001'))
#O1_x = np.array(sp_ivol.get("/Output001/xgal"))
#O1_y = np.array(sp_ivol.get("/Output001/ygal"))
#O1_z = np.array(sp_ivol.get("/Output001/zgal"))
## open Output 5 to 8 for 0.1 < z < 0.2 bin

subvols_dir = np.sort(glob('/data/dega1/ddmn39/galform_out/r504/lightcone_input/Gonzalez13.PAU.MillGas/ivol_*'))
#CFHTLS_i_obs = None #sp_O1['magMio_tot_ext'] observed frame
#CFHTLS_i_res = None #sp_01magMir_tot_ext] rest frame
S1 = np.zeros(0,dtype=object)
S2 = np.zeros(0,dtype=object)
S3 = np.zeros(0,dtype=object)
S4 = np.zeros(0,dtype=object)
Ss = [S1,S2,S3,S4]
print('reading snapshots...\n')
for sn,sv in enumerate(subvols_dir):
    spName = sv + '/galaxies.hdf5'
    sp_ivol = h5py.File(spName,'r')
    for i in range(4):
        O_x = np.array(sp_ivol.get("/Output00"+str(i+5)+"/xgal"))
        O_y = np.array(sp_ivol.get("/Output00"+str(i+5)+"/ygal"))
        O_z = np.array(sp_ivol.get("/Output00"+str(i+5)+"/zgal"))
        Mi_obs = np.array(sp_ivol.get("/Output00"+str(i+5)+"/magMio_tot_ext"))
        if Mi_obs.any() == None: break
        pos_Mi = np.array([O_x,O_y,O_z,Mi_obs]).T
        if sn == 0: 
            Ss[i] = pos_Mi
        else:
            Ss[i] = np.append(Ss[i],pos_Mi,axis=0)
#
boxsize=500.
Vbox = boxsize**3.

Nd = 5 #per side, so the number of subdiv is Nd**3.
Nsubd = Nd**3

sboxsize = boxsize / Nd
Vsubbox = sboxsize**3.

#take a snapshot compute the luminosity function for Nd**3 subvols
print('creating IDs subvols...\n')
IDs = []
for i in range(len(Ss)):
    N = len(Ss[i])
    id_cell = np.zeros(N,dtype=int)
    for p in range(N):
        lp = (Ss[i][p,:3] / sboxsize).astype(int)
        id_cell[p] = lp[0]*Nd**2 + lp[1]*Nd + lp[2]
    IDs.append(id_cell)
# ==== luminosity function parameters === #
b = 32 #N of bins LF
M_min = -16.
M_max = -24.
Mbins= np.linspace(M_max,M_min,b+1,endpoint=True)
dM = abs(M_max-M_min)/(b) # Nb = 32 between (-16,-24) --> dM = 0.25 (Nb = 40 --> dM = 0.2 for the same range)
#
print('computing LF snapshots...\n')
L_LF = []
for i in range(len(Ss)):
    Ox = Ss[i]
    LFs = np.zeros((125,32),dtype = float)
    for c in range(Nsubd):
        s1 = Ox[IDs[i] == c]
        LFs[c] = np.histogram(s1,bins=Mbins)[0]
    L_LF.append(LFs)
print('saving in: /cosma5/data/dp004/dc-armi2/PAU/PAU_test/data/sLF_125subvols_Vsbox100Mpc_h_s57.npy')
for l in range(len(Ss)):
    np.save('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/data/sLF_125subvols_Vsbox100Mpc_h_s_'+str(57-l)+'.npy',L_LF[l])
    
print('\n end of program.\n')

################ output checking #################
import numpy as np
from mpi4py import MPI
import h5py,time
from glob import glob
import matplotlib.pyplot as plt
# load info
sLF_files = glob('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/data/sLF_125subvols_Vsbox100Mpc_h_s_*')
#details:
    # s57: z = 0.089288
    # s56: z = 0.115883
    # s55: z = 0.144383
    # s54: z = 0.174898
z57 = 0.089288
z56 = 0.115883
z55 = 0.144383
z54 = 0.174898
#
sLFs = []
sLF_mean = []
for i in range(len(sLF_files)):
    LFs = np.load(sLF_files[i])
    sLFs.append(LFs)
    sLF_mean.append(np.mean(LFs,axis=0))
#
b = 32 #N of bins LF
M_min = -16.
M_max = -24.
Mbins= np.linspace(M_max,M_min,b+1,endpoint=True)
dM = abs(M_max-M_min)/(b) # Nb = 32 between (-16,-24) --> dM = 0.25 (Nb = 40 --> dM = 0.2 for the same range)
bb = Mbins[:-1] + np.diff(Mbins)[0]/2.
#

#
