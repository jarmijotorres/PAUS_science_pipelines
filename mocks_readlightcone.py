# coding: utf-8
# %load mocks_readlightcone.py
#SRIPT REQUIRES NUMPY AND H5PY LIBRARIES
import numpy as np
import h5py

#LIGHTCONE FILENAMES
#lcCoreFilename = "/cosma5/data/durham/jarmijo/mocks_Durham/Gonzalez13.PAU.MillGas.field1.core.0.hdf5"
lcCoreFilename = '/cosma5/data/dp004/dc-armi2/LC_PAUS/Gonzalez13.PAU.MillGas.field1.core.0.hdf5'
lcPhotometryFilename = "/cosma5/data/dp004/dc-armi2/LC_PAUS/Gonzalez13.PAU.MillGas.field1.photometry.0.hdf5"
#lcEmlineFilename = "Gonzalez13.PAU.MillGas.field1.emlines.0.hdf5"
#lcKinematicsFilename = "Gonzalez13.PAU.MillGas.field1.kinematics.0.hdf5"

#LOAD LIGHTCONE HDF5 FILES
coreFile = h5py.File(lcCoreFilename, "r")
photometryFile = h5py.File(lcPhotometryFilename, "r")
#emlineFile = h5py.File(lcEmlineFilename, "r")
#kinematicsFile = h5py.File(lcKinematicsFilename, "r")

#READ IN PARAMETERS FROM HEADER
Omega0 =  np.array(coreFile.get("Header/" + "Omega0"))
print( "Omega0 = " + str(Omega0))

#READ IN DATASETS FROM APPROPRIATE LC FILE INTO NUMPY ARRAY
z_obs =  np.array(coreFile.get("Data/" + "z_obs"))
z_cos = np.array(coreFile.get("Data/" + "z_cos"))
mag = np.array(photometryFile.get("Data/" + "app10o_tot_ext"))
