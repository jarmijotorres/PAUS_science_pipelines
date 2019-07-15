import numpy as np
import h5py,sys
from astropy.cosmology import Planck15 as cosmo
from astropy.table import Table

icut = sys.argv[1]
lcCoreFilename = "/cosma5/data/durham/jarmijo/mocks_Durham/Gonzalez13.PAU.MillGas.field1.core.0.hdf5"
lcPhotometryFilename = "/cosma5/data/durham/jarmijo/mocks_Durham/Gonzalez13.PAU.MillGas.field1.photometry.0.hdf5"
#lcEmlinesFilename = "/cosma5/data/durham/jarmijo/mocks_Durham/Gonzalez13.PAU.MillGas.field1.emlines.0.hdf5"
coreFile = h5py.File(lcCoreFilename, "r")
photometryFile = h5py.File(lcPhotometryFilename, "r")
#emlinesFile = h5py.File(lcEmlinesFilename, "r")
#
z_obs =  np.array(coreFile.get("Data/" + "z_obs"))
ra = np.array(coreFile.get("Data/" + "ra"))
dec = np.array(coreFile.get("Data/" + "dec"))
#LHalpha = np.array(emlinesFile.get("Data/" + "L_tot_Halpha_ext"))
mstardot = np.array(coreFile.get("Data/" + 'mstardot'))
SDSS_i = np.array(photometryFile.get("Data/" + "appiSo_tot_ext"))
SDSS_r = np.array(photometryFile.get("Data/" + "apprSo_tot_ext"))
SDSS_g = np.array(photometryFile.get("Data/" + "appgSo_tot_ext"))
SDSS_u = np.array(photometryFile.get("Data/" + "appuSo_tot_ext"))
#
SDSS_I = np.array(photometryFile.get("Data/" + "magiSr_tot_ext"))
SDSS_R = np.array(photometryFile.get("Data/" + "magrSr_tot_ext"))
SDSS_G = np.array(photometryFile.get("Data/" + "maggSr_tot_ext"))
SDSS_U = np.array(photometryFile.get("Data/" + "maguSr_tot_ext"))
#
table = np.array([ra,dec,z_obs,SDSS_g,SDSS_i,SDSS_r,SDSS_u,SDSS_G,SDSS_I,SDSS_R,SDSS_U,mstardot])
NBs = np.arange(0,len(table.T)).reshape(1,len(table.T)) # Narrow band magnitudes (observed frame)
nb = 40
print ("reading mocks magnitudes... \n")
for i in range(nb):
    s = str(10 + i)
    mag = np.array(photometryFile.get("Data/" + "app"+s+"o_tot_ext"))
    NBs = np.vstack([NBs,mag])
#
############## selection ##################
table = table.T
NBs = NBs.T
table = table[SDSS_i<float(icut)]
NBs = NBs[SDSS_i<float(icut)]
##########################################
print ("done.\n")
lamb = np.zeros(nb)
for l in range(nb):
    ni = 100*l
    lamb[l] = 4550 + ni
xi = 4050
xf = 4450
P_B = []
print ("summing fluxes...\n")
for i in range(len(table)):
    xs = lamb/(1+table[i,2])
    dx0 = np.diff(lamb)[0]
    dx = np.diff(xs)[0]
    dx2 = 0.5*dx
    mags = NBs[i,1:]
    Ix = xs[(xs>xi - dx2)&(xs<xf + dx2)]
    IF = mags[(xs>xi - dx2)&(xs<xf + dx2)]
    F_in = 10**(((IF - cosmo.h)+48)/-2.5)
#    F_in /= (1+table[i,2])
    for c in range(len(Ix)):
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
    P_B.append(np.sum(F_in*dx)/dx0)
PB = np.array(P_B)
PAUS_BLUE = PB[(table[:,2]<0.9)&(table[:,2]>0.11)]
new_table = table[(table[:,2]<0.9)&(table[:,2]>0.11)]
mB = -2.5*np.log10(PAUS_BLUE) - 48 # from flux to PAUS magnitud                                                                                                                                                                                                                                                                                                                                                                                                 e
dL = cosmo.luminosity_distance(new_table[:,2]).value * cosmo.h
MB = mB - 25 - 5*np.log10(dL)
new_table = np.vstack([new_table.T,MB]).T
t = Table(new_table,names=['RA','DEC','Z','SDSS_g','SDSS_i','SDSS_r','SDSS_u','SDSS_G','SDSS_I','SDSS_R','SDSS_U','SFR Mdot','MB'],)
print "saving data from "+str(len(MB))+" of " +str(i) + "galaxies in fits format... at /cosma/home/durham/jarmijo/PAU_test/catalogs/mocks_MBlue_SDSSphot.fits\n"
t.write('/cosma/home/dp004/dc-armi2/PAU_test/catalogs/mocks_radecz_MIMB_SFRHaplha_'+icut+'cut.fits',format='fits')
print "end of program.\n"

