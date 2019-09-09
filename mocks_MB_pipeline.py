# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:14:00 2019

@author: jarmijo
"""

import numpy as np
import h5py,sys
from astropy.cosmology import Planck15 as cosmo
from astropy.table import Table

icut = sys.argv[1]
localdir = '/home/jarmijo/Documents/mocks/'
lcCoreFilename = localdir+"Gonzalez13.PAU.MillGas.field1.core.0.hdf5"
lcPhotometryFilename = localdir+"Gonzalez13.PAU.MillGas.field1.photometry.0.hdf5"
#lcEmlinesFilename = localdir+"Gonzalez13.PAU.MillGas.field1.emlines.0.hdf5"
coreFile = h5py.File(lcCoreFilename, "r")
photometryFile = h5py.File(lcPhotometryFilename, "r")
#emlinesFile = h5py.File(lcEmlinesFilename, "r")

z_obs =  np.array(coreFile.get("Data/" + "z_obs"))
ra = np.array(coreFile.get("Data/" + "ra"))
dec = np.array(coreFile.get("Data/" + "dec"))
#LHalpha = np.array(emlinesFile.get("Data/" + "L_tot_Halpha_ext"))
mstardot = np.array(coreFile.get("Data/" + 'mstardot'))
magI = np.array(coreFile.get("Data/" + "appMio_tot_ext"))#
MagI = np.array(coreFile.get("Data/" + "magMir_tot_ext"))
table = np.array([ra,dec,z_obs,magI,MagI,mstardot])
nb = 40
print ("reading mocks magnitudes... \n")
for i in range(nb):
    s = str(10 + i)
    mag = np.array(photometryFile.get("Data/" + "app"+s+"o_tot_ext"))
    table = np.vstack([table,mag])
table = table.T
#
############## selection ##################
table = table[magI<float(icut)]
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
    mags = table[i,7:]
    Ix = xs[(xs>xi - dx2)&(xs<xf + dx2)]
    IF = mags[(xs>xi - dx2)&(xs<xf + dx2)]
    F_in = 10**((IF+48)/-2.5)
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
new_table = table[:,(0,1,2,3,4,5,6)][(table[:,2]<0.9)&(table[:,2]>0.11)]
mB = -2.5*np.log10(PAUS_BLUE) - 48
MB = mB - 25 - 5*np.log10(cosmo.luminosity_distance(new_table[:,2]).value)
new_table = np.vstack([new_table.T,MB]).T
t = Table(new_table,names=['RA','DEC','Z','m_i','M_I','FHalpha','SFR','MB'])
print "saving data from "+str(len(MB))+" of " +str(i) + "galaxies in fits format... at /cosma/home/durham/jarmijo/PAU_test/catalogs/mocks_MBlue.fits\n"
t.write('/cosma/home/dp004/dc-armi2/PAU_test/catalogs/mocks_radecz_MIMB_SFRHaplha_'+icut+'cut.fits',format='fits')
print "end of program.\n"
