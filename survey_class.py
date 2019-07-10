import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
#

survey_path = "/cosma/home/dp004/dc-armi2/PAU_test/OLD_cats/3835.csv"

survey_data = pd.read_csv(survey_path,sep=",",comment="#")
PAUS_data = pd.DataFrame(survey_data)

#print PAUS colums: 
    #print(PAUS_data.columns)
#

a = PAUS_data['ref_id'].values
OID = list(set(a))
NID = len(OID)
ra = PAUS_data['ra'].values
dec = PAUS_data['dec'].values
i_AB = PAUS_data['i_auto'].values
zspec = PAUS_data['zspec'].values
zphot = PAUS_data['zb_mean'].values
zpw =PAUS_data['pz_width'].values
nf = PAUS_data['band'].values
Nnf = PAUS_data['n_band'].values
fl = PAUS_data['flux'].values
fle = PAUS_data['flux_error'].values
zcosmos = PAUS_data['zp_gal'].values
iodds = PAUS_data['odds'].values
conf = PAUS_data['conf'].values
c1 = PAUS_data['chi2'].values
c2 = PAUS_data['zu99_gal'].values - PAUS_data['zl99_gal'].values


##### get spectra pipeline #####
band_count = Counter(nf)
df = pd.DataFrame.from_dict(band_count, orient='index')
bands = df.axes[0]
ID_sorted = np.sort(a)
nf_sorted = nf[np.argsort(a)]
flux_sorted = fl[np.argsort(a)]
bandss = np.sort(bands)
#
flx_obj = []
nf_obj =[]
err_obj = []
for i in range(NID):
    nf_obj.append(nf[np.where(a == OID[i])[0]])
    flx_obj.append(fl[np.where(a == OID[i])[0]])
    err_obj.append(fle[np.where(a == OID[i])[0]])
#
flux_obj_srt_nf = []
band_obj_srt_nf = []
fluxerr_obj_srt_nf = []
for i in range(NID):
    flux_obj_srt_nf.append(flx_obj[i][np.argsort(nf_obj[i])])
    band_obj_srt_nf.append(np.sort(nf_obj[i]))
    fluxerr_obj_srt_nf.append(err_obj[i][np.argsort(nf_obj[i])])
    
file_spectra = open('/cosma/home/dp004/dc-armi2/PAU_test/catalogs/3835_ID_narrow_fluxes.dat','a')
file_bands = open('/cosma/home/dp004/dc-armi2/PAU_test/catalogs/3835_ID_narrow_bands.dat','a')
file_errors = open('/cosma/home/dp004/dc-armi2/PAU_test/catalogs/3835_ID_narrow_fluxerror.dat','a')

for i in range(NID):
    for l in range(len(flux_obj_srt_nf[i])):
        file_spectra.write('%.5lf '% flux_obj_srt_nf[i][l])
        file_bands.write('%s '%band_obj_srt_nf[i][l])
        file_errors.write('%s '%fluxerr_obj_srt_nf[i][l])
    file_spectra.write('\n')
    
file_spectra.close()
file_bands.close()
file_errors.close()

############### One-object catalog ################
ID_PAU = []
ra_PAU = []
dec_PAU = []
iauto_PAU = []
zs = []
zp = []
zp_width = []
zc = []
CONF = []
ODDS = []
Qz = []
Nb = []
cu = a[0]
cd = a[-1]
C1 = []
C2 = []
for i in range(len(a)):
    if a[i] != cu:
        #Qz.append(qz[i-1])
        ODDS.append(iodds[i-1])
        CONF.append(conf[i-1])
        ID_PAU.append(a[i-1])
        iauto_PAU.append(i_AB[i-1])
        zs.append(zspec[i-1])
        zp.append(zphot[i-1])
        zp_width.append(zpw[i-1])
        zc.append(zcosmos[i-1])
        ra_PAU.append(ra[i-1])
        dec_PAU.append(dec[i-1])
        Nb.append(Nnf[i-1])
        C1.append(c1[i-1])
        C2.append(c2[i-1])
    cu = a[i]
    if a[i] == cd:
        #Qz.append(qz[i])
        ODDS.append(iodds[i])
        CONF.append(conf[i])
        ID_PAU.append(a[i])
        iauto_PAU.append(i_AB[i])
        zs.append(zspec[i])
        zp.append(zphot[i])
        zp_width.append(zpw[i-1])
        zc.append(zcosmos[i])
        ra_PAU.append(ra[i])
        dec_PAU.append(dec[i])
        Nb.append(Nnf[i])
        C1.append(c1[i])
        C2.append(c2[i])
        break
        
S = np.array([ID_PAU,ra_PAU,dec_PAU,zs,zp,zp_width,zc,iauto_PAU,CONF,ODDS,C1,C2,Nb])
head = 'ID_PAU\t ra\t dec\t zspec\t zphot\t zp_width\t zcosmos\t imag\t confidence_flag\t Odds\t chi2\t zu99 - zl99\t N_bands\n'
fm = '%i %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.1lf %.5lf %.5lf %.5lf %.5lf %i'
np.savetxt('/cosma/home/dp004/dc-armi2/PAU_test/catalogs/cat_3835.dat',S.T,header= head,fmt=fm)