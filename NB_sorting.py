import numpy as np
from glob import glob
#
with open('/cosma/home/dp004/dc-armi2/PAU_test/mask/filters_unique.dat','r') as f:
    responses = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each 
#
PAUS_filters = []
for i in range(len(responses)):
    PAUS_filters.append(np.array(responses[i].split('PAU-')))
ini = 63522
deltaNB = 212
PAUS_filters_ascii = []
PAUS_filters_ascii.append(responses[63312:63522])
for i in range(40):
    PAUS_filters_ascii.append(responses[ini:ini+deltaNB])
    ini += deltaNB
#
### save each band in ascii file
l_cen = 4550
for nb in range(1,41):
    file_filter = open('/cosma/home/dp004/dc-armi2/PAU_test/filters/PAUS_'+str(l_cen)+'_response.dat','a')
    for i in range(len(PAUS_filters_ascii[nb])):
        file_filter.write(PAUS_filters_ascii[nb][i])
    file_filter.close()
    l_cen += 100
#
NB_names = glob('/cosma/home/dp004/dc-armi2/PAU_test/filters/PAUS_*')
PAUS_NBs = []
for i in range(len(NB_names)):
    PAUS_NBs.append(np.loadtxt(NB_names[i]))
c_nm = np.zeros(40)
for i in range(40):
    c_nm[i-1] = np.trunc(np.round(np.average(PAUS_NBs[i][:,0],weights=PAUS_NBs[i][:,1])/10.)/10.)
#### sort and save bands by central wavelength 
PAUS_NBs_sorted = np.array(PAUS_NBs)[1:][np.argsort(c_nm)]
#
lini=4550
for i in range(40):
    np.savetxt('/cosma/home/dp004/dc-armi2/PAU_test/filters/PAUS_'+str(lini)+'_band.dat',PAUS_NBs_sorted[i],fmt='%.6lf',header='wavelength response')
    lini+=100
#
##### plot to test
#NB_names = glob('/cosma/home/dp004/dc-armi2/PAU_test/filters/PAUS_*')
#PAUS_NBs = []
#for i in range(len(NB_names)):
#    PAUS_NBs.append(np.loadtxt(NB_names[i]))
#f,ax = plt.subplots(1,1,figsize=(7,6))
#for i in range(40):
#    ax.plot(PAUS_NBs[i][:,0],PAUS_NBs[i][:,1],'-',c=np.random.random(3))
#ax.set_xlim(4450,8550)
#plt.show()