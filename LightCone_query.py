import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import os
#import sys
from astropy.table import Table
from optparse import OptionParser
#import time


import astropy.io.fits as fits

import eagleSqlTools as sql # do "pip3 install eagleSqlTools --user"


# =============================================================================
# #===============SQL CONNECTION================#
# con = sql.connect('nxr908',password='hmi20hwk')
# #=============================================#
# 
# SEARCH_RADIUS = 1.0 # deg
# NSIDE = 1 # 12*NSIDE**2 points
# 
# parser = OptionParser()
# parser.add_option('--sql',dest='sql_force',action='store_true',default=False)
# 
# (options,args) = parser.parse_args()
# 
# ra_vals = [45.0]
# dec_vals = [45.0]
# 
# DIR = './EuclidGonzalez14SQL/'
# if os.path.exists(DIR) is False:
#     os.mkdir(DIR)
# 
# 
# zbins = np.arange(0,0.3,0.005) 
# fig,ax = plt.subplots()
# 
# for ii,(ra,dec) in enumerate(zip(ra_vals,dec_vals)):
#     cat_path = DIR + './Euclid_LC_DEEP_z_03_i23.cat'
#     
#     if True: #(os.path.exists(cat_path) is False) | (options.sql_force): # ignore this if statement...
# 
#         
# 
#         print('\033[33;1mStart query %i of %i...\033[0m' %(ii+1,1))
# 
#         opt_mag_cols = ','.join('DES_%s_obs_app as %s_kronMag'%(b,b.lower()) for b in 'g r i z Y'.split())
#         nir_mag_cols = ','.join('mag_%s_obs_app as %s_kronMag'%(b,b[0]) for b in 'J H K'.split())
# 
# 
#         query = (
#             ''' select GalaxyID,ra as NIR_ra,dec as NIR_dec,%s,%s,''' %(opt_mag_cols,nir_mag_cols)
#             +''' redshift_cos,redshift_obs''' 
#             +''' from EUCLID_v1..LC_DEEP_Gonzalez2014a ''' # LC_DEEP_Gonzalez2014a
#             +''' where (redshift_cos < 0.3) '''
#             +''' and ((ra-(%.3f))*(ra-(%.3f)) + (dec-(%.3f))*(dec-(%.3f))<%.2f)''' %(ra,ra,dec,dec,SEARCH_RADIUS)
#             +''' and ( DES_i_obs_app < 23.0 ) '''
#             )
# 
#         print(query)
# 
#         data = sql.execute_query(con,query)
# 
#         print("\033[31;1m%i rows.\033[0m"%len(data))
#         
#         df = pd.DataFrame(data)
#         #df.to_csv(csv_path,index=False)
# 
#         cat = Table.from_pandas(df)
# 
# 
#         #J_stell = Column(np.zeros(len(cat)),name='J_stell')
#         #K_stell = Column(np.zeros(len(cat)),name='K_stell')
# 
#         #cat.add_columns([J_stell,K_stell])
# 
#         cat.write(cat_path,format='fits',overwrite=True)
# 
# 
#     else:
#         print('cat exists!')
#         cat = fits.open(cat_path)[1].data
# 
#     print(np.histogram(cat['redshift_cos'],bins=zbins))
# 
#     ax.hist(cat['redshift_cos'],bins=zbins,histtype='step')
# 
# plt.show()
# 
# =============================================================================
con = sql.connect('nxr908',password='hmi20hwk')

ra = 180.0
dec = 40.0
a = 75
b = 45
#
SEARCH_RADIUS = 1.0 # deg
#opt_mag_cols = 'dasd'
#nir_mag_cols = 'asds'
query = (
    ''' select ra, dec, z_obs, app_mag, abs_mag'''
    +''' from Smith2017a..Galaxies ''' # LC_DEEP_Gonzalez2014a
    +''' where (z_cos < 0.5) and app_mag < 23'''
    +''' and ((ra-(%.3f))*(ra-(%.3f))/%.2lf/%.2lf + (dec-(%.3f))*(dec-(%.3f))/%.2lf/%.2lf<%.2f)''' %(ra,ra,a,a,dec,dec,b,b,SEARCH_RADIUS)
    )

data = sql.execute_query(con,query)
df = pd.DataFrame(data)
         #df.to_csv(csv_path,index=False)
# 
cat = Table.from_pandas(df)
cat_path = '../test.fits'
cat.write(cat_path,format='fits',overwrite=True)
cat.write('/home/jarmijo/Dropbox/test_LC.fits',format='fits',overwrite=True)
print("cata created saved in /home/jarmijo/Dropbox/test_LC.fits...\n")
