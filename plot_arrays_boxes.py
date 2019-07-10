# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch 
bl = 0.
bu = 1.
Nb = 11
BINS = np.linspace(bl,bu,Nb)#
b = BINS[:-1] + np.diff(BINS)[0]/2.
z0 = 0.05
sz = 0.1
for i in range(1):
     zi = z0 + sz*i
     zf = zi + sz
     xx = S_int(Qz_f_no[(zp_f_no>zi)&(zp_f_no<zf)])
     yy = dz_no[(zp_f_no>zi)&(zp_f_no<zf)]
     bmedian = np.zeros_like(b)
     p84 = np.zeros_like(b)
     p16 = np.zeros_like(b)
     p5 = np.zeros_like(b)
     p95 = np.zeros_like(b)
     N = np.zeros_like(b)
     Ni = np.zeros_like(b)
     for i in range(Nb-1):
         bmedian[i] = np.median(yy[(xx>BINS[i])&(xx<BINS[i+1])])
         p84[i] = np.percentile(yy[(xx>BINS[i])&(xx<BINS[i+1])],84)
         p16[i] = np.percentile(yy[(xx>BINS[i])&(xx<BINS[i+1])],16)
         p5[i] = np.percentile(yy[(xx>BINS[i])&(xx<BINS[i+1])],5)
         p95[i] = np.percentile(yy[(xx>BINS[i])&(xx<BINS[i+1])],95)
         N[i] = len(yy[(xx>BINS[i])&(xx<BINS[i+1])])
         Ni[i] = len(yy[(xx>BINS[i])&(xx<BINS[i+1])&(abs(yy)<0.05)])
     N_out = (N-Ni)/N
     f,ax = plt.subplots(2,1,figsize=(7,8),sharex=True,gridspec_kw={'height_ratios': [1,3]})
     ax[1].text(0.1,0.04,'%.2lf < z < %.2lf'%(zi,zf))
     ax[1].plot(xx,(yy),'k.',ms=1.)
     md, = ax[1].plot(b,bmedian,'m*',ms=12.,label='Median data')
     first_legend = ax[1].legend(handles=[md], loc=3,prop={'size':15})
     ax[1].add_artist(first_legend)
     ax[1].hlines(0,-0.01,1.3,color='b',linestyles='--',alpha=0.7)
     for i in range(len(b)):
         ax[1].fill_between([BINS[i],BINS[i+1]],p5[i],p95[i],facecolor='r',alpha=0.4,edgecolor='r')
         ax[1].fill_between([BINS[i],BINS[i+1]],p16[i],p84[i],facecolor='g',alpha=0.4,edgecolor='g')
     cr = [Patch(facecolor='green',edgecolor='green',label='68% confidence',alpha=0.4),
       Patch(facecolor='red',edgecolor='red',label='95% confidence',alpha=0.4)]
     second_legend = ax[1].legend(handles=[cr[0],cr[1]], loc=(0.02,0.1),prop={'size':15})
     ax[1].add_artist(second_legend)
     ax[0].plot(b,N_out,'k+-',ms=12.,label='Fraction of outliers')
     ax[0].legend(loc=2,prop={'size':15})
     ax[0].vlines(BINS,-0.01,1.01,linestyle='--',alpha=0.4)
     ax[1].set_ylabel(r'$(z_{ph} - z_s)/(1+z_s)$')
     ax[1].set_xlabel(r'n(<Qz)')
     ax[1].set_ylim(-.05,.05)
     ax[0].set_ylim(0.,0.3)
     ax[0].set_yticks([0.0,0.1,0.2,0.3])
     ax[1].set_xlim(-0.01,1.01)
     plt.tight_layout()
     plt.subplots_adjust(hspace = 0.05)
     #plt.savefig('figures_png/dz_nQz_z%.2lf_%.2lf.png'%(zi,zf))
     plt.close()
