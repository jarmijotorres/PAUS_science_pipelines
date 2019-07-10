import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
from astroML.stats import binned_statistic_2d
import matplotlib as mpl

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='inout')
plt.rc('ytick',direction='inout')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='sans-serif')
plt.rc('font',size=18)

f1 = fits.open('/Users/jarmijo/Documents/Mocks/mocks_radecz_MIMB_SFRHaplha_23cut_fix_nofrf.fits')
data1 = f1[1].data
dL = cosmo.luminosity_distance(data1['Z'])

Ns = 9
zmin = 0.11
zmax = 0.9
M_min = -16.
M_max = -23.5
b = 30 #N of bins LF
dM = abs(M_max-M_min)/b
zbins = np.linspace(0.11,0.9,Ns,endpoint=True)

# =============================================================================
# f,ax = plt.subplots(1,1,figsize=(7,6))
# ax.hist(data1['M_I'],bins=b,color='r',histtype='step',label='Mock')
# ax.set_yscale('log')
# ax.vlines(M_min,1,len(data1['M_I']),linestyles='--',colors='k')
# ax.vlines(M_max,1,len(data1['M_I']),linestyles='--',colors='k')
# ax.set_ylim(1e1,5*len(data1['M_I'])/(b))
# ax.set_xlabel(r'$M_i - 5\log h$')
# ax.set_ylabel('Counts')
# ax.legend()
# plt.show()
# 
# =============================================================================
#Mbins limits given the M_i luminosity function (see MI_LF.py)
M_faint_cut = -1*np.array([16.,16.5,17.25,18.,19.,19.5,19.75,20.5])
M_bright_cut = -23.
# =============================================================================

I = data1['M_I']
B = data1['MB']
Z = data1['Z']
Bs,Is = [],[]
nf = 2
nc = 2
f,ax = plt.subplots(nf,nc,figsize=(4*nf,4.),sharex=True,
                        sharey=True,gridspec_kw={'width_ratios': [1]*nc, 'height_ratios': [1]*nf})
for s in range(Ns-1):
    zi = zbins[s]
    zf = zbins[s+1]
    xx = I[(Z>zi)&(Z<zf)]
    yy = B[(Z>zi)&(Z<zf)]
    cc = Z[(Z>zi)&(Z<zf)]
    N, xedges, yedges = binned_statistic_2d(xx,yy,xx,'count', bins=50,range=[(-26,-15),(-26,-15)])
    xb = xedges[:-1]+np.diff(xedges)[0]/2.;yb = yedges[:-1]+np.diff(yedges)[0]/2.
    if s < (Ns - 1)/2:
        i = 0
        j = s
        ax[i,j].contourf(xb,yb,np.log10(N.T),cmap='PuBu_r',origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1],],levels=[1,2,3,4,5],zorder=10)
        ax[i].scatter(xx,(yy),c='b',s=1.,label='%.2lf < z < %.2lf'%(zi,zf),zorder=0)
        ax[i].legend(loc=4,prop={'size':12})
        ax[i].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
    elif s >= (Ns - 1)/2:
        i = 1
        j = s - i*(Ns - 1)/2
        ax[i,j].contourf(xb,yb,np.log10(N.T),cmap='PuBu_r',origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1],],levels=[1,2,3,4,5],zorder=10)
        ax[i].scatter(xx,(yy),c='b',s=1.,label='%.2lf < z < %.2lf'%(zi,zf),zorder=0)
        ax[i].legend(loc=4,prop={'size':12})
        ax[i].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
ax[0,0].set_xticks([-16,-18,-20,-22,-24])
ax1 = f.add_axes([0.99, 0.23, 0.01, 0.72])
cmap = mpl.cm.PuBu_r
norm = mpl.colors.Normalize(vmin=1., vmax=5.)

cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical',ticks=[1,2,3,4,5])
cb1.set_label(r'$\log\ N$')
ax[0].set_ylabel(r'$M_{Blue}$')
ax[1].set_xlabel(r'$M_{i}$')
ax[1].xaxis.set_label_coords(1.1,-0.15)
plt.tight_layout()
plt.subplots_adjust(wspace=0.00)
plt.show()
    