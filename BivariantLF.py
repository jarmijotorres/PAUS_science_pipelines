
# coding: utf-8

# # PAUS Blue, I-auto K-corrected

# In[320]:


from astroML.stats import binned_statistic_2d
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='inout')
plt.rc('ytick',direction='inout')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='sans-serif')
plt.rc('font',size=18)


# In[117]:


D = np.loadtxt('../catalogs/Imag_MB.dat')
Imag = D[:,0]; MB = D[:,1]; zB = D[:,2]


# In[118]:


z_mean, xe, ye = binned_statistic_2d(Imag, MB, zB,'mean', bins=50)
N, xedges, yedges = binned_statistic_2d(Imag, MB, zB,'count', bins=50)


# In[48]:


fig = plt.figure(figsize=(15, 8))
fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95,
                    bottom=0.07, top=0.95)

#--------------------
# First axes:
plt.subplot(121)
plt.imshow(np.log10(N.T), origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', interpolation='nearest', cmap='jet')
plt.xlim(-27.5,-15)
plt.ylim(-32.5,-15)
plt.xlabel('M$_i - 5\log h $')
plt.ylabel(r'M$_{Blue} $')

cb = plt.colorbar(ticks=[0, 1, 2, 3],
                  format=r'$10^{%i}$', orientation='horizontal')
cb.set_label('number in pixel')
#plt.clim(0, 3)

#--------------------
# Third axes:
plt.subplot(122)
plt.imshow(z_mean.T, origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', interpolation='nearest', cmap='jet')
plt.yticks([])
plt.xlim(-27.5,-15)
plt.ylim(-32.5,-15)
plt.xlabel('M$_i - 5\log h $')
cb = plt.colorbar(format=r'$%.1f$', orientation='horizontal')
cb.set_label(r'$z$')
#plt.clim(-2.5, 0.5)

# Draw density contours over the colors
levels = np.linspace(0, np.log10(N.max()), 7)[3:]
CR = plt.contour(np.log10(N.T), levels, colors='k', linewidths=2.,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.tight_layout()
#plt.savefig('../MagB_Magi_z.png')
plt.show()


# In[119]:


red_data = []
for i in range(len(N)):
    for j in range(len(N)):
        if N[i,j]>10.:
            a = Imag[(Imag>xe[i])&(Imag<xe[i+1])&((MB>ye[j])&(MB<ye[j+1]))]
            b = MB[(Imag>xe[i])&(Imag<xe[i+1])&((MB>ye[j])&(MB<ye[j+1]))]
            c = zB[(Imag>xe[i])&(Imag<xe[i+1])&((MB>ye[j])&(MB<ye[j+1]))]
            red_data.append([a,b,c])


# In[120]:


temp_data = []
for i in red_data:
    temp_data.append(np.array(i).T)


# In[121]:


gdata = np.concatenate(temp_data)


# In[122]:


f,ax = plt.subplots(1,1,figsize=(8,6))
c=ax.scatter(gdata[:,0],gdata[:,1],c=gdata[:,2],s=0.5,cmap='jet')

#ax.plot(zrange,lc,'r-')
#ax.set_ylim(16,25)
#plt.gca().invert_yaxis()
ax.set_xlabel('M$_i - 5\log h $')
ax.set_ylabel(r'M$_{Blue} $')
#ax.set_ylim(-14,-24)
clb = plt.colorbar(c)
clb.set_label(r'$z$')
#plt.savefig('/cosma/home/durham/jarmijo/PAU_test/MagB_Magi_z_scatter.png')
#plt.gca().invert_yaxis()
plt.show()


# In[123]:


f,ax = plt.subplots(1,1,figsize=(8,6))
c=ax.scatter(gdata[:,0],gdata[:,1]-gdata[:,0],c=gdata[:,2],s=0.5,cmap='jet')

#ax.plot(zrange,lc,'r-')
#ax.set_ylim(16,25)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
ax.set_xlabel('M$_i - 5\log h $')
ax.set_ylabel(r'M$_{Blue}$ - M$_i$')
#ax.set_ylim(-14,-24)
clb = plt.colorbar(c)
clb.set_label(r'$z$')
#plt.savefig('/cosma/home/durham/jarmijo/PAU_test/MagB-Magi_magi_zc_1.png')
#plt.gca().invert_yaxis()
plt.show()


# In[124]:


f,ax = plt.subplots(1,1,figsize=(8,6))
c=ax.scatter(gdata[:,2],gdata[:,1]-gdata[:,0],c='k',s=0.5)

#ax.plot(zrange,lc,'r-')
#ax.set_ylim(16,25)
#plt.gca().invert_yaxis()
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'M$_{Blue}$ - M$_i$')
#ax.set_ylim(-14,-24)
#clb = plt.colorbar(c)
#clb.set_label(r'$z$',fontsize=fs)
#plt.savefig('/cosma/home/durham/jarmijo/PAU_test/MagB-Magi_z_1.png')
#plt.gca().invert_yaxis()
plt.show()


# In[126]:


z = gdata[:,2]
I = gdata[:,0]
B = gdata[:,1]
B_I = gdata[:,1]-gdata[:,0]


# In[135]:


zi


# In[328]:


# %load ../codes/plot_arrays_boxes.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch 
bl = 0.11
bu = 0.9
Nb = 4
BINS = np.linspace(bl,bu,Nb)#
b = BINS[:-1] + np.diff(BINS)[0]/2.
z0 = 0.11
sz = (bu-bl)/float(Nb)
f,ax = plt.subplots(1,4,figsize=(3.5*Nb,4.),sharex=True,
                        sharey=True,gridspec_kw={'width_ratios': [1,1,1,1]})
for i in range(Nb):
    zi = z0 + sz*i
    zf = zi + sz
    xx = I[(z>zi)&(z<zf)]
    yy = B[(z>zi)&(z<zf)]
    N, xedges, yedges = binned_statistic_2d(xx,yy,xx,'count', bins=50,range=[(-25,-16),(-30,-18)])
    xb = xedges[:-1]+np.diff(xedges)[0]/2.;yb = yedges[:-1]+np.diff(yedges)[0]/2.
    #
    ax[i].contourf(xb,yb,np.log10(N.T),cmap='PuBu_r',origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1],],levels=[1.,1.5,2.0,2.6],zorder=10)
    ax[i].scatter(xx,(yy),c='b',s=1.,label='%.2lf < z < %.2lf'%(zi,zf),zorder=0)
    #ax[i].imshow(np.log10(N.T), origin='lower',
    #       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
     #      aspect='auto', interpolation='nearest', cmap='viridis',vmin=1.3,vmax=2.5)
    ax[i].set_ylim(-19,-29)
    ax[i].set_xlim(-15.5,-25)
    ax[i].legend(loc=4,prop={'size':12})
    ax[i].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
ax[1].set_xticks([-16,-18,-20,-22,-24])
ax1 = f.add_axes([0.99, 0.23, 0.01, 0.72])

cmap = mpl.cm.PuBu_r
norm = mpl.colors.Normalize(vmin=1.3, vmax=2.6)

cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical',ticks=[1.0,1.5,2.,2.5])
cb1.set_label(r'$\log\ N$')
ax[0].set_ylabel(r'$M_{Blue}$')
ax[1].set_xlabel(r'$M_{i}$')
ax[1].xaxis.set_label_coords(1.1,-0.15)
plt.tight_layout()
plt.subplots_adjust(wspace=0.00)
#plt.savefig('/cosma/home/durham/jarmijo/PAU_test/Notebooks/MB_MI_4zbins_conts.png',bbox_inches='tight')
plt.show()

