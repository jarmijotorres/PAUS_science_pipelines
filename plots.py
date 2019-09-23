import matplotlib.pyplot as plt
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('xtick',direction='inout')
plt.rc('ytick',direction='inout')
plt.rc('axes',linewidth=1.5)
plt.rc('font',family='sans-serif')
plt.rc('font',size=16)
#
nf = 3
nc = 4
f,ax = plt.subplots(nf,nc,figsize=(4*nf,4.),sharex=True,
                        sharey=True,gridspec_kw={'width_ratios': [1]*nc, 'height_ratios': [1]*nf})
for f in range(nf):
    for c in range(nc):
        Vr = L_Vgal[nc*f+c]/L_Vmax[nc*f+c]
        zi = zbins[nc*f+c]
        zf = zbins[nc*f + (c+1)]
        ax[f,c].hist(Vr,bins=20,range=(0,1),histtype='step',label = "%.2f < z < %.2f"%(zi,zf))
        ax[f,c].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
        ax[f,c].legend(prop={'size':10})
plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0.0)
plt.show()

# =============================================================================
nf = 3
nc = 4
f,ax = plt.subplots(nf,nc,figsize=(4*nf,4.),sharex=True,
                        sharey=True,gridspec_kw={'width_ratios': [1]*nc, 'height_ratios': [1]*nf})
for f in range(nf):
    for c in range(nc):
        Vr = L_Vgal[nc*f+c]/L_Vmax[nc*f+c]
        zi = zbins[nc*f+c]
        zf = zbins[nc*f + (c+1)]
        NVr,_ = np.histogram(Vr,bins=20,range=(0,1))
        N_bar = np.mean(NVr)
        bc_Vr = _[:-1] + np.diff(_)[0]/2.
        ax[f,c].step(bc_Vr,NVr/N_bar,where='post',color='b',linestyle='-',label='%.2lf < z < %.2lf'%(zi,zf))
        ax[f,c].hlines(1.0,0.,1.,linestyle='--',color='k',linewidth=1.5)
        #ax[f,c].hist(Vr,bins=20,range=(0,1),histtype='step',label = "%.2f < z < %.2f"%(zi,zf))
        ax[f,c].tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
ax[f,c].legend(prop={'size':10})

ax[0,0].set_xlim(0,1)
plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0.0)
plt.show()
# =============================================================================
############################
f,ax = plt.subplots(1,1,figsize=(7,6))
for c in range(len(L_M)):
    zi = zbins[c]
    zf = zbins[c+1]
    ax.hist(L_M[c],bins=30,range=(M_max,M_min),histtype='step',label = "%.2f < z < %.2f"%(zi,zf),linewidth=2.5)
ax.tick_params(direction='inout', length=8, width=2, colors='k',
               grid_color='k', grid_alpha=0.5)
ax.set_xlabel(r'$M_{i} - 5\log_{10}h$')
ax.set_ylabel('Counts')
ax.legend(prop = {'size':10},loc=2)
plt.tight_layout()
plt.show()
#==============================================

f,ax = plt.subplots(1,1,figsize=(7,6))
for c in range(len(L_LF)):
    zi = zbins[c]
    zf = zbins[c+1]
    ax.plot(bb,np.log10(L_LF[c]),'o-',c=np.random.random(3),ms=5.,label = "%.2f < z < %.2f"%(zi,zf))
ax.legend(prop={'size':12})
# =============================================================================
ax.set_xticks(np.arange(-24,-15,1))
ax.set_xlim(-24,-16)
ax.set_yticks(np.arange(-7,-1,1))
ax.set_ylim(-6.,-1.8)
ax.legend(prop = {'size':12})
ax.set_xlabel('$M_{i} - 5\log_{10}h$')
ax.set_ylabel('$\log$ [$\Phi$ Mpc$^3$/$h^{-3}$ (0.25 mag)$^{-1}]$')
#plt.savefig('./Dropbox/PhD/Durham/Projects/PAU/figures/July/Mblue_LF_mocks_mi23cut_8zbins_maxcut_2.png',bbox_inches='tight')
plt.tight_layout()
plt.show()

# =============================================================================
f,ax = plt.subplots(1,1,figsize=(7,6))
for i in range(1):
#    ax.plot(bb,np.log10(L_LF[i]/(1e2)**3),'-')
   c1 = ax.plot(bb,np.log10(sLF_mean[i]/(1e2)**3),'-',color=np.random.random(3),linewidth=1.5,label='$z = %.3lf$'%zs[i])
c2 = ax.plot(bb,np.log10(lcLF[0]),'b--')
ax.set_xticks(np.arange(-24,-15,1))
ax.set_xlim(-24,-16)
ax.set_yticks(np.arange(-7,-1,1))
ax.set_ylim(-6.,-1.8)
ax.set_xlabel('$M_{i} - 5\log_{10}h$')
ax.set_ylabel('$\log$ [$\Phi$ Mpc$^{-3}$/$h^{3}$ (0.25 mag)$^{-1}]$')
lines = ax.get_lines()
l1 = ax.legend(c1 ,[r"$z$ = %.3lf"%zs[i]],loc=4)
l2 = ax.legend(c2,["$0.1 < z < 0.2$"])
ax.add_artist(l1)
ax.add_artist(l2)
plt.tight_layout()
plt.savefig('/cosma5/data/dp004/dc-armi2/PAU/PAU_test/figures/LF_s_lc_Mm16m24_dM0.25_z57_zbin0_mi23cut.png',bbox_inches= 'tight')
plt.show()