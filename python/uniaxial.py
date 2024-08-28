#%%
import numpy as np
import matplotlib.pyplot as plt
import cmocean 
import matplotlib.colors as colors

plt.style.use('figure_style.mplstyle')

exp = '57'
ExpDir = '/aos/home/fstdenis/SIM/Experiments/'

# dates = ['1900_01_01_00_00_00', '1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04','1990_01_01_00_00_05',  \
#     '1990_01_01_00_00_06', '1990_01_01_00_00_07', '1990_01_01_00_00_08', '1990_01_01_00_00_09', '1990_01_01_00_00_10', '1990_01_01_00_00_11']

dates = ['1990_01_01_00_15_00', '1990_01_01_00_30_00', '1990_01_01_00_45_00', '1990_01_01_01_00_00', '1990_01_01_01_15_00','1990_01_01_01_30_00',  \
    '1990_01_01_01_45_00', '1990_01_01_02_00_00', '1990_01_01_02_05_00', '1990_01_01_02_10_00', '1990_01_01_02_15_00', '1990_01_01_02_20_00', '1990_01_01_02_25_00']

# dates = [ '1900_01_01_00_00_00', '1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04','1990_01_01_00_00_05']
# outputdir = "/aos/home/fstdenis/SIM-Mu-Phi/output/"
outputdir = "/storage/fstdenis/output_sim/"


dx = 1e3
dt = 1
# Nx, Ny = 100, 250
Nx, Ny = 200, 500
# Nx, Ny = 102, 402

X = np.arange(0, Nx+2)*dx/1e3
Y = np.arange(0, Ny+2)*dx/1e3


# VP = True
MuPhi = True

k = 1
t = 0
for date in dates:
    
    file_divergence = outputdir+'div'+ date +'_k{:04d}.'.format(k)+exp
    file_shear = outputdir+'shear'+ date +'_k{:04d}.'.format(k)+exp
    file_h = outputdir+'h' + date + "." +exp
    file_A = outputdir+'A' + date + "." +exp
    file_p = outputdir+'p' + date + "." +exp
    file_sig_1 = outputdir+'sigI'+ date +'_k{:04d}.'.format(k)+exp
    file_sig_2 = outputdir+'sigII'+ date +'_k{:04d}.'.format(k)+exp
    
    h = np.loadtxt(file_h, dtype=None)
    A = np.loadtxt(file_A, dtype=None)
    shear = np.loadtxt(file_shear, dtype=None)
    divergence = np.loadtxt(file_divergence, dtype=None)
    p = np.loadtxt(file_p, dtype=None)
    sig_I = np.loadtxt(file_sig_1)
    sig_II = np.loadtxt(file_sig_2)
    
    shear[shear == -999] = np.nan
    sig_I[sig_I == -999] = 0
    sig_II[sig_II == -999] = 0
    
    if MuPhi:
        
        file_mu = outputdir+'mu_I' + date + "." +exp
        file_I = outputdir+'I' + date + "." +exp
        file_phi = outputdir+'phi_I' + date + "." +exp
        file_shear = outputdir+'e_II' + date + "." +exp
        
# 
        mu = np.loadtxt(file_mu, dtype=None)
        I = np.loadtxt(file_I, dtype=None)
        phi = np.loadtxt(file_phi, dtype=None)
    

    
    
    if k>0:
        sigII_norm = np.zeros_like(sig_I)   
        sigI_norm = np.zeros_like(sig_I) 
        
        # print(Nx)
        for i in range(0, Ny+2):
            for j in range(0, Nx+2):
                # print(i,j)
                    
                if p[i, j] < 1e-12:
                    sigI_norm[i,j] = np.nan
                    sigII_norm[i,j] = np.nan
                    
                else:
                    sigI_norm[i, j] = sig_I[i, j]/p[i, j]
                    sigII_norm[i, j] = sig_II[i, j]/p[i, j]

    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 6))
    
    axes = [ax1, ax2, ax3, ax4]
    
    for ax in axes:
        ax.set_aspect('equal', adjustable='box')
 

    
    pc = ax1.pcolormesh(X, Y, 1-A, cmap = cmocean.cm.ice, norm = colors.LogNorm(vmin=1e-8, vmax=1e-6))
                        
                        # vmin = -1e8, vmax = 1e-8)
    fig.colorbar(pc, ax = ax1, label = r'$1-A$ ')
    ax1.set_ylabel('y (km)')
    
    pc = ax2.pcolormesh(X, Y, 1-h, cmap = cmocean.cm.ice, norm=colors.LogNorm(vmin=1e-8, vmax=5e-6))
    fig.colorbar(pc, ax = ax2, label = r'$1 - h$ (m)')
    
    pc = ax3.pcolormesh(X, Y, divergence, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=5e-4, vmax=1e-2))
    fig.colorbar(pc, ax = ax3, label = r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$ ')
    ax3.set_ylabel('y (km)')
    ax3.set_xlabel('x (km)')
    
    pc = ax4.pcolormesh(X, Y, shear, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=5e-4, vmax=1e-2) )
    fig.colorbar(pc, ax = ax4,label = r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$')
    ax4.set_xlabel('x (km)')
    
    # fig.supxlabel('x (km)')
    # fig.supylabel('y (km)')
    
    fig.suptitle('$t = {}$ s'.format(t))
    fig.tight_layout()
    
    plt.savefig(ExpDir+exp+'/mu_uniaxial_{}.png'.format(date))
    
    
    fig = plt.figure(figsize = (4,5))
    ax = plt.axes()
    ax.set_aspect('equal', adjustable='box')
    pc = ax.pcolormesh(X, Y, 1-h, cmap = cmocean.cm.ice, norm=colors.Normalize(vmin=1e-8, vmax=1e-6))
    fig.colorbar(pc, ax = ax, label = r'$1 - h$ (m)')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    plt.savefig(ExpDir+exp+'/h_{}.png'.format(date))
    
    zeta = np.zeros_like(I)
    for i in range(0, Ny):
        for j in range(0, Nx):
            
            if shear[i, j] < 1e-13:
                zeta[i, j] = np.nan
            
            else:
                
                zeta[i, j] = p[i, j]*mu[i, j]/(shear[i, j]/3600*24)
    
    plt.figure()
    plt.scatter(I.flatten(), zeta.flatten()/np.nanmax(zeta.flatten()))
    plt.grid()
    plt.xlabel('I')
    plt.ylim(0, 1e-5)
    plt.ylabel(r'$\frac{\mu}{\dot{\epsilon}_{\mathrm{II}}}p_{eq}$')
    plt.savefig(ExpDir+exp+'/zeta_tes_{}.png'.format(date))
    
    if k>0:
        fig = plt.figure()
        ax = plt.axes()
        # ax.set_aspect('equal', adjustable='box')
        ax.scatter(sigI_norm.flatten(), sigII_norm.flatten(), s = 2, color = 'k')
        plt.grid('--')
        # plt.xlim(-800, -1050)
        ax.set_xlabel(r'$\sigma_{I}/P$')
        ax.set_ylabel(r'$\sigma_{II}/P$')
        plt.savefig(ExpDir+exp+'/stresses_{}.png'.format(date))
    
    
    plt.clf()
    
    
    
    k+=1
    t+=1
    
# %%
