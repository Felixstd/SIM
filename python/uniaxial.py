#%%
import numpy as np
import matplotlib.pyplot as plt
import cmocean 
import matplotlib.colors as colors

plt.style.use('figure_style.mplstyle')

exp = '34'
dates = ['1900_01_01_00_00_00', '1990_01_01_00_00_01']

# outputdir = "/aos/home/fstdenis/SIM-Mu-Phi/output/"
outputdir = "/storage/fstdenis/output_sim/"


dx = 25
Nx, Ny = 400, 1000

# Nx, Ny = 102, 402

X = np.arange(0, Nx+2)*dx/1e3
Y = np.arange(0, Ny+2)*dx/1e3


# VP = True
MuPhi = False

k = 0
for date in dates:
    
    file_divergence = outputdir+'div'+ date +'_k{:04d}.'.format(k)+exp
    file_shear = outputdir+'shear'+ date +'_k{:04d}.'.format(k)+exp
    file_h = outputdir+'h' + date + "." +exp
    file_A = outputdir+'A' + date + "." +exp
    
    
    h = np.loadtxt(file_h, dtype=None)
    A = np.loadtxt(file_A, dtype=None)
    shear = np.loadtxt(file_shear, dtype=None)
    divergence = np.loadtxt(file_divergence, dtype=None)
    
    if MuPhi:
        
        file_mu = outputdir+'mu_I' + date + "." +exp
        file_I = outputdir+'I' + date + "." +exp
        file_phi = outputdir+'phi_I' + date + "." +exp
        file_shear = outputdir+'e_II' + date + "." +exp
# 
        mu = np.loadtxt(file_mu, dtype=None)
        I = np.loadtxt(file_I, dtype=None)
        phi = np.loadtxt(file_phi, dtype=None)
    
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 6))
    
    axes = [ax1, ax2, ax3, ax4]
    
    for ax in axes:
        ax.set_aspect('equal', adjustable='box')
 
    shear[shear == -999] = np.nan
    
    
    pc = ax1.pcolormesh(X, Y, A, cmap = cmocean.cm.ice, vmin = 0, vmax = 1)#norm=colors.LogNorm(vmin=0, vmax=1))
    fig.colorbar(pc, ax = ax1, label = r'$A$ ')
    ax1.set_ylabel('y (m)')
    
    pc = ax2.pcolormesh(X, Y, 1-h, cmap = cmocean.cm.ice)#, norm=colors.LogNorm(vmin=0.9, vmax=1))
    fig.colorbar(pc, ax = ax2, label = r'$1 - h$ (m)')
    
    pc = ax3.pcolormesh(X, Y, divergence, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=1e-3, vmax=1e-1))
    fig.colorbar(pc, ax = ax3, label = r'$\dot{\epsilon}_{\mathrm{I}} \text{ s}^{-1}$ ')
    ax3.set_ylabel('y (m)')
    ax3.set_xlabel('x (m)')
    
    pc = ax4.pcolormesh(X, Y, shear, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=1e-3, vmax=1e-1) )
    fig.colorbar(pc, ax = ax4,label = r'$\dot{\epsilon}_{\mathrm{II}} \text{ s}^{-1}$')
    ax4.set_xlabel('x (m)')
    
    # fig.supxlabel('x (km)')
    # fig.supylabel('y (km)')
    
    
    plt.savefig('mu_uniaxial_{}.png'.format(date))
    k+=1
    
# %%
