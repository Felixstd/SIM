import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import scienceplots

plt.style.use('science')

def uniaxial(dates, expno, data_dict, dx, figdir, MuPhi = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, muI_tot, phi_tot, I_tot = \
        data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot = \
            data_dict.values()
        
    Nx, Ny = np.shape(divergence_tot[0])
    print(Ny)
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    
    print(Nx, Ny)
    
    for k, date in enumerate(dates):
        
        divergence, shear = divergence_tot[k], shear_tot[k]
        h, A, p = h_tot[k], A_tot[k], p_tot[k]
        sigI, sigII = sigI_tot[k], sigII_tot[k]
        muI, phi, I = muI_tot[k], phi_tot[k], I_tot[k]
        
        print(np.shape(p))
        sigII_norm = np.zeros_like(sigI)   
        sigI_norm = np.zeros_like(sigI) 
        
        for i in range(0, Nx):
            for j in range(0, Ny):

                    
                if p[i, j] < 1e-12:
                    sigI_norm[i,j] = np.nan
                    sigII_norm[i,j] = np.nan
                    
                else:
                    sigI_norm[i, j] = sigI[i, j]/p[i, j]
                    sigII_norm[i, j] = sigII[i, j]/p[i, j]
        
        
        #----- Deformation Plots -----#
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 6))
    
        axes = [ax1, ax2, ax3, ax4]
    
        for ax in axes:
            ax.set_aspect('equal', adjustable='box')
    
        pc = ax1.pcolormesh(Y, X, 1-A, cmap = cmocean.cm.ice, norm = colors.LogNorm(vmin=1e-8, vmax=1e-6))
                            
                            # vmin = -1e8, vmax = 1e-8)
        fig.colorbar(pc, ax = ax1, label = r'$1-A$ ')
        ax1.set_ylabel('y (km)')
        
        pc = ax2.pcolormesh(Y, X, 1-h, cmap = cmocean.cm.ice, norm=colors.LogNorm(vmin=1e-8, vmax=5e-6))
        fig.colorbar(pc, ax = ax2, label = r'$1 - h$ (m)')
        
        pc = ax3.pcolormesh(Y, X, divergence, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=5e-4, vmax=1e-2))
        fig.colorbar(pc, ax = ax3, label = r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$ ')
        ax3.set_ylabel('y (km)')
        ax3.set_xlabel('x (km)')
        
        pc = ax4.pcolormesh(Y, X, shear, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=5e-4, vmax=1e-2) )
        fig.colorbar(pc, ax = ax4,label = r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$')
        ax4.set_xlabel('x (km)')
        
        fig.tight_layout()
        
        plt.savefig(figdir+expno+'/deformations_{}.png'.format(date))
        
        #---- Mu-Phi variables ----#
        if MuPhi:
            
            fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(8, 6))
    
            axes = [ax1, ax2, ax3]
        
            for ax in axes:
                ax.set_aspect('equal', adjustable='box')
        
            pc = ax1.pcolormesh(Y, X, I, cmap = cmocean.cm.amp, norm = colors.LogNorm(vmin=1e-8, vmax=1e-6))
                                
            fig.colorbar(pc, ax = ax1, label = r'$I$ ')
            ax1.set_ylabel('y (km)')
            
            pc = ax2.pcolormesh(Y, X, phi, cmap = cmocean.cm.ice, norm=colors.LogNorm(vmin=1e-8, vmax=5e-6))
            fig.colorbar(pc, ax = ax2, label = r'$\Phi(I)$')
            
            pc = ax3.pcolormesh(Y, X, muI, cmap = cmocean.cm.thermal, norm=colors.LogNorm(vmin=5e-4, vmax=1e-2))
            fig.colorbar(pc, ax = ax3, label = r'$\mu(I)$ ')
            ax3.set_ylabel('y (km)')
            ax3.set_xlabel('x (km)')
            
            fig.tight_layout()
            
            plt.savefig(figdir+expno+'/mu_phi{}.png'.format(date))
            
        
        
        
        
    
    
    
    