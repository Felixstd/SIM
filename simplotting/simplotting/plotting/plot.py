import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import scienceplots

plt.style.use('science')

def normalize_stresses(sigI, sigII, p, eps = 1e-12):
    
    sigI_norm = np.full_like(sigI, np.nan)
    sigII_norm = np.full_like(sigII, np.nan)

    # Create a mask where p >= 1e-12
    # print(np.where(p >= 1e-20))
    mask = p >= 1e-12

    # Apply the mask to compute the normalized values
    sigI_norm[mask] = sigI[mask] / p[mask]
    sigII_norm[mask] = sigII[mask] / p[mask]
    
    return sigI_norm, sigII_norm


def plot_colormesh(ax, fig, X, Y, Field, cmap, norm, label, xlabel, ylabel):
    
        ax.set_aspect('equal', adjustable='box')
        pc = ax.pcolormesh(Y, X, Field, cmap = cmap, norm = norm)
        fig.colorbar(pc, ax = ax, label = label)
        
        if xlabel:
            ax.set_xlabel(xlabel)
        
        if ylabel:
            ax.set_ylabel(ylabel)
            
def histograms(data_dict, mu_0, mu_infty, dates, figdir, expno):
        
    _, _, _, _, _, _, _, mu_tot, _, I_tot = \
    data_dict.values()

    for k, date in enumerate(dates):
        
        mu = mu_tot[k]
        I = I_tot[k]

        
        min_I, max_I = np.min(I), np.max(I)
        bins_I = np.linspace(0, 1, 1000)
        
        bins_mu = np.linspace(mu_0, mu_infty, 200)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6))
        

        ax1.hist(I.flatten(), bins_I, color = 'b')
        ax1.set_yscale('log')
        # ax1.set_xscale('log')
        ax1.set_xlabel(r'$I$')
        ax1.set_ylabel('Counts')
        
        ax2.hist(mu.flatten(), bins_mu, color = 'b')
        ax2.set_xlabel(r'$\mu(I)$')

        plt.savefig(figdir+expno+'/histogram_mu_I{}.png'.format(date))
    
        
def uniaxial(dates, expno, data_dict, dx, figdir, MuPhi = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, muI_tot, phi_tot, I_tot = \
        data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot = \
            data_dict.values()
        
    Nx, Ny = np.shape(divergence_tot[0])
    
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        
        
        divergence, shear = divergence_tot[k], shear_tot[k]
        h, A, p = h_tot[k], A_tot[k], p_tot[k]
        sigI, sigII = sigI_tot[k], sigII_tot[k]
        # print(p)
        
        sigI_norm, sigII_norm = normalize_stresses(sigI, sigII, p)
        
        #----- Deformation Plots -----#
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 6))
        plot_colormesh(ax1, fig, X, Y, 1-A, cmocean.cm.ice, colors.LogNorm(vmin=1e-8, vmax=1e-6), r'$1-A$', None, 'y (km)')
        plot_colormesh(ax2, fig, X, Y, 1-h, cmocean.cm.ice, colors.LogNorm(vmin=1e-8, vmax=1e-6), r'$1-h$', None, None)
        plot_colormesh(ax3, fig, X, Y, divergence, cmocean.cm.thermal, colors.Normalize(vmin=1e-3, vmax=1e-1), 
            r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', 'x (km)', 'y (km)')
        plot_colormesh(ax4, fig, X, Y, shear, cmocean.cm.thermal, colors.LogNorm(vmin=1e-5, vmax=1e-2), 
            r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', 'x (km)', None)       
        fig.tight_layout()
        plt.savefig(figdir+expno+'/deformations_{}.png'.format(date))
        
        #------ Stresses plots ------#
        
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(9, 6))
        plot_colormesh(ax1, fig, X, Y, sigI_norm, cmocean.cm.ice, colors.Normalize(vmin=-1, vmax=1), r'$\sigma_{\mathrm{I}}/P$', 'x (km)', 'y (km)')
        plot_colormesh(ax2, fig, X, Y, sigII_norm, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1), r'$\sigma_{\mathrm{II}}/P$', 'x (km)', None)      
        fig.tight_layout()
        plt.savefig(figdir+expno+'/stresses_invariant{}.png'.format(date))
        # print(h[:10, :10], p[:10, :10])
        #---- Mu-Phi variables ----#
        if MuPhi:
            
            muI, phi, I = muI_tot[k], phi_tot[k], I_tot[k]
            
            
            fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
            
            plot_colormesh(ax1, fig, X, Y, I, cmocean.cm.amp, colors.LogNorm(vmin=1e-5, vmax=1), r'$I$', 'x (km)', 'y (km)')
            plot_colormesh(ax2, fig, X, Y, phi, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1), r'$\Phi(I)$', 'x (km)', None)
            plot_colormesh(ax3, fig, X, Y, muI, cmocean.cm.ice, colors.LogNorm(vmin=0.1, vmax=0.9), r'$\mu(I)$ ', 'x (km)', None)

            fig.tight_layout()
            
            plt.savefig(figdir+expno+'/mu_phi{}.png'.format(date))
            
            
            fig = plt.figure()
            ax = plt.axes()
            plt.axis('equal')
            sc = ax.scatter(sigI_norm.flatten(), sigII_norm.flatten(), s = 2, c = I.flatten(), cmap = plt.cm.magma, norm = colors.LogNorm(vmin=1e-5, vmax=1e-3))
            plt.grid('--')
            fig.colorbar(sc, label = 'I')
            ax.set_xlabel(r'$\sigma_{I}/P$')
            ax.set_ylabel(r'$\sigma_{II}/P$')
            plt.savefig(figdir+expno+'/stress_space{}.png'.format(date))
            
        plt.close()
            
        
        
        
        
    
    
    
    