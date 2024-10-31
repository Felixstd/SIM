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
        # ax.invert_yaxis()
        if xlabel:
            ax.set_xlabel(xlabel)
        
        if ylabel:
            ax.set_ylabel(ylabel)
            
def histograms(ax, var, mu0, mu_infty,varname, I = False, I_0 = 0):
        

    min, max = np.min(var), np.max(var)
    bins = np.linspace(min, max, 1000)
        
        
    # fig, ax1 = plt.subplots(1, 1, figsize = (12, 6))
    
    if I: 
        var_nonmax=var[var < max] 
        max = np.max(var_nonmax)
        bins = np.linspace(min, max, 1000)
        ax.axvline(I_0, color = 'k', label = r'$I_0$')
        plt.legend()
        

    ax.hist(var.flatten(), bins, color = 'b')
    ax.set_xlim(mu0-0.05, 0.15)
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlabel(varname)
    ax.set_ylabel('Counts')

    # plt.savefig(figdir+expno+'/'+filename+'{}.png'.format(date))
    
        
def uniaxial(dates, expno, data_dict, dx, figdir, mu_0, mu_infty, MuPhi = True, log = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, uair_tot, vair_tot = \
            data_dict.values()
        
    Nx, Ny = np.shape(divergence_tot[0])
    
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        divergence, shear = divergence_tot[k], shear_tot[k]
        h, A, p = h_tot[k], A_tot[k], p_tot[k]
        sigI, sigII = sigI_tot[k], sigII_tot[k]
        zeta = zeta_tot[k]
        uair, vair = uair_tot[k], vair_tot[k]
        
        sigI_norm, sigII_norm = normalize_stresses(sigI, sigII, p)
        
        #----- Deformation Plots -----#
        x, y = np.meshgrid(Y, X)
        # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 6))
        quiver_step = max(x.shape[1] // 10, 1)
        quiver_step_2 = max(x.shape[0] // 10, 1)
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 4))
        plot_colormesh(ax1, fig, X, Y, 1-A, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1e-7), r'$1-A$', None, 'y (km)')
        ax1.quiver(x[::quiver_step_2, ::quiver_step], y[::quiver_step_2, ::quiver_step], uair[::quiver_step_2, ::quiver_step], vair[::quiver_step_2, ::quiver_step], color = 'r')
        plot_colormesh(ax2, fig, X, Y, 1-h, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1e-7), r'$1-h$', None, None)
        if log:
            plot_colormesh(ax3, fig, X, Y, divergence, cmocean.cm.thermal, colors.LogNorm(vmin=1e-4, vmax=1e-1), 
            r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', 'x (km)', 'y (km)')
            plot_colormesh(ax4, fig, X, Y, shear, cmocean.cm.thermal, colors.LogNorm(vmin=1e-4, vmax=1e-1), 
            r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', 'x (km)', None) 
        else:
            plot_colormesh(ax3, fig, X, Y, divergence, cmocean.cm.thermal, colors.Normalize(vmin=1e-4, vmax=1e-2), 
                r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', 'x (km)', 'y (km)')
            plot_colormesh(ax4, fig, X, Y, shear, cmocean.cm.thermal, colors.Normalize(vmin=1e-4, vmax=1e-2), 
                r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', 'x (km)', None) 
        
              
        fig.tight_layout()
        plt.savefig(figdir+expno+'/deformations_{}.png'.format(date))
        
        #------ Stresses plots ------#
        
        fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
        plot_colormesh(ax1, fig, X, Y, sigI_norm, cmocean.cm.ice, colors.Normalize(vmin=-1, vmax=1e-2), r'$\sigma_{\mathrm{I}}/P$', 'x (km)', 'y (km)')
        plot_colormesh(ax2, fig, X, Y, sigII_norm, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1e-2), r'$\sigma_{\mathrm{II}}/P$', 'x (km)', None)     
        plot_colormesh(ax3, fig, X, Y, zeta, cmocean.cm.ice, colors.Normalize(vmin=np.min(1e8), vmax=np.max(1e12)), r'$\zeta$', 'x (km)', None)   
        fig.tight_layout()
        plt.savefig(figdir+expno+'/stresses_invariant{}.png'.format(date))
        
        # print(h[:10, :10], p[:10, :10])
        #---- Mu-Phi variables ----#
        if MuPhi:
            
            muI, phi, I = muI_tot[k], phi_tot[k], I_tot[k]
            p_VP, p_muphi = Pmax_tot[k], Peq_tot[k]
            
            
            fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 6))
            
            # plot_colormesh(ax1, fig, X, Y, I, cmocean.cm.amp, colors.LogNorm(vmin=1e-5, vmax=1), r'$I$', 'x (km)', 'y (km)')
            # plot_colormesh(ax1, fig, X, Y, 1-phi, cmocean.cm.ice, colors.Normalize(vmin=1e-7, vmax=1e-6), r'$1-\Phi(I)$', 'x (km)', None)
            plot_colormesh(ax1, fig, X, Y, muI, cmocean.cm.ice, colors.LogNorm(vmin=mu_0, vmax=mu_infty), r'$\mu(I)$ ', 'x (km)', None)
            histograms(ax2, muI, mu_0, mu_infty, r'$\mu(I)$', I = False, I_0 = 1e-3)
            fig.tight_layout()
            
            plt.savefig(figdir+expno+'/mu_phi{}.png'.format(date))
            
            
            fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
            
            plot_colormesh(ax1, fig, X, Y, p_VP, cmocean.cm.amp, colors.Normalize(vmin=27.4e3, vmax=27.6e3), r'$P_{VP} (N/m)$', 'x (km)', 'y (km)')
            plot_colormesh(ax2, fig, X, Y, p_muphi, cmocean.cm.ice, colors.LogNorm(vmin=50e3, vmax=100e3), r'$P$ (N/m)', 'x (km)', None)
            plot_colormesh(ax3, fig, X, Y, np.abs(p_VP - p_muphi), cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1e3), r'$|P_{VP} - P|$ (N/m)', 'x (km)', None)
            fig.tight_layout()
            
            plt.savefig(figdir+expno+'/pressures_{}.png'.format(date))
        
            
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
            
        
        
        
def plot_single(Nx, Ny, dx, Var, label, figname):
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    x, y = np.meshgrid(Y, X)
    
    fig = plt.figure(figsize=(5, 4))
    ax1 = plt.axes()
    plot_colormesh(ax1, fig, Y, X, Var, cmocean.cm.ice, colors.LogNorm(vmin=1e-4, vmax=1e-2), label, 'x (km)', 'y (km)')
    plt.savefig(figname)


    