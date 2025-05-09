import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import cmocean
import scienceplots
plt.style.use('science')


R_EARTH = 6370.0  # Earth's radius (is smaller for better looking plots)
BETA = 32.0  # Angle between domain and Greenwich
BETA_RGPS = 45.0  # Angle between domain and Greenwich
RES_RGPS = 12.5  # rgps resolution
PLANE_RGPS = 70.0  # rgps sterographic plane latitude

def dilatation(mu, angle_phi):
    
    return (mu - np.tan(angle_phi))/(1+mu*angle_phi) 

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

def plot_colormesh(ax, dx, fig, X, Y, Field, cmap, norm, label, xlabel, ylabel):
        Nx, Ny = np.shape(Field)
        ax.set_aspect('equal', adjustable='box')
        pc = ax.pcolormesh(Y, X, Field, cmap = cmap, norm = norm)
        # ax.axvline(20*dx/1e3, color = 'pink', zorder = 3)
        # ax.axvline(180*dx/1e3, color = 'pink', zorder = 3)
        ax.set_xticks(np.arange(0, Ny, 100)*dx/1e3)
        # cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        cbar = fig.colorbar(pc, ax = ax)
        cbar.set_label(label = label, fontsize = 16)
        # ax.invert_yaxis()
        if xlabel:
            ax.set_xlabel(xlabel)
        
        if ylabel:
            ax.set_ylabel(ylabel)
            
def histograms(ax, var, mu0, mu_infty,varname, I = False, I_0 = 0):
        

    min, max = np.min(var), np.max(var)
    bins = np.linspace(mu0, mu_infty, 100)
        
        
    # fig, ax1 = plt.subplots(1, 1, figsize = (12, 6))
    
    if I: 
        var_nonmax=var[var < max] 
        max = np.max(var_nonmax)
        bins = np.linspace(min, max, 1000)
        ax.axvline(I_0, color = 'k', label = r'$I_0$')
        plt.legend()
        

    ax.hist(var.flatten(), bins, color = 'b')
    ax.set_xlim(mu0-0.05, mu_infty+0.05)
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlabel(varname)
    ax.set_ylabel('Counts')

    # plt.savefig(figdir+expno+'/'+filename+'{}.png'.format(date)) 
        
def uniaxial(dates, expno, data_dict, dx, figdir, mu_0, mu_infty,angle_phi, MuPhi = True, log = True, arrows = True, Shear = True):
    
    # ------ Reading the data ------# 
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot,sigI_tot, sigII_tot, zeta_tot, eta_tot, uair_tot, vair_tot = \
            data_dict.values()
        
    Nx, Ny = np.shape(divergence_tot[0])
    
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    
    Nx2, Ny2 = np.shape(u_tot[0])
    

    X2, Y2 =  np.arange(0, Nx2)*dx/1e3,  np.arange(0, Ny2-1)*dx/1e3
    
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        divergence = divergence_tot[k]
        if MuPhi:
            shear = shearI_tot[k]*86400
        else:
            shear = shear_tot[k]
        # shear = shearI_tot[k]*86400
        h, A, p = h_tot[k], A_tot[k], p_tot[k]
        sigI, sigII = sigI_tot[k], sigII_tot[k]
        zeta, eta = zeta_tot[k], eta_tot[k]
        uair, vair = uair_tot[k], vair_tot[k]
        u, v = u_tot[k], v_tot[k]
        

        divergence[divergence < -100] = np.nan
        shear[shear < 0 ] = np.nan
        sigI[sigI < -500 ] = np.nan
        sigII[sigII < -500 ] = np.nan
        
        
        #----- DEFORMATION Plots -----#
        
        x, y = np.meshgrid(Y, X)
        x2, y2 = np.meshgrid(Y2, X2)
        
        # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(5, 6))
        quiver_step = max(x2.shape[1] // 20, 1)
        quiver_step_2 = max(x2.shape[0] // 20, 1)
        
        
        if Shear:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(10, 6))
        else:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(6, 7))
        plot_colormesh(ax1, dx, fig, X, Y, 1-A, cmocean.cm.ice, colors.LogNorm(vmin=1e-6, vmax=1e-1), r'$1-A$', None, None)
        plot_colormesh(ax2, dx, fig, X, Y, h-1, cmocean.cm.balance, colors.SymLogNorm(linthresh  = 1e-4, vmin=-1e-2, vmax=1e-2), r'$h-1$', None, None)
    
        if log:
            plot_colormesh(ax3, dx, fig, X, Y, divergence, cmocean.cm.balance, colors.LogNorm(vmin=1e-4, vmax=1e-0), 
            r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', None, None)
            plot_colormesh(ax4, dx, fig, X, Y, shear, cmocean.cm.amp, colors.LogNorm(vmin=1e-4, vmax=1e-0), 
            r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', None, None) 
        else:
            plot_colormesh(ax3, dx, fig, X, Y, divergence, cmocean.cm.balance, colors.SymLogNorm(1e-6, vmin=-1e-1, vmax=1e-1), 
                r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', None, None)
            plot_colormesh(ax4, dx, fig, X, Y, shear, cmocean.cm.amp, colors.LogNorm(vmin=1e-7, vmax=1e-1), 
                r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', None, None) 
        
        plt.subplots_adjust(hspace = 0.1)
        fig.supxlabel('x (km)')
        fig.supylabel('y (km)')
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/deformations_{}.png'.format(date))
        
        divergence = divergence[1:-1, 1:-1]
        shear = shear[1:-1, 1:-1]
        h = h[1:-1, 1:-1]
        A = A[1:-1, 1:-1]
        uair = uair[:, :-1]
        u = u[:, :-1]
        vair = vair[:-1, :]
        v = v[:-1, :]
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex = True, sharey = True, figsize = (8, 4))
        
        
        plot_colormesh(ax1, dx, fig, X2, Y2, divergence, cmocean.cm.balance, colors.SymLogNorm(1e-4, vmin=-1e-1, vmax=1e-1), 
                r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', None, None)
        ax1.quiver(x2[::quiver_step_2, ::quiver_step], y2[::quiver_step_2, ::quiver_step], u[::quiver_step_2, ::quiver_step], v[::quiver_step_2, ::quiver_step], color = 'r')
        
        plot_colormesh(ax2, dx, fig, X2, Y2, shear, cmocean.cm.amp, colors.LogNorm(vmin=1e-4, vmax=1e-0), 
            r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', None, None) 
        ax2.quiver(x2[::quiver_step_2, ::quiver_step], y2[::quiver_step_2, ::quiver_step], u[::quiver_step_2, ::quiver_step], v[::quiver_step_2, ::quiver_step], color = 'r')
        
        
        plot_colormesh(ax3, dx, fig, X2, Y2, divergence, cmocean.cm.balance, colors.SymLogNorm(1e-4, vmin=-1e-1, vmax=1e-1), 
                r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', None, None)
        ax3.quiver(x2[::quiver_step_2, ::quiver_step], y2[::quiver_step_2, ::quiver_step], uair[::quiver_step_2, ::quiver_step], vair[::quiver_step_2, ::quiver_step], color = 'g')
        
        plot_colormesh(ax4, dx, fig, X2, Y2, shear, cmocean.cm.amp, colors.LogNorm(vmin=1e-4, vmax=1e-0), 
            r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', None, None) 
        ax4.quiver(x2[::quiver_step_2, ::quiver_step], y2[::quiver_step_2, ::quiver_step], uair[::quiver_step_2, ::quiver_step], vair[::quiver_step_2, ::quiver_step], color = 'g')
       
        fig.supylabel('Y (km)')
        fig.supxlabel('X (km)')
        plt.savefig(figdir+expno+'/velocity_{}.png'.format(date))
        
        #------ Stresses plots ------#
        if Shear:
            fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 3))
        else:
            fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
        # plot_colormesh(ax1, fig, X, Y, sigI_norm, cmocean.cm.ice, colors.Normalize(vmin=-1, vmax=1e-2), r'$\sigma_{\mathrm{I}}/P$', 'x (km)', 'y (km)')
        # plot_colormesh(ax2, fig, X, Y, sigII_norm, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1e-2), r'$\sigma_{\mathrm{II}}/P$', 'x (km)', None)   
        plot_colormesh(ax1, dx, fig, X, Y, sigI, cmocean.cm.ice, colors.Normalize(vmin=-1, vmax=1e-2), r'$\sigma_{\mathrm{I}}/P$', 'x (km)', 'y (km)')
        plot_colormesh(ax2, dx, fig, X, Y, sigII, cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1), r'$\sigma_{\mathrm{II}}/P$', 'x (km)', None)       
        plot_colormesh(ax3, dx, fig, X, Y, zeta, cmocean.cm.ice, colors.Normalize(vmin=np.min(1e8), vmax=np.max(1e12)), r'$\zeta$', 'x (km)', None)   
        fig.tight_layout()
        plt.savefig(figdir+expno+'/stresses_invariant{}.png'.format(date))
        
        # print(h[:10, :10], p[:10, :10])
        #---- Mu-Phi variables ----#
        if MuPhi:
            
            muI, phi, I = muI_tot[k], phi_tot[k], I_tot[k]
            p_VP, p_muphi = Pmax_tot[k], Peq_tot[k]
            
            
            if Shear:
                fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
                plot_colormesh(ax1, dx, fig, X, Y, muI, cmocean.cm.ice, colors.LogNorm(vmin=mu_0, vmax=mu_infty), r'$\mu(I)$ ', 'x (km)', 'y (km)')
                plot_colormesh(ax2, dx, fig, X, Y, I, cmocean.cm.amp, colors.LogNorm(vmin=1e-6, vmax=1e-2), r'$I$ ', 'x (km)', None)
                histograms(ax3, muI, mu_0, mu_infty, r'$\mu(I)$', I = False, I_0 = 1e-3)
            else:
                fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
                plot_colormesh(ax1, dx, fig, X, Y, muI, cmocean.cm.ice, colors.LogNorm(vmin=mu_0, vmax=mu_infty), r'$\mu(I)$ ', 'x (km)', None)
                plot_colormesh(ax2, dx, fig, X, Y, dilatation(muI, angle_phi), cmocean.cm.ice, colors.SymLogNorm(linthresh = 1e-3, vmin=dilatation(mu_0, angle_phi), 
                                vmax=dilatation(mu_infty, angle_phi)), r'$\tan \psi$ ', 'x (km)', None)
                histograms(ax3, muI, mu_0, mu_infty, r'$\mu(I)$', I = False, I_0 = 1e-3)
            
            fig.tight_layout()
            
            plt.savefig(figdir+expno+'/mu_phi{}.png'.format(date))
            
            
            if Shear:
                fig, (ax1) = plt.subplots(1,1, figsize=(5, 2))
            else:
                fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 6))
                
            
            plot_colormesh(ax1, dx, fig, X, Y, p, cmocean.cm.amp, colors.Normalize(vmin=27.4e3, vmax=27.6e3), r'$P (N/m)$', 'x (km)', 'y (km)')
            # plot_colormesh(ax1, dx, fig, X, Y, p_VP, cmocean.cm.amp, colors.Normalize(vmin=27.4e3, vmax=27.6e3), r'$P_{VP} (N/m)$', 'x (km)', 'y (km)')
            # plot_colormesh(ax2, dx, fig, X, Y, p_muphi, cmocean.cm.ice, colors.LogNorm(vmin=50e3, vmax=100e3), r'$P$ (N/m)', 'x (km)', None)
            # plot_colormesh(ax3, dx, fig, X, Y, np.abs(p_VP - p_muphi), cmocean.cm.ice, colors.Normalize(vmin=0, vmax=1e3), r'$|P_{VP} - P|$ (N/m)', 'x (km)', None)
            fig.tight_layout()
            
            plt.savefig(figdir+expno+'/icestrenght{}.png'.format(date))
            
            sigI = sigI[1:-1, 1:-1]
            sigII = sigII[1:-1, 1:-1]
            I = I[1:-1, 1:-1]
            fig = plt.figure()
            ax = plt.axes()
            if (np.shape(sigI)[0] == 1002):

            # plt.axis('equal')
            # sc = ax.scatter(sigI_norm.flatten(), sigII_norm.flatten(), s = 2, c = I.flatten(), cmap = plt.cm.magma, norm = colors.LogNorm(vmin=1e-5, vmax=1e-3))
                sc = ax.scatter(sigI[502:, :].flatten(), sigII[502:, :].flatten(), s = 2, c = 'r', cmap = plt.cm.magma, norm = colors.Normalize(vmin=1e-6, vmax=1e-3))
                sc = ax.scatter(sigI[:502, :].flatten(), sigII[:502, :].flatten(), s = 2, c = 'b', cmap = plt.cm.magma, norm = colors.Normalize(vmin=1e-6, vmax=1e-3))
                plt.plot(np.unique(sigI[502:, :]).ravel(), np.abs(np.unique(sigI[502:, :]).ravel()*mu_0))
                plt.plot(np.unique(sigI[502:, :]).ravel(), np.abs(np.unique(sigI[502:, :]).ravel()*mu_infty))
            
            else:
                sc = ax.scatter(sigI.flatten(), sigII.flatten(), s = 2, c = I.flatten(), cmap = plt.cm.magma, norm = colors.Normalize(vmin=1e-6, vmax=1e-3))
                # sc = ax.scatter(sigI.flatten(), sigII.flatten(), s = 2, c = 'b', cmap = plt.cm.magma, norm = colors.Normalize(vmin=1e-6, vmax=1e-3))
                plt.plot(np.unique(sigI).ravel(), np.abs((np.unique(sigI).ravel())*mu_0))
                plt.plot(np.unique(sigI).ravel(), np.abs((np.unique(sigI).ravel())*mu_infty)) 
                
                if arrows:
                    qpfac = 10
                    fac = 1
                    eu=np.hypot(divergence[::qpfac,::qpfac],fac*shear[::qpfac,::qpfac])
                    ax.quiver(sigI[::qpfac,::qpfac],fac*sigII[::qpfac,::qpfac],divergence[::qpfac,::qpfac]/eu,fac*shear[::qpfac,::qpfac]/eu, scale=10, color='red')
            # print(sigI.flatten())
            plt.grid('--')
            fig.colorbar(sc, label = 'I')
            ax.set_xlabel(r'$\sigma_{I}/P$')
            ax.set_ylabel(r'$\sigma_{II}/P$')
            plt.savefig(figdir+expno+'/stress_space{}.png'.format(date))
        

        fig, ((ax1, ax2)) = plt.subplots(2,1, figsize=(8, 6))

        plot_colormesh(ax1, dx, fig, X2, Y2, 1-A, cmocean.cm.ice, colors.LogNorm(vmin=1e-6, vmax=1e-1), r'$1-A$', None, None)
        plot_colormesh(ax2, dx, fig, X2, Y2, h-1, cmocean.cm.balance, colors.SymLogNorm(linthresh  = 1e-4, vmin=-1e-2, vmax=1e-2), r'$h-1$', None, None)

        plt.subplots_adjust(hspace = 0.1)
        fig.supxlabel('x (km)', fontsize = 16)
        fig.supylabel('y (km)', fontsize = 16)
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/tracers_{}.png'.format(date))
        
        
        # fig, ((ax1, ax2)) = plt.subplots(1,2, figsize=(8, 6))
        fig, ((ax1, ax2)) = plt.subplots(2,1, figsize=(8, 6))
        plot_colormesh(ax1, dx, fig, X2, Y2, divergence, cmocean.cm.balance, colors.SymLogNorm(1e-6, vmin=-1e-1, vmax=1e-1), 
            r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', None, None)
        plot_colormesh(ax2, dx, fig, X2, Y2, shear, cmocean.cm.amp, colors.LogNorm(vmin=1e-4, vmax=1e-2), 
            r'$\dot{\epsilon}_{\mathrm{II}} \text{ (day}^{-1})$', None, None) 
       
        plt.subplots_adjust(hspace = 0.1)
        fig.supxlabel('x (km)', fontsize = 16)
        fig.supylabel('y (km)', fontsize = 16)
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/shear_div_{}.png'.format(date))
        
        fig = plt.figure(figsize=(5, 3))
        ax1 = plt.axes()
        plot_colormesh(ax1, dx, fig, X2, Y2, divergence, cmocean.cm.balance, colors.SymLogNorm(1e-6, vmin=-1e-1, vmax=1e-1), 
            r'$\dot{\epsilon}_{\mathrm{I}} \text{ (day}^{-1})$', None, None)
       
        plt.subplots_adjust(hspace = 0.1)
        ax1.set_xlabel('x (km)', fontsize = 16)
        ax1.set_ylabel('y (km)', fontsize = 16)
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/div_{}.png'.format(date))

        
        plt.close()
            
def totdef_uniaxial(dates, expno, data_dict, dx, figdir, mu_0, mu_infty,angle_phi, MuPhi = True, log = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot,sigI_tot, sigII_tot, zeta_tot, eta_tot, uair_tot, vair_tot = \
            data_dict.values()
        
    Nx, Ny = np.shape(divergence_tot[0])
    
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        divergence, shear = divergence_tot[k], shear_tot[k]
        
        h = h_tot[k]

        tot_def = np.sqrt(divergence**2+shear**2)
        
        fig = plt.figure()
        ax = plt.axes()
        plot_colormesh(ax, dx, fig, X, Y, tot_def, cmocean.cm.balance, colors.LogNorm(vmin=1e-4, vmax=1e-0), 
        r'$\dot{\epsilon} \text{ (day}^{-1})$', 'x (km)', 'y (km)')
        
        plt.subplots_adjust(hspace = 0.1)
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/tot_def_{}.png'.format(date))
        
        fig = plt.figure()
        ax = plt.axes()
        plot_colormesh(ax, dx, fig, X, Y, divergence, cmocean.cm.balance, colors.SymLogNorm(linthresh=1e-3, vmin=-1e-2, vmax=1e-2), 
        r'$\dot{\epsilon}_\mathrm{I} \text{ (day}^{-1})$', 'x (km)', 'y (km)')
        
        plt.subplots_adjust(hspace = 0.1)
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/div_{}.png'.format(date))
        
        fig = plt.figure()
        ax = plt.axes()
        plot_colormesh(ax, dx, fig, X, Y, h-1, cmocean.cm.balance, colors.SymLogNorm(linthresh=1, linscale=1,
                                              vmin=-1e-3, vmax=1e-3, base=10), 
        r'$h-1$ (m)', 'x (km)', 'y (km)')
        
        # plot_colormesh(ax, dx, fig, X, Y, h-1, cmocean.cm.balance, colors.SymLogNorm(linthresh=0.5, linscale=1,
        #                                       base=10), 
        # r'$h-1$ (m)', 'x (km)', 'y (km)')
        
        plt.subplots_adjust(hspace = 0.1)
        # fig.tight_layout()
        plt.savefig(figdir+expno+'/h_{}.png'.format(date))
                 
def plot_single(Nx, Ny, dx, Var, label, figname):
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    x, y = np.meshgrid(Y, X)
    
    fig = plt.figure(figsize=(5, 4))
    ax1 = plt.axes()
    plot_colormesh(ax1, fig, Y, X, Var, cmocean.cm.ice, colors.LogNorm(vmin=1e-4, vmax=1e-2), label, 'x (km)', 'y (km)')
    plt.savefig(figname)

def plot_mean(Var, time, label, figname, expno, outputdir):
    
    plt.figure()
    plt.plot(time, Var, color = 'b')
    plt.xlabel('Time (min)')
    plt.ylabel(label)
    plt.grid()
    plt.savefig(outputdir+expno+'/'+figname)





def coordinates( x0, y0, RGPS = False):
    """
    Function that computes the latitude and longitude of the grid data using a moving cone inside (see book 1, page 147-148) and (see book 2, page 22 for RGPS).

    Args:
        x0 (np.ndarray): x distance in the tangent plane from the north pole
        y0 (np.ndarray): y distance in the tangent plane from the north pole
        
    This is taken from Antoine Savard SIM-Plots package. 

    Returns:
        np.ndarray: returns lon, lat in degrees
    """

    # convert to matrix
    x = np.broadcast_to(x0, (len(y0), len(x0)))
    y = np.broadcast_to(y0, (len(x0), len(y0))).T

    # polar coordinates on the plane
    r = np.sqrt((x) ** 2 + (y) ** 2)

    if not RGPS:
        lon = np.degrees(np.arctan2(y, x)) + BETA

        # angle of the cone
        tan_theta = r / (2 * R_EARTH)
        # short radius on sphere
        rs = 2 * R_EARTH * tan_theta / (1 + tan_theta ** 2)

        lat = np.degrees(np.arccos(rs / R_EARTH))

    elif RGPS:
        lon = np.degrees(np.arctan2(y, x)) + BETA_RGPS

        # small radius corresponding to plane at phi = 70
        rhat = R_EARTH * np.cos(np.pi / 2 - np.radians(PLANE_RGPS))
        lat = np.degrees(np.pi / 2 - np.arctan(r / rhat))

    return lon, lat

def pan_arctic_plot(ax,fig, data, bounds, lat, lon, fig_shape, colormap, datalabel):
    """
    This is a function to plot pan-arctic simulations taken from Antoine Savard SIM-Plots.  


    Args:
        data (_type_): _description_
        nx (_type_): _description_
        ny (_type_): _description_
        resolution (_type_): _description_
        fig_shape (_type_): _description_
    """
    
# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
    # if fig_shape == "round":
    #     theta = np.linspace(0, 2 * np.pi, 100)
    #     center, radius = [0.5, 0.5], 0.5
    #     verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    #     circle = mpath.Path(verts * radius + center)
    #     ax.set_boundary(circle, transform=ax.transAxes)

    # Limit the map to 65 degrees latitude and above.
    # ax.set_extent([-180, 180, 65, 90], crs=ccrs.PlateCarree())
    ax.set_aspect('equal', adjustable='box')
    # ax.add_feature(cfeature.OCEAN, color="black", zorder=0)
    cf = ax.pcolormesh(
        data,
        norm=colors.Normalize(vmin=bounds[0], vmax=bounds[1]),
        cmap=colormap)
    #     transform=ccrs.PlateCarree(),
    #     zorder=1,
    # )
    
    fig.colorbar(cf, ax = ax, label = datalabel, fraction=0.046, pad=0.04)
  
def plot_data_panarctic(figdir, datadict, dates, nx, ny, resolution, MuPhi = True):
    
    
    x0 = np.arange(nx + 2) * resolution - 2500
    y0 = np.arange(ny + 2) * resolution - 2250

    lon, lat = coordinates(x0, y0)
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()
            
    for k, date in enumerate(dates):
        
        divergence, shear = divergence_tot[k], shear_tot[k]
        h, A, p = h_tot[k], A_tot[k], p_tot[k]
        sigI, sigII = sigI_tot[k], sigII_tot[k]
        zeta, eta = zeta_tot[k], eta_tot[k]
        uair, vair = uair_tot[k], vair_tot[k]
        u, v = u_tot[k], v_tot[k]
        
        fig = plt.figure(figsize = (8, 5))
        ax1 = fig.add_subplot(1, 2, 1)#, projection = ccrs.NorthPolarStereo())
        ax2 = fig.add_subplot(1, 2, 2)#, projection = ccrs.NorthPolarStereo())

        
        pan_arctic_plot(ax1, fig, h, [0, 3], lat, lon, "round",cmocean.cm.ice, r'$h$ (m)')
        pan_arctic_plot(ax2, fig, A, [0.95, 1], lat, lon, "round",cmocean.cm.ice,r'$A$')
        plt.savefig(figdir+'thickness_concentration_{}.png'.format(date))
        
        
        fig = plt.figure(figsize = (8, 5))
        ax1 = fig.add_subplot(1, 2, 1)#, projection = ccrs.NorthPolarStereo())
        ax2 = fig.add_subplot(1, 2, 2)#, projection = ccrs.NorthPolarStereo())

        
        pan_arctic_plot(ax1, fig, divergence*60*60*24, [-5e-1, 1e-1], lat, lon, "round",cmocean.cm.amp, r'$\dot{\epsilon}_I$ (day$^{-1}$)')
        pan_arctic_plot(ax2, fig, shear*60*60*24, [0, 1e1], lat, lon, "round",cmocean.cm.amp, r'$\dot{\epsilon}_{II}$ (day$^{-1}$)')
        plt.savefig(figdir+'deformation_{}.png'.format(date))