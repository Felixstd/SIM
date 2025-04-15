import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.integrate as inte
import cmocean
import scienceplots

from random import randint
plt.style.use('science')

def wind_forcing(data_dict, N, dy, Ny, time, figdir, MuPhi = True):
    """
    
    The shape of the data is such that 
        shape(variable) = [Time, Y, X]

    Args:
        data_dict (_type_): _description_
        N (int): idx of the transect (in the y direction)
    """
    
    Y = np.arange(0, Ny*dy, dy)
    rhoair   =  1.3
    Cdair = 1.2e-3
    # print(Y)
    
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, uair_tot, vair_tot = \
            data_dict.values()
            
    uair_tot = np.array(uair_tot)
    vair_tot = np.array(vair_tot)


    uair_transect =  uair_tot[:, :, N]
    vair_transect =  vair_tot[:, :, N]
    
    wind_integral = np.trapz(abs(vair_transect) * vair_transect * rhoair*Cdair, x=Y, dx = dy, axis=1)
    
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (11, 5))
    
    ax1.scatter(time, np.mean(vair_transect, axis = 1), color = 'r', label = r'$v_{air}$')
    ax1.scatter(time, np.mean(uair_transect, axis = 1), color = 'k', label = r'$u_{air}$')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Velocity (m/s)')
    ax1.legend()
    ax1.grid()
    ax1.set_box_aspect(aspect=1)
    
    ax2.set_box_aspect(aspect=1)
    ax2.scatter(time, abs(wind_integral*Ny*dy))
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel( r'$L\int \tau_a dy$ (N)')
    ax2.grid()
    # ax1.legend()
    
    plt.savefig(figdir+'wind.png') 
    
def invariants(dates, expno, data_dict, nx, ny, deltax, maskC, figdir, mu_0, mu_infty):
        
    divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
        uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    
    
    for k, date in enumerate(dates):
        
        
        print('Plotting stresses: ', date)
        u, v = u_tot[k], v_tot[k]
        zeta, eta = zeta_tot[k], eta_tot[k]
        shear, mu = shearI_tot[k], muI_tot[k]
        P = p_tot[k]
        I = I_tot[k]
        
        
        sigI_norm, sigII_norm = compute_invariant_stresses(u, v, zeta, eta, shear, mu, maskC, P, nx, ny, deltax, mu0 = np.tan(20*np.pi/180))
        
        
        fig, (ax) = plt.subplots(1, 1)
        
        sc = ax.scatter(sigI_norm.flatten(), sigII_norm.flatten(), s = 2, c = I.flatten(), cmap = plt.cm.magma, norm = colors.Normalize(vmin=1e-6, vmax=1e-3))
        # sc = ax.scatter(sigI.flatten(), sigII.flatten(), s = 2, c = 'b', cmap = plt.cm.magma, norm = colors.Normalize(vmin=1e-6, vmax=1e-3))
        plt.plot(np.unique(sigI_norm).ravel(), np.abs(np.unique(sigI_norm).ravel()*mu_0))
        plt.plot(np.unique(sigI_norm).ravel(), np.abs(np.unique(sigI_norm).ravel()*mu_infty)) 
    # print(sigI.flatten())
        plt.grid('--')
        fig.colorbar(sc, label = 'I')
        # ax.set_xscale('symlog')
        # ax.set_yscale('symlog')
        ax.set_xlabel(r'$\sigma_{I}/P$')
        ax.set_ylabel(r'$\sigma_{II}/P$')
        plt.savefig(figdir+expno+'/stress_space_2_{}.png'.format(date))
        plt.close()
        
        # ax1.scatter(time, sigI_mean, color = 'r')
        # # ax1.scatter(time, sigII_mean, color = 'k', label = r'$\sigma_{II}/P$')
        # ax1.set_xlabel('Time (s)')
        # ax1.set_ylabel(r'$\overline{\sigma_{I}/P}$')
        # # ax1.legend()
        # ax1.grid()
        # ax1.set_box_aspect(aspect=1)
        
        # ax2.scatter(time, sigII_mean, color = 'k')
        # ax2.set_xlabel('Time (s)')
        # ax2.set_ylabel(r'$\overline{\sigma_{II}/P}$')
        # # ax2.legend()
        # ax2.grid()
        # ax1.set_xlabel('Time (s)')
        # ax2.set_box_aspect(aspect=1)
        
        plt.savefig(figdir+'invariant_time.png')
    
def dilatation(mu, mu0):
        
        tan_psi = (mu-mu0)/(1+mu*mu0)
        
        return tan_psi       

def compute_invariant_stresses(u, v, zeta, eta, shear, mu, maskC, P, nx, ny, deltax, mu0 = np.tan(20*np.pi/180)):
    
    sig11 = np.zeros((nx, ny))
    sig22 = np.zeros((nx, ny))
    sig12 = np.zeros((nx, ny))
    sig21 = np.zeros((nx, ny))
    sig1norm = np.zeros((nx, ny))
    sig2norm = np.zeros((nx, ny))
    sigInorm = np.zeros((nx, ny))
    sigIInorm = np.zeros((nx, ny))
    
    tan_psi = dilatation(mu, mu0)

    for j in range(0, ny-2):
        for i in range(0, nx-2):
            
            dudx = 0
            dudy = 0
            dvdx = 0
            dvdy = 0
            # print(i, j)
            #--- dudx ---#
            dudx = ( u[i+1,j] - u[i,j] ) / deltax
            
            #--- dvdy ---#
            dvdy = ( v[i,j+1] - v[i,j] ) / deltax
            
            #--- dvdx ---#
            if maskC[i + 1, j] + maskC[i - 1, j] == 2:
                    dvdx = ((v[i + 1, j] + v[i + 1, j + 1]) -
                            (v[i - 1, j] + v[i - 1, j + 1])) / (4.0 * deltax)
            
            elif maskC[i + 1, j] - maskC[i - 1, j] == 1:
                dvdx = (1.0 * (v[i + 1, j] + v[i + 1, j + 1]) +
                        3.0 * (v[i, j] + v[i, j + 1])) / (6.0 * deltax)
            
            elif maskC[i + 1, j] - maskC[i - 1, j] == -1:
                dvdx = (-1.0 * (v[i - 1, j] + v[i - 1, j + 1]) -
                        3.0 * (v[i, j] + v[i, j + 1])) / (6.0 * deltax)
            
            #--- dudy ---#
            if maskC[i, j + 1] + maskC[i, j - 1] == 2:
                dudy = ((u[i, j + 1] + u[i + 1, j + 1]) -
                            (u[i, j - 1] + u[i + 1, j - 1])) / (4.0 * deltax)
            
            elif maskC[i, j + 1] - maskC[i, j - 1] == 1:
                dudy = (1.0 * (u[i, j + 1] + u[i + 1, j + 1]) +
                        3.0 * (u[i, j] + u[i + 1, j])) / (6.0 * deltax)
            
            elif maskC[i, j + 1] - maskC[i, j - 1] == -1:
                dudy = (-1.0 * (u[i, j - 1] + u[i + 1, j - 1]) -
                        3.0 * (u[i, j] + u[i + 1, j])) / (6.0 * deltax)
            
            
            ep = dudx + dvdy
            
            e12 = 1/2*(dudy + dvdx)
            e21 = e12
            
            sig11[i, j] = zeta[i, j]*ep+eta[i, j]*dvdx-(P[i, j]+zeta[i, j]*shear[i, j]*tan_psi[i,j])
            sig22[i, j] = zeta[i, j]*ep+eta[i, j]*dvdy-(P[i, j]+zeta[i, j]*shear[i, j]*tan_psi[i,j])
            sig12[i, j] = 2*e21*eta[i, j]
            sig21[i, j] = 2*e12*eta[i, j]
            
            sigp = sig11[i, j]+sig22[i, j]
            sigm=sig11[i, j]-sig22[i, j]

            sigTmp=np.sqrt(np.abs(sigm**2+4*sig12[i, j]*sig21[i, j]))
            
            
            sig1norm[i, j]=0.5*(sigp+sigTmp)/P[i, j]
            sig2norm[i, j]=0.5*(sigp-sigTmp)/P[i, j]
            
            sigInorm[i,j] = 0.5*(sig1norm[i, j]+sig2norm[i, j])
            sigIInorm[i,j] = 0.5*(sig1norm[i, j]-sig2norm[i, j])
            
    return sigInorm, sigIInorm
    
               
def mean_values(dates, data_dict, dx, dt, MuPhi = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot,sigI_tot, sigII_tot, zeta_tot, eta_tot, uair_tot, vair_tot = \
            data_dict.values()
        
    Nx, Ny = np.shape(divergence_tot[0])
    
    X = np.arange(0, Nx)*dx/1e3
    Y = np.arange(0, Ny)*dx/1e3
    
    mean_shear_time = np.zeros(len(dates))
    mean_div_time   = np.zeros(len(dates))
    mean_mu_time    = np.zeros(len(dates))
    mean_sigI_time  = np.zeros(len(dates))
    mean_sigII_time = np.zeros(len(dates))
    
    for k, date in enumerate(dates):
        
        print('Analysing mean variables: ', date)
        
        mu_t    = muI_tot[k]
        div_t   = divergence_tot[k]
        shear_t = shear_tot[k]
        
        div_t[div_t < 100] = np.nan
        shear_t[shear_t < 0 ] = np.nan
        
        
        mean_mu_time[k]    = np.mean(mu_t)
        mean_div_time[k]   = np.nanmean(div_t)
        mean_shear_time[k] = np.nanmean(shear_t)
        
        
        
    return mean_mu_time, mean_div_time, mean_shear_time
    
def velocity_transects(data_dict, dates, dx,figdir, expno, MuPhi = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot,sigI_tot, sigII_tot, zeta_tot, eta_tot, uair_tot, vair_tot = \
            data_dict.values()
    
    Nx2, Ny2 = np.shape(u_tot[0])
    
    X2, Y2 =  np.arange(0, Nx2)*dx/1e3,  np.arange(0, Ny2-1)*dx/1e3

    
    Half_domain = Nx2//2
    
    
    colors = ['0C5DA5', '00B945', 'FF9500', 'FF2C00', '845B97', 'forestgreen', 'orangered']
    fig, ((ax1, ax2), (ax3, ax4)) =  plt.subplots(2, 2, sharey = True, figsize = (5, 4))
    axs = [ax1, ax2, ax3, ax4]
    for ax in axs:
        ax.grid()

    date_list = []
    lines = []
    for k, date in enumerate(dates):
        if k % 10 == 0 :
            u, v       = u_tot[k], v_tot[k]
            uair, vair = uair_tot[k],  vair_tot[k]
            
            uair = uair[:, :-1]
            u = u[:, :-1]
            vair = vair[:-1, :]
            v = v[:-1, :]
            
            color = '#%06X' % randint(0, 0xFFFFFF)
            line = ax1.plot(u[:, Half_domain], X2, color = color, linestyle = '-',)
            ax3.plot(uair[:, Half_domain], X2, color = color, linestyle = '--')
            
            ax2.plot(v[:, Half_domain], X2, color = color, linestyle = '-', label = 'v') 
            ax4.plot(vair[:, Half_domain], X2, color = color, linestyle = '--')
            date_list.append(date)
            lines.append(line[0])
            
    ax3.set_xlabel(r'$u$ (m/s)')
    ax4.set_xlabel(r'$v$ (m/s)')
    plt.suptitle('Sea Ice')
    
    ax1.set_ylabel(r'$y$ (km)')
    ax3.set_ylabel(r'$y$ (km)')
    
    plt.subplots_adjust(hspace=0.5)

# Add text in figure coordinates
    plt.figtext(0.5, 0.5, 'Wind', ha='center', va='center')

    fig.legend(lines, date_list, loc = 'outside center right', bbox_to_anchor = (1.3, 0.5))
    plt.savefig(figdir+expno+'/velocitytransect.png')
    
    
        
    
       
    
    
    
    
    
     
            
