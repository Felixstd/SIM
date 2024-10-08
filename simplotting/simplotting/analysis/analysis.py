import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as inte

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
    
def invariants(data_dict, N, time, figdir, MuPhi = True):
    
    if MuPhi:
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, \
            uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = data_dict.values()
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, uair_tot, vair_tot = \
            data_dict.values()
    
    sigI_tot_norm  = np.array(sigI_tot)/np.array(p_tot)
    sigII_tot_norm = np.array(sigII_tot)/np.array(p_tot)
    
    
    sigI_transect  =  sigI_tot_norm[:, :, N]
    sigII_transect =  sigII_tot_norm[:, :, N]
    
    sigI_mean = np.nanmean(sigI_transect, axis = 1)
    sigII_mean = np.nanmean(sigII_transect, axis = 1)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (11, 5))
    
    ax1.scatter(time, sigI_mean, color = 'r')
    # ax1.scatter(time, sigII_mean, color = 'k', label = r'$\sigma_{II}/P$')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel(r'$\overline{\sigma_{I}/P}$')
    # ax1.legend()
    ax1.grid()
    ax1.set_box_aspect(aspect=1)
    
    ax2.scatter(time, sigII_mean, color = 'k')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel(r'$\overline{\sigma_{II}/P}$')
    # ax2.legend()
    ax2.grid()
    ax1.set_xlabel('Time (s)')
    ax2.set_box_aspect(aspect=1)
    
    plt.savefig(figdir+'invariant_time.png')
            
            
            
            
