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
            

def compute_invariant_stresses(u, v, zeta, eta, maskC, P, nx, ny, deltax):
    
    print(np.shape(eta))
    sig11 = np.zeros((nx, ny))
    sig22 = np.zeros((nx, ny))
    sig12 = np.zeros((nx, ny))
    sig21 = np.zeros((nx, ny))
    sig1norm = np.zeros((nx, ny))
    sig2norm = np.zeros((nx, ny))
    sigI = np.zeros((nx, ny))
    sigII = np.zeros((nx, ny))
    sigInorm = np.zeros((nx, ny))
    sigIInorm = np.zeros((nx, ny))

    print(ny-3)
    for j in range(0, ny-2):
        for i in range(0, nx-1):
            
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
            em = dudx - dvdy
            
            e12 = 1/2*(dudy + dvdx)
            e21 = e12
            
            sig11[i, j] = zeta[i, j]*ep+eta[i, j]*em-P[i, j]
            sig22[i, j] = zeta[i, j]*ep-eta[i, j]*em-P[i, j]
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
    
               
            
            
