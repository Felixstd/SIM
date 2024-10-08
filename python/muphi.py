import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

def mu(A, I_0, mu_0, mu_infty, A_max = 1):

    delta_mu = mu_infty - mu_0
    dilantancy = A_max-A

        
    friction = mu_0 + delta_mu /(I_0 / dilantancy + 1)
    friction[dilantancy < 1e-12] = mu_0
    
    return friction

def p(A, h, d_average, e_II, rho_i = 900):

    return rho_i * h *(d_average * e_II /(1 - A)) **2

def zeta(p, mu, e_II, mu_b ): 
        
    zeta = np.minimum((mu_b + mu/2) * p/e_II, 1e12)
    
    return zeta

def eta(p, mu, e_II): 
    
    eta = np.minimum((mu/2) * p/e_II, 2e8*p)
    
    return eta


def plotting(h, A, mu, p, zeta, eta, x, I_0, d_average,mu_0, mu_infty, mu_b, filename):

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize =(7,7), sharex = True)


    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    
    for ax in axes:
        ax.grid()
        ax.set_box_aspect(aspect=1)
    # ax1.set_aspect('equal', adjustable='box')
    ax1.plot(x, h, color = 'g')
    ax1.set_ylabel('Ice thickness (m)')

    ax2.plot(x, A, color = 'g')
    ax2.set_ylabel('Ice Concentration')


    ax3.plot(x, mu, color = 'g')
    ax3.set_ylim(mu_0, mu_infty)
    ax3.set_yticks(np.arange(mu_0, mu_infty+0.2, 0.1))
    ax3.set_ylabel(r'$\mu$')
    ax3.set_xlabel('x (m)')

    ax4.plot(x, p, color = 'g')
    ax4.set_ylabel(r'$p$ (N/m)')
    ax4.set_xlabel('x (m)')
    
    ax5.plot(x, zeta, color = 'g')
    ax5.set_ylabel(r'$\zeta$ (N/m)')
    ax5.set_xlabel('x (m)')

    ax6.plot(x, eta, color = 'g')
    ax6.set_ylabel(r'$\eta$ (N/m)')
    ax6.set_xlabel('x (m)')
    
    fig.suptitle(r'$\mu-\Phi$({},{},{},{}, {})'.format(I_0, d_average, mu_0, mu_infty, mu_b))

    plt.savefig(filename)



x = np.arange(0, 20000, 1)
L = 2000
I_0 = 1e-3
mu_0 = 0.1
mu_infty = 0.9
mu_b = 1
d_average = 1e3
e_II = 1e-2*1/(3600*24)

I_0_tot = np.linspace(1e-5, 1, 10)

h = (1-np.tanh((x-10000)/L)) /2
A = h

friction_mu = mu(A, I_0, mu_0, mu_infty)
pressure = p(A, h, d_average, e_II)
eta_x = eta(pressure, friction_mu, e_II)
zeta_x = zeta(pressure, friction_mu, e_II, mu_b)

plotting(h, A, friction_mu, pressure, eta_x, zeta_x, x, I_0, d_average, mu_0, mu_infty,mu_b, 'idealized_var.png')





