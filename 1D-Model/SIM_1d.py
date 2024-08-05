#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

from numba import njit

#%%
def first_derive_taylor_h2(variable, deltax ) : 
    """
    Function to compute the derivative in x

    Args:
        variable (_type_): _description_
        a (_type_): _description_
        b (_type_): _description_
        n (_type_): _description_

    Returns:
        _type_: _description_
    """

    derive = np.zeros_like(variable)
    n = len(derive)-1


    for i in range(0, n+1) : 

        #Non-centered formula for the first and last point
        if i == 0 : 
            derive[i] = (-3*variable[i] + 4*variable[i+1] - variable[i+2])/(2*deltax)

        elif i == n : 
            derive[i] = (3*variable[i] - 4*variable[i-1] + variable[i-2])/(2*deltax)
            

        else : 
            derive[i] = (variable[i+1] - variable[i-1])/(2*deltax)


    return derive

@njit()
def inertial_number(shear, p) :

    
    
    
    I = np.zeros_like(shear)
    P_min = 1e-6
    
    for i in range(len(p)): 
        
        if p[i] < P_min:
            I[i] = 1e6
        else:
            I[i] = d_mean * abs(shear[i]) * np.sqrt(H * rhoice / p[i])

    return I

def dilatancy(inertial_number):
    
    return Phi_0 - c_phi*inertial_number

@njit()
def pressure(shear, phi):
    p = np.zeros_like(shear)
    
    # p_max = 27e3
    
    p_max = 27e3*H
    # p_max  = 2.7 
    
    for i in range(len(shear)): 
        
        if (phi[i] - Phi_0) < 1e-7:
            p[i] = p_max
        else:
            p[i] = np.nanmin([p_max, rhoice * H * (shear[i] / (phi[i] - Phi_0))**2 * d_mean**2])
    return p

@njit()
def mu(I): 
    
    min_I = 1e-6
    mu_array = np.zeros_like(I)
    
    for i in range(len(I)):
        if I[i] < min_I:
            mu_array[i] = mu_infty
            
        else:
            mu_array[i] = mu_0 + (mu_infty - mu_0)/(I_0/I[i] + 1)

    
    return mu_array
    
@njit()
def zeta_muphi(mu_I, mu_b, shear, P):
    
    zeta = np.zeros(shear.shape[0]-2)
    
    P = P*H

    for i in range(1,len(zeta)-1): 
        
        if shear[i] < 1e-8: 
            zeta[i] = 2.5e8*P[i]
        
        else:
            zeta[i] = np.nanmin([(mu_b+mu_I[i]/2)*P[i]/shear[i], 2.5e8*P[i]])
        
    return zeta 


def zeta_vp(shear, P = 27e3): 
    
    P = P*H
    zeta = np.zeros(shear.shape[0]-2)
    
    for i in range(1,len(zeta)-1): 
        
        if shear[i] < 1e-8: 
            zeta[i] = 2.5e8*P
        
        else:
            delta = (shear[i]**2* (1 + 1/e_square))
            zeta[i] = np.nanmin([P/(delta*2), 2.5e8*P])
        
    return zeta



#%%
#---------- Constants ---------#
L = 50e3
Nx = 200
delta_x = L /(Nx -1)
delta_x_2 = delta_x**2

Nt = 100
# delta_t = 3
# T_tot = delta_t*Nt
T_tot = 50000

delta_t = T_tot/Nt

delta_t = 1



H = 0.1
rhoice = 900
d_mean = delta_x/2
Phi_0 = 1
c_phi = 0.5
mu_0 = 0.1
mu_infty = 0.43
I_0 = 6.7e-3
tau = 0.09
mu_b = 0.3

x = np.arange(0, L + delta_x, delta_x)
t_tot = np.arange(0, T_tot, delta_t)


#---------- VP -----------#

#%%
zeta = 1e5
e = 2
e_square = e**2

u_vp = np.zeros((Nx, Nt+1))
shear_vp = np.zeros((Nx, Nt+1))

#initial conditions
u_vp[:, 0] = np.zeros_like(x)
shear_vp[:, 0] = np.zeros_like(x)


#Matrix A
# alpha = zeta * delta_t / delta_x_2**2
# beta = rhoice * H + 2 * alpha

# A = np.zeros((3, Nx-2))
# A[0, 1:] = -alpha  # Upper diagonal
# A[1, :] = beta     # Main diagonal
# A[2, :-1] = -alpha # Lower diagonal

# Time-stepping loop
for n in range(Nt):
    
    zeta = zeta_vp(shear_vp[:, n])
    
    alpha = zeta * delta_t / delta_x_2**2
    beta = rhoice * H + 2 * alpha

    A = np.zeros((3, Nx-2))
    A[0, 1:] = -alpha[1:]  # Upper diagonal
    A[1, :] = beta     # Main diagonal
    A[2, :-1] = -alpha[:-1] # Lower diagonal


    b = u_vp[1:-1, n] * rhoice*H + delta_t * tau
    u_vp[1:-1, n + 1] = solve_banded((1, 1), A, b)
    # Boundary conditions
    u_vp[0, n + 1] = 0
    u_vp[-1, n + 1] = 0
    
    shear_vp[:, n+1] = first_derive_taylor_h2(u_vp[:, n+1], delta_x)
    
    

print(u_vp)
plt.figure(figsize=(8, 6))
for n in range(0, Nt + 1, Nt//10):
    # print(n)
    u_vp[1, n] = 0
    u_vp[-2, n] = 0
    plt.plot(x, u_vp[:, n], label=f't = {(n * delta_t):.2f} s')
plt.xlabel('x (m)')
plt.ylabel('u (m/s)')
plt.title('Vicous-Plastic')
plt.legend()
plt.savefig('u_vp.png', dpi = 500, bbox_inches = 'tight')




# %%

#%%
#-------- mu - Phi -----------#


u_muphi = np.zeros((Nx, Nt+1))
I_muphi = np.zeros((Nx, Nt+1))
p_muphi = np.zeros((Nx, Nt+1))
phi_muphi = np.zeros((Nx, Nt+1))
mu_muphi = np.zeros((Nx, Nt+1))
shear_muphi = np.zeros((Nx, Nt+1))

# u_muphi[:, 0] = u_vp[:, -1]

u_muphi[:, 0] = np.ones_like(x)*1e-8


p_init = np.ones_like(x)*1e-8


p_muphi[:, 0] = p_init
shear_muphi[:, 0] = first_derive_taylor_h2(u_muphi[:, 0], delta_x)
I_muphi[:, 0] = inertial_number(shear_muphi[:, 0], p_muphi[:, 0])
mu_muphi[:, 0] = mu(I_muphi[:, 0])
phi_muphi[:, 0] = dilatancy(I_muphi[:, 0])


i = 0
for n in range(Nt):
    
    p_muphi[:, n] = pressure(shear_muphi[:, n], phi_muphi[:, n])
    
    
    zeta = zeta_muphi(mu_muphi[:, n], mu_b, shear_muphi[:, n], p_muphi[:, n])
    
    alpha = zeta * delta_t / delta_x_2**2
    beta = rhoice * H + 2 * alpha

    A = np.zeros((3, Nx-2))
    A[0, 1:] = -alpha[1:]  # Upper diagonal
    A[1, :] = beta     # Main diagonal
    A[2, :-1] = -alpha[:-1] # Lower diagonal
    
    
    b = u_muphi[1:-1, n] * rhoice*H + delta_t * tau
    u_muphi[1:-1, n + 1] = solve_banded((1, 1), A, b)
    # Boundary conditions
    u_muphi[0, n + 1] = 0
    u_muphi[-1, n + 1] = 0
    
    print('iteration : ', i)

    
    shear_muphi[:, n+1] = first_derive_taylor_h2(u_muphi[:, n+1], delta_x)
    I_muphi[:, n+1] = inertial_number(shear_muphi[:, n+1], p_muphi[:, n])
    phi_muphi[:, n+1] = dilatancy(I_muphi[:, n+1])
    mu_muphi[:, n+1] = mu(I_muphi[:, n+1])
    
    
    i+=1
    
    
plt.figure(figsize=(8, 6))
for n in range(0, Nt+1 , Nt//10):
    plt.plot(x, u_muphi[:, n], label=f't = {(n * delta_t):.2f}s')
plt.xlabel('x (m)')
plt.ylabel('u (m/s)')
plt.title(r'$\mu(I) - \Phi(I)$')
plt.legend()

plt.savefig('u_muphi.png', dpi = 500, bbox_inches = 'tight')

plt.figure(figsize=(8, 6))
for n in range(0, Nt + 1, Nt//10):
    plt.plot(x, I_muphi[:, n], label=f't = {(n * delta_t):.2f} s')
plt.xlabel('x')
plt.ylabel('I')
plt.legend()

plt.savefig('I_muphi.png', dpi = 500, bbox_inches = 'tight')

#---- Comparaison between the 2 ----# 

plt.figure()
plt.plot(x, u_vp[:, -1], color = 'r', label = 'VP')
plt.plot(x, u_muphi[:, -1], color = 'k', label = r'$\mu(I) - \Phi(I)$', alpha = 0.5)
plt.xlabel('x (m)')
plt.ylabel('u (m/s)')
plt.legend()
plt.grid()
plt.savefig('comp_mu_VP.png', dpi = 500, bbox_inches = 'tight')

    

    
    











    
    
    


# %%
