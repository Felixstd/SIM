import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

def wavespeed_mu(A_0, mu_0, mu_infty, P_star, C, I_0, mu_b, rhoi = 900):
    
    Gamma = mu_b + mu_0 + (mu_infty - mu_0)/(I_0/(1-A_0)+1)
    P_0 = P_star*np.exp(-C*(1-A_0))
    gamma = (mu_infty - mu_0)*I_0/((I_0+1-A_0)**2)
    
    cp_plus = np.sqrt((1/rhoi)*(Gamma *P_0 * (1+C*A_0) - P_0 * gamma * A_0))

    return cp_plus

def wavespeed_vp(A_0, C, P_star, e, rho_i = 900):
    
    lamda = (-np.sqrt(1+e**(-2))-1)/2
    lamda_star = lamda * P_star*np.exp(-C*(1-A_0))
    cp = np.sqrt(-lamda_star / rho_i * (1+C*A_0))
    
    return cp


A_0 = np.linspace(0.5, 1, 1000)
mu_0 = 0.4
mu_infty = 0.8
I_0 = 1e-3
P_star = 27.5e3
C = 20
mu_b = 1
e=2

deltat = np.linspace(0, 300, 1000)
deltax = np.linspace(0, 10, 1000)*1e3

cp = wavespeed_mu(A_0, mu_0, mu_infty, P_star,C, I_0, mu_b)
cp_vp = wavespeed_vp(A_0, C, P_star, e )

# CFL = cp * deltat

cp_typical = wavespeed_mu(0.99, mu_0, mu_infty, P_star,C, I_0, mu_b)
CFL = deltax/cp_typical
# 
print('Typical cp: ', cp_typical, CFL[-1])


plt.figure()
plt.plot(A_0, cp, label = 'mu')
plt.plot(A_0, cp_vp, label = 'vp')
plt.grid()
plt.legend()
plt.xlabel(r'$A_0$')
plt.ylabel(r'$\omega / k$ (m/s)')
plt.savefig('wavespeed_plastic.png')


plt.figure()
plt.plot(deltax/1e3, CFL)
plt.fill_between(deltax/1e3, CFL)
plt.grid()
plt.xlabel(r'$\Delta x$ (km)')
plt.ylabel(r'$\Delta t$ (s)')
plt.title(r'For $A_0 = 0.99$')
plt.savefig('CFL.png')
    
    
    
    
    
    
    
    