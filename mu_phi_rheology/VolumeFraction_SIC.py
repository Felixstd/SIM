import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

phi_min = 0.0
phi_max = 0.9

A = np.arange(0,1+1e-4,1e-4)

def phiofA(A,phi_min,phi_max,  tau = 0.2):
    
    phi_s = (1-phi_min) * A + phi_min
    
    phi_g = - (1-phi_max) * np.tanh((A-1)/tau) + (phi_max)
     
    phi = phi_s * phi_g
    
    return phi, phi_s, phi_g

phi, phi_s, phi_g = phiofA(A,phi_min,phi_max, tau = 0.5)

plt.figure(figsize=(4,4))
plt.plot(A,phi,label='corrected')
plt.plot(A,phi_s, ls='dotted', label='linear')
plt.plot(A,phi_g,label='correction')
plt.plot()
plt.axhline(phi_max,color='k',ls='dashed')

plt.axhline(0,color='k')
plt.axvline(1,color='k')
plt.axvline(0,color='k')
plt.axhline(1,color='k')
plt.grid()
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0.5))
plt.axis('equal')
plt.ylim([0-0.1,1+0.1])
plt.xlim([0-0.1,1+0.1])
plt.xlabel('$A$')
plt.ylabel(r'$\phi$')


plt.savefig('phi_SIC_correspondance.png')