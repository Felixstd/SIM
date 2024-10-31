#%%
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

#%%
def friction_coefficient_jop(A, mu_infty, mu_0, I_0): 
    
    delta_mu = mu_infty - mu_0
    I = 1-A
    
    mu = mu_0 + delta_mu/(I_0/I + 1)
    
    return mu

def friction_coefficient_barker(A, mu_infty, mu_0, I_0, I_shift, alpha = 1.99): 
    
    I = 1-A
    
    A_minus = I_shift * np.exp(alpha/friction_coefficient_jop(1-I_shift, mu_infty, mu_0, I_0))
    
    mu = np.zeros_like(I)
    
    mu[I<=I_shift] = np.sqrt(alpha/np.log(A_minus/I[I<=I_shift]))
    mu[I>I_shift] = friction_coefficient_jop(1-I[I>I_shift], mu_infty, mu_0, I_0)
    
    
    return mu

#%%
A = np.linspace(0.97, 1, 10000)
mu_0 = 0.2
mu_infty = 0.9
I_0 = 0.2
I_shift = 0.01



mu_jop = friction_coefficient_jop(A, mu_infty, mu_0, I_0)
mu_barker = friction_coefficient_barker(A, mu_infty, mu_0, I_0, I_shift)

#%%
plt.figure()
plt.plot(1-A, mu_jop, label = r'$\mu_J$')
plt.plot(1-A, mu_barker, label = r'$\mu_B$')
plt.grid()
plt.xlabel('1-A')
plt.ylabel(r'$\mu$')
plt.legend()
plt.savefig('mu.png')

# %%
