"""
Code for figure 03 in St-Denis et al. 2025.
This code plots the friction coefficient in terms of I
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
plt.rcParams.update({'font.size': 15})
plt.rc('axes', labelsize=15) 

#%%
def friction_coefficient_jop(I, mu_infty, mu_0, I_0): 
    """
    This function computes the friction coefficient as 
    defined in (MIDI 2004, DaCruz 2005, and Jop. 2006)

    Args:
        I (array): inertial numbers
        mu_infty (float): dynamic friction coefficient
        mu_0 (float): static friction coefficient
        I_0 (float): regime frontiers 

    Returns:
        mu (array): friction coefficient for a range of 
        I. 
    """
    
    #compute range of mu
    delta_mu = mu_infty - mu_0
    
    #compute mu
    mu = mu_0 + delta_mu/(I_0/(I+1e-20) + 1)
    
    return mu

def friction_coefficient_cohesion(I, mu_infty, mu_0, I_0, C, beta, alpha, K=1e-3, b = 1/2, I_1 = 1e-2):
        
        mu_jop = friction_coefficient_jop(I, mu_infty, mu_0, I_0)
        mu_c   = 1.31*C/(1-beta*np.log(1-I/(1+alpha*C)**(1/2)))
        
        
        W = C*K**(1/2)/I
        g = 1+  b*W/(1+I_1/I)

        return mu_jop*g

def dilatation(mu, mu0):
        
        tan_psi = (mu-mu0)/(1+mu*mu0)
        
        return tan_psi

def concentration(I):
    
    return np.maximum(1-I, 0)


#%%
#--- Constants ---# 
I = np.linspace(0, 1, 1000000)
mu_0 = 0.1
mu_infty = 0.9
I_0 = 1e-3
C = 5e-1
beta = 1
alpha = 1

#--- Main Computations ---#

#friction coefficient
mu_jop = friction_coefficient_jop(I, mu_infty, mu_0, I_0)

#friction coefficient at I_0
mu_j_i_0 = friction_coefficient_jop(I_0, mu_0, mu_infty, I_0)

mu_cohesive = friction_coefficient_cohesion(I, mu_infty, mu_0, I_0, C,alpha, beta )

mu_0_dilat = np.tan(10*np.pi/180)
# tan_psi = dilatation(mu_jop, mu_0_dilat)


tan_psi = dilatation(mu_jop, mu_cohesive)

phi = concentration(I)

#%%

#--- Figure ---#

# --  Code to make figure 3 of the paper -- #

plt.figure(figsize = (5.5, 4))
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.axvline(I_0, ymax = mu_j_i_0,
            color = 'r')
plt.axvline(5e-5,linestyle = '--',
            color = 'k')
plt.axvline(5e-2,linestyle = '--',
            color = 'k')

plt.plot(I, mu_jop, color = 'k')
plt.scatter(I_0, mu_j_i_0, color = 'r', zorder = 2)


ticks_y = [mu_0, mu_infty]
ax.set_yticks(ticks_y)
ax.set_yticklabels([r'$\mu_0$', r'$\mu_\infty$'])

ax.text(0.03, 1.05, 'Quasi-Static', transform=ax.transAxes, fontsize=15,
        verticalalignment='top')

ax.text(0.37, 1.05, 'Dense-Inertial', transform=ax.transAxes, fontsize=15,
        verticalalignment='top')

ax.text(0.8, 1.05, 'Collisional', transform=ax.transAxes, fontsize=15,
        verticalalignment='top')

ax.text(0.52, 0.5, r'$(I_0, \mu(I_0))$', transform=ax.transAxes, fontsize=15,
        verticalalignment='top')


plt.xscale('log')
plt.xlabel('$I$')
plt.ylabel(r'$\mu$')
plt.savefig('mu.png')#, dpi = 500, bbox_inches = 'tight')


plt.figure()
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.plot(I, phi, color = 'k')
# plt.grid()
plt.xlabel(r'I')
plt.ylabel(r'$\phi$')
plt.savefig('conc.png')


# plt.figure()
# ax = plt.axes()
# plt.plot(mu_jop, tan_psi, color = 'r')
# plt.grid()
# plt.axvline(mu_0_dilat,
#             color = 'k', label = r'$\mu_0$')
# plt.legend()
# plt.xlabel(r'$\mu$')
# plt.ylabel(r'$\tan \psi$')
# plt.savefig('dilatation.png')


#--- Comparaison between cohesion or not ---# 

# plt.figure()
# plt.plot(I, mu_jop, label = r'$\mu_j$')
# plt.plot(I, mu_cohesive, label = r'$\mu_c$')
# plt.xlabel('I')
# plt.ylabel(r'$\mu$')
# plt.legend()
# plt.grid()
# plt.xscale('log')
# plt.savefig('mu_cohesive.png')


# %%
