#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean as cm
import matplotlib.lines as mlines


import scienceplots
plt.style.use('science')

#%%

def pressure_hibler(phi, h, Pmax = 27.5e3, C = 20):
    """
    Hibler pressure

    Args:
        phi (_type_): _description_
        h (_type_): _description_
        Pmax (_type_, optional): _description_. Defaults to 27.5e3.
        C (int, optional): _description_. Defaults to 20.

    Returns:
        _type_: _description_
    """
    
    return Pmax*h*np.exp(-C*(1-phi))


def I_shear(press0, shear, h, rhoi = 910,d = 1e3 ):
    
    return shear*d*np.sqrt(rhoi*h/press0+1e-20)

def pressure_friction(phi, h, tan_psi, Ishear, c_1 = 1e-2, c_2 = 1/4):
    """
    Friction pressure following Shi et al. 2021

    Args:
        phi (_type_): _description_
        h (_type_): _description_
        tan_psi (_type_): _description_
        Ishear (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    press0 = pressure_hibler(phi, h)
    
    press_fric = press0*np.exp(phi*c_1*Ishear**c_2*tan_psi)
    
    return press_fric

def press_collision(phi, h, shear, d = 1e3, rhoi = 910, Pmax = 27.5e3) : 
    
    """
    Collisional pressure

    Args:
        phi (_type_): _description_
        h (_type_): _description_
        shear (_type_): _description_
        d (_type_, optional): _description_. Defaults to 1e3.
        rhoi (int, optional): _description_. Defaults to 910.

    Returns:
        _type_: _description_
    """

    # P_max = pressure_hibler(phi, h)
    
    P_max = Pmax*h
    
    press_c = P_max * np.tanh(rhoi*h*(d*shear/(1-phi+1e-20))**2/(P_max+1e-20))
    # press_c[h == 0] = 0
    # 
    press_c_phi = P_max* np.tanh(rhoi*h*(d*shear/(1-phi+1e-20))**2/(P_max+1e-20))
    
    return press_c, press_c_phi


def dilatancy(mu, mu_0):
    
    return (mu-mu_0)/(1+mu*mu_0)


def friction_coefficient(phi, mu0, mu_infty, I_0 = 1e-3):
    
    delta_mu = mu_infty - mu0
    
    return mu0 + delta_mu/(I_0/(1-phi+1e-20) + 1)



#--- Constants ---# 
phi = np.linspace(0, 1, 100000)

shear_norm = 10**(np.arange(-5, -1, dtype = float))
shear = np.tile(shear_norm[:, None], (1, len(phi)))
    


h = 1

mu0 = 0.1
mu_infty = 0.8


#--- Compute pressure ---#


#-- Hibler --#
press_h = pressure_hibler(phi, h)


#-- Friction --#
mu = friction_coefficient(phi, mu0, mu_infty)
tan_psi = dilatancy(mu, 0.36)
Ishear = I_shear(press_h, shear, h)
press_f = pressure_friction(phi, h, tan_psi, Ishear, c_1 = 1)

#-- Collisional --#
press_c, press_c_phi = press_collision(phi, h, shear)


fig = plt.figure()
ax = plt.axes()

norm = mpl.colors.LogNorm(vmin=shear_norm.min(), vmax=shear_norm.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=plt.cm.coolwarm)
cmap.set_array([])
for i, norm in enumerate(shear_norm):
    f = plt.plot(1-phi, press_f[i],c=cmap.to_rgba(norm), linestyle = '--', label = r'$p_f$')
    c = plt.plot(1-phi, press_c[i],c=cmap.to_rgba(norm), linestyle = ':', label = r'$p_c$')
    # c = plt.plot(phi, press_c_phi[i]+press_f[i],c=cmap.to_rgba(norm), linestyle = '-.', label = r'$p_c$')
h = plt.plot(1-phi, press_h,color = 'k', linestyle = '-', label = r'$p_h$')

fig.colorbar(cmap, ax = ax, label = r'$\dot{\epsilon}_\mathrm{II}$ (s$^{-1}$)')


dashed_line = mlines.Line2D([], [], color='black', linestyle='--', label='$p^f$')
dot = mlines.Line2D([], [], color='black', linestyle=':', label='$p^c$')
full = mlines.Line2D([], [], color='black', linestyle='-', label='$p_H$')

plt.xlabel(r'$1-\phi$')
plt.ylabel(r'$P$ (N/m)')
secax = ax.secondary_xaxis('top')
secax.set_xlabel('$1-h$ (m)')
# plt.yscale('log')
plt.grid()
ax.legend([dashed_line, dot, full], [r'$p^f$', r'$p^c$', r'$p_H$'], loc='upper right')
plt.savefig('pressure.png')


fig = plt.figure()
ax = plt.axes()
i=2

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

f = plt.plot(1-phi, press_f[i]/1e3, linestyle = '--', label = r'$p_f$')
c = plt.plot(1-phi, press_c[i]/1e3, linestyle = ':', label = r'$p_c$')
h = plt.plot(1-phi, (press_f[i]+press_c[i])/1e3, color = 'r', linestyle = '-', label = r'$p_t$')
h = plt.plot(1-phi, (press_h)/1e3, color = 'k', linestyle = '-', label = r'$p_t$')

dashed_line = mlines.Line2D([], [], color='black', linestyle='--', label='$p^f$')
dot = mlines.Line2D([], [], color='black', linestyle=':', label='$p^c$')
full = mlines.Line2D([], [], color='r', linestyle='-', label='$p_t$')
# full = mlines.Line2D([], [], color='black', linestyle='.', label='$p_H$')

plt.xlabel(r'$1-\phi$')
plt.ylabel(r'$P$ (kN/m)')
# secax = ax.secondary_xaxis('top')
# secax.set_xlabel('$1-h$ (m)')
# plt.yscale('log')
# plt.grid()
ax.legend([dashed_line, dot, full], [r'$p^f$', r'$p^c$', r'$p_t$'], loc='upper right')
plt.savefig('pressure_simplified.png')





# %%
