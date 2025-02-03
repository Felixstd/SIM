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

def press_collision(phi, h, shear, d = 1e4, rhoi = 910, Pmax = 27.5e3) : 
    
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

    P_max = pressure_hibler(phi, h)
    
    # P_max = Pmax*h
    
    press_c = P_max * np.tanh(rhoi*h*(d*shear/(1-phi+1e-20))**2/(P_max+1e-20))
    
    # press_c= np.minimum(5e3, rhoi*h*(d*shear/(1-phi+1e-20))**2)
    # press_c[h == 0] = 0
    # 
    press_c_phi = P_max* np.tanh(rhoi*h*(d*shear/(1-phi+1e-20))**2/(P_max+1e-20))
    
    return press_c, press_c_phi


def dilatancy(mu, mu_0):
    
    return (mu-mu_0)/(1+mu*mu_0)


def friction_coefficient(phi, mu0, mu_infty, I_0 = 1e-3):
    
    delta_mu = mu_infty - mu0
    
    return mu0 + delta_mu/(I_0/(1-phi+1e-20) + 1)

# def pressure_sutherland(A, h, dilat, gamma = 1/2, Pstar = 27.5e3, C = 20):

    
#     return Pstar*h*np.tanh(dilat*h)*np.exp(-C*(1-A))

def pressure_sutherland(A, h, dilat, gamma = 5, Pstar = 27.5e3, C = 20):

    # return Pstar*h*np.exp(-(1-dilat))
    return Pstar*h*np.exp(-(C)*(1-A))*np.exp(dilat)

def press_rothrcok(A, h, mu, C = 20, g = 9.80, rhow = 1026, rhoi = 900, tanphi = 0.8, k = 5):
    
    
    cp = 1/2*g*(rhow - rhoi)*rhoi/(rhow)
    print('cp', cp)
    
    cf = mu*(rhow - rhoi)*g*(rhoi*(k-1)/rhow)**2/(2*tanphi)
    print((rhoi*(k-1)/rhow)**2)
    
    print('cf', cf)
    
    p = (k*cp + k/(k-1)*cf)*h**2/6
    
    return p
    
#--- Constants ---# 
phi = np.linspace(0, 1, 100000)

shear_norm = 10**(np.arange(-5, -1, dtype = float))
shear = np.tile(shear_norm[:, None], (1, len(phi)))
    


h = 1

mu0 = 0.1
mu_infty = 0.9


#--- Compute pressure ---#


#-- Hibler --#
press_h = pressure_hibler(phi, h)



#-- Friction --#
mu = friction_coefficient(phi, mu0, mu_infty, I_0=1e-3)
tan_psi = dilatancy(mu, np.tan(20*np.pi/180))
Ishear = I_shear(press_h, shear, h)
press_f = pressure_friction(phi, h, tan_psi, Ishear, c_1 = 1)

press_suth = pressure_sutherland(phi, h, tan_psi , Pstar= 27.5e3)
press_roth = press_rothrcok(phi, 3, mu, g = 9.80, rhow = 1026, rhoi = 900, tanphi = 0.8, k = 5)
print(tan_psi, press_suth)
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
i=1

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# f = plt.plot(1-phi, press_f[i]/1e3, color = 'g', linestyle = '--', label = r'$p_f$')
h = plt.plot(1-phi, (press_h)/1e3, color = 'k', linestyle = '-', label = r'$p_H$')
for i in range(len(shear)):
    c = plt.plot(1-phi, press_c[i]/1e3, linestyle = ':', label = r'$p_c$')

plt.plot(1-phi, press_roth/1e3, color = 'orange', label = 'p_s')
full = mlines.Line2D([], [], color='black', linestyle='-', label='$p_H$')
dot = mlines.Line2D([], [], color='black', linestyle=':', label='$p_\mu$')
full2 = mlines.Line2D([], [], color='orange', linestyle='-', label='$p_s$')
# full = mlines.Line2D([], [], color='black', linestyle='.', label='$p_H$')

ax.legend([dot, full, full2], [ r'$p^c$', r'$p_H$', r'$p_s$'], loc='upper right')
plt.xlabel(r'$1-A$')
plt.ylabel(r'$P$ (kN/m)')
plt.xlim(0, 0.3)
# plt.legend()

plt.savefig('pressure_simplified.png')

plt.figure()
plt.plot(1-phi, press_roth/1e3)
plt.xlabel('dilat')
plt.ylabel('p')
plt.savefig('press_dilat.png')


fig = plt.figure()
ax = plt.axes()
i=1

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# f = plt.plot(1-phi, press_f[i]/1e3, linestyle = '--', label = r'$p_f$')
c = plt.plot(1-phi,  press_c[i]*(press_h)/1e3, linestyle = ':', label = r'$p_c p_h$')
# h = plt.plot(1-phi, (press_f[i]+press_c[i])/1e3, color = 'r', linestyle = '-', label = r'$p_t$')
# h = plt.plot(1-phi, (press_h)/1e3, color = 'k', linestyle = '-', label = r'$p_t$')

# dashed_line = mlines.Line2D([], [], color='black', linestyle='--', label='$p^f$')
# dot = mlines.Line2D([], [], color='black', linestyle=':', label='$p^c$')
# full = mlines.Line2D([], [], color='r', linestyle='-', label='$p_t$')
# full = mlines.Line2D([], [], color='black', linestyle='.', label='$p_H$')

plt.xlabel(r'$1-A$')
plt.ylabel(r'$P$ (kN/m)')
# secax = ax.secondary_xaxis('top')
# secax.set_xlabel('$1-h$ (m)')
# plt.yscale('log')
# plt.grid()
# ax.legend([dashed_line, dot, full], [r'$p^f$', r'$p^c$', r'$p_t$'], loc='upper right')
plt.savefig('pressure_multiplied.png')





# %%
