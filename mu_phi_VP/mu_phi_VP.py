#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
def inertial_number(h,shear, P, rhoice, d_average):
    
    return d_average * shear * np.sqrt( ( rhoice * h ) / P)

def friction_coefficient(eta, shear, P): 
    
    return 2 * eta * shear / P
    

def load_params_muphi(filename, mask):

    var = np.loadtxt(filename, dtype=None)

    varplot = np.zeros(var.shape)
    
    if var.shape != mask.shape:
        mask = mask[1:, 1:-1]


    varplot[:,:] = var[:,:]
    varplot[mask == 0] = np.nan
    var_nonmasked = np.copy(varplot)
    varplot = np.ma.masked_invalid(varplot)

    return varplot, var_nonmasked

#%%

outputdir = '../output/'
exp = "21"
dates = ["1991_01_01_00_00"]
dx = "20"

filemask="../python/MASKS/mask"+dx+"_1_0_1.dat"
mask = np.genfromtxt(filemask, dtype=None)

rhoice = 900
d_mean = 20e3*7


for date in dates:
    file_shear = outputdir + "shear" + date + '_k0001' + '.' + exp
    file_p     = outputdir + "p" + date + '.' + exp
    file_eta     = outputdir + "eta" + date + '.' + exp
    file_h     = outputdir + "h" + date + '.' + exp
    file_A     = outputdir + "A" + date + '.' + exp
    

    shear, shear_nonmask = load_params_muphi(file_shear, mask)
    h, h_nonmask = load_params_muphi(file_h, mask)
    eta, eta_nonmask = load_params_muphi(file_eta, mask)
    p, p_nonmask = load_params_muphi(file_p, mask)
    
    shear = shear[1:, 1:-1]/(60*60*24)
    eta = eta[1:, 1:-1]
    h  = h[1:, 1:-1]
    
    
    I = inertial_number(h.flatten(), shear.flatten(), p.flatten(), rhoice, d_mean)
    mu = friction_coefficient(eta.flatten(), shear.flatten(), p.flatten())
    
    
plt.figure()
plt.scatter(I, mu, color = 'k', marker = 'o')
plt.xscale('log')
plt.xlim(1e-5, 1)
plt.xlabel(r'$I$')
plt.ylabel(r'$\mu(I)$')
plt.title('With V-P')
plt.savefig('mu_I_VP_20km.png', dpi = 500, bbox_inches = 'tight')
# %%
