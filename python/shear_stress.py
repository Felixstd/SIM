import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import datetime
import cmocean as cm

from numba import njit

dx = "20"

if int(dx) == 20:
    nx, ny, Deltax = 258, 218, 20e03

# @njit()
def shear(u, v, mask): 
    
    shear = np.zeros_like(u)
    
    
    for i in range(0, nx): 
        for j in range(0, ny): 
            
            dudx = ( u[i+1,j] - u[i,j] ) / Deltax
            dvdy = ( v[i,j+1] - v[i,j] ) / Deltax
            dudy       = 0e0
            dvdx       = 0e0
            
    
            if    (mask[i+1,j] + mask[i-1,j]) == 2 :
                
                dvdx = ( ( v[i+1,j] + v[i+1,j+1] ) - \
                        ( v[i-1,j] + v[i-1,j+1] ) ) / \
                        ( 4e0 * Deltax )
                
            if ( mask[i+1,j] - mask[i-1,j] ) == 1  :
                
                dvdx = ( 1e0 * ( v[i+1,j] + v[i+1,j+1] ) +  \
                        3e0 * ( v[i,j]   + v[i,j+1] ) ) /   \
                        ( 6e0 * Deltax )
                
            if ( mask[i+1,j] - mask[i-1,j]) == -1  :
                
                dvdx = ( -1e0 * ( v[i-1,j] + v[i-1,j+1] ) - \
                        3e0 * ( v[i,j]   + v[i,j+1] ) ) /   \
                        ( 6e0 * Deltax )
            

            if     (mask[i,j+1] + mask[i,j-1]) == 2  :
                
                dudy = ( ( u[i,j+1] + u[i+1,j+1] ) -        \
                        ( u[i,j-1] + u[i+1,j-1] ) ) /       \
                        ( 4e0 * Deltax )
                
            if (mask[i,j+1] - mask[i,j-1]) == 1  :
                
                dudy = ( 1e0 * ( u[i,j+1] + u[i+1,j+1] ) +  \
                        3e0 * ( u[i,j]   + u[i+1,j] ) ) / \
                        ( 6e0 * Deltax )
                
            if ( mask[i,j+1] - mask[i,j-1] == -1 ) :
                
                dudy = ( -1e0 * ( u[i,j-1] + u[i+1,j-1] ) - \
                        3e0 * ( u[i,j]   + u[i+1,j] ) ) / \
                        ( 6e0 * Deltax )
    
    
            shear[i,j] = np.sqrt(( dudx - dvdy )**2e0 \
                       + ( dudy + dvdx )**2e0 )
    
    
    
    return shear
    






outputdir = "/aos/home/fstdenis/SIM/output/"

exp = "16"
dates = ["1994_01_01_00_00", "1994_01_02_00_00", "1994_01_03_00_00"]

filemask="./MASKS/mask"+dx+"_1_0_1.dat"
mask = np.genfromtxt(filemask, dtype=None)

    
k = 1
for date in dates:

    file_shear = outputdir+'shear'+ date + '_k000'+str(k)+'.' + exp

    shear_arctic = np.loadtxt(file_shear, dtype=None)/(60*60*24)
    
    shear_arctic_plot = np.zeros(shear_arctic.shape)
    shear_arctic_plot[:,:] = shear_arctic[:,:]
    shear_arctic_plot[mask == 0] = np.nan
    var_nonmasked = np.copy(shear_arctic_plot)
    shear_arctic_plot = np.ma.masked_invalid(shear_arctic_plot)
    
    datetime_format = '%Y_%m_%d_%H_%M'
    mydates = datetime.datetime.strptime(date, datetime_format)

    fileout="shear" + date +"_" + exp +".png"


    cmap = mpl.cm.get_cmap("inferno").copy()
    cmap.set_bad(color = 'gray', alpha = 1.)

    cmap_ice = cm.cm.ice.copy()
    cmap_ice.set_bad(color = 'gray', alpha = 1.)

    fig = plt.figure(figsize = (10, 8))


    pc = plt.pcolormesh(shear_arctic_plot, cmap=cmap, vmin = 0, vmax = 1e-6)
    fig.colorbar(pc, label = r'$\dot{\epsilon}_{II}$')
    plt.savefig(fileout, dpi = 500, bbox_inches = 'tight')
    
    k+=1
    # varplot[:,:] = var[:,:]
    # varplot[mask == 0] = np.nan
# var_nonmasked = np.copy(varplot)
# varplot = np.ma.masked_invalid(varplot)




# print(var)

# varplot = np.zeros(var.shape)

