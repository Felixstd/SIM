#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import netCDF4 as nc
import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean 
import matplotlib.colors as colors


import os 


def inertial_number(shear, p, rho, h, d_average):
    
    I = d_average * shear * np.sqrt(rho*h/p)
    
    return I 
    
#%%
DataDir = '/storage/fstdenis/RGPS_Data/'
P = 30e3
rho = 910
h = 1
d_average = 1e3


i = 0
for file in os.listdir(DataDir):

    if i == 1:
        rgps = nc.Dataset(DataDir + file ,"r", format="NETCDF4")

        reference_date = rgps.variables['time'].units[14:]
        ref_date = datetime.datetime.strptime(reference_date, "%Y-%m-%dT%H:%M:%SZ")
        
        seconds_since_start = rgps.variables['time'][:]
        
        rgps_x = rgps.variables['latitude'][:]
        rgps_y = rgps.variables['longitude'][:]
        
        for j, time in enumerate(seconds_since_start):
            
            start_date = ref_date + datetime.timedelta(seconds=int(time))
    
            shear = rgps.variables['shear'][j][:]/(60*60*24)
            I = inertial_number(shear, P, rho, h, d_average)
            
            end_date = start_date + datetime.timedelta(days=3)
            
    
            fig = plt.figure(figsize = (10, 4), dpi = 500)
    
            ax1 = fig.add_subplot(121, projection = ccrs.NorthPolarStereo())
            ax2 = fig.add_subplot(122, projection = ccrs.NorthPolarStereo())
    

            ax1.add_feature(cfeature.LAND, color = 'peru')

            cf = ax1.pcolormesh(
                rgps_y,
                rgps_x,
                shear[:-1, :-1],
                cmap=cmocean.cm.thermal,
                norm=colors.Normalize(vmin=np.min(shear), vmax=np.max(shear)),
                transform=ccrs.PlateCarree(),
                zorder=1,
            )
            ax1.add_feature(cfeature.LAND, color = 'peru')
            fig.colorbar(cf, ax = ax1, label = r'$\dot{\epsilon}_{II}$ s$^{-1}$')
            
            ax2.add_feature(cfeature.LAND, color = 'peru')

            cf = ax2.pcolormesh(
                rgps_y,
                rgps_x,
                I[:-1, :-1],
                cmap=cmocean.cm.thermal,
                transform=ccrs.PlateCarree(),
                vmin=np.min(I), vmax=np.max(I),
                zorder=1,
            )
            ax2.add_feature(cfeature.LAND, color = 'peru')
            fig.colorbar(cf, ax = ax2, label = r'$I$')
            
            fig.suptitle('3d average: '+str(start_date)[:-9]+' to '+str(end_date)[:-9])
            plt.savefig('deformation_I_{}.png'.format(str(start_date)[:-9]), bbox_inches = 'tight')


        
        


    i+=1
# %%
