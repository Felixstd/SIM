import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore")


dates = ['1990_01_01_00_00_01', 
        '1990_01_01_00_00_02', 
        '1990_01_01_00_00_03', 
        '1990_01_01_00_00_04',
        '1990_01_01_00_00_05',
        '1990_01_01_00_00_06',  
        '1990_01_01_00_00_07',
        '1990_01_01_00_00_08', 
        '1990_01_01_00_00_09',
        '1990_01_01_00_00_10', 
        '1990_01_01_00_00_11', 
        '1990_01_01_00_00_12', 
        '1990_01_01_00_00_13',
        '1990_01_01_00_00_14',
        '1990_01_01_00_00_15',
        '1990_01_01_00_00_16']
        # '1990_01_01_00_00_17', 
        # '1990_01_01_00_00_18', 
        # '1990_01_01_00_00_19',
        # '1990_01_01_00_00_20',
        # '1990_01_01_00_00_21', 
        # '1990_01_01_00_00_22', 
        # '1990_01_01_00_00_23', 
        # '1990_01_01_00_00_24', 
        # '1990_01_01_00_00_25', 
        # '1990_01_01_00_00_26', 
        # '1990_01_01_00_00_27', 
        # '1990_01_01_00_00_28', 
        # '1990_01_01_00_00_29']

expno = '09'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/aos/home/fstdenis/SIM/Experiments/MuPhi_DilatationLaw/'


Hibler_P = 0
MuPhi_P  = 0
Hibler   = 0
compression = 0
phi = 0


# if Hibler_P:
#     figdir = figdir+'Hibler_P/'

# if MuPhi_P:
#     figdir = figdir+'MuPhi_P/'
    
# if Hibler:
#     figdir = figdir+'VP/'
    
# if compression:
#     figdir = '/aos/home/fstdenis/SIM/Experiments/Ice_top_Ice/'
    
# if phi:
#     figdir = figdir+'Phi/'
    


# dx = 1e3/2
dx = 2e3
dy = dx
# Nx = 502
# Ny = 502
# Nx = 202

# Ny = 252
# Nx = 102

Ny = 252
Nx = 102

N_transect = 51
muphi = 1
log = 1
time = np.arange(1, 14, 1)
# time = np.arange(1, 4, 1)

# mu_0 = np.tan(5*np.pi/36)

# mu_infty = np.tan(13*np.pi/36)

mu_0 = 0.1

mu_infty = 0.9

#_--------- READING DATA ----------#
datadict = read_data.read_data(expno, dates, outputdir, MuPhi = muphi)
# print(datadict)
# #---------- Analysing Wind Forcing----------#
# analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', muphi)

# # #---------- Analysing Invariants ----------#
# analysis.invariants(datadict, N_transect, time, figdir+expno+'/', muphi)

# #---------- Plotting ----------#
plot.uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, MuPhi = muphi, log = log)

# if muphi:
#     divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, \
#         zeta_tot, uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()
# else: 
        
#         divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, uair_tot, vair_tot = \
#             datadict.values()

# k = -1
# tot_deformation =np.sqrt(divergence_tot[k]**2 + shear_tot[k]**2)



# plot.plot_single(Nx, Ny, dx, tot_deformation, r'Total Deformation (1/day)', figdir+expno+'/totdeformation_{}.png'.format(dates[k]))




# dates = ['1990_01_01_00_15_00',
#         '1990_01_01_00_30_00',
#         '1990_01_01_00_45_00',
#          '1990_01_01_01_00_00',
#          '1990_01_01_01_15_00',
#          '1990_01_01_01_30_00',
#          '1990_01_01_01_45_00',
#          '1990_01_01_02_00_00',
#          '1990_01_01_02_15_00',
#          '1990_01_01_02_20_00',
#          '1990_01_01_02_25_00',
#          '1990_01_01_02_28_00',
#          '1990_01_01_02_29_00']