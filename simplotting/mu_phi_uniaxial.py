import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore")


date_range = np.arange(1, 4)
dates = []
for data in date_range:
        dates.append('1990_01_01_00_00_{:02d}'.format(data))


# dates = ['1990_01_01_00_05_00',
#         '1990_01_01_00_15_00', 
#         '1990_01_01_00_20_00',
#         '1990_01_01_00_30_00', 
#         '1990_01_01_00_40_00', 
#         '1990_01_01_00_45_00', 
#         '1990_01_01_00_50_00',
#         '1990_01_01_01_00_00',
#         '1990_01_01_01_10_00',
#         '1990_01_01_01_15_00']
        # '1990_01_01_01_25_00',    
        # '1990_01_01_01_30_00',   
        # '1990_01_01_01_40_00'] 
        # '1990_01_01_01_45_00',   
        # '1990_01_01_02_00_00',   
        # '1990_01_01_02_15_00',   
        # '1990_01_01_02_20_00',
        # '1990_01_01_02_25_00',   
        # '1990_01_01_02_28_00',   
        # '1990_01_01_02_29_00']  

dates = ['1990_01_01_00_00_30',
        '1990_01_01_00_01_00',
        '1990_01_01_00_01_30',
        '1990_01_01_00_02_00',
        '1990_01_01_00_02_30',
        '1990_01_01_00_03_00',
        '1990_01_01_00_03_30',
        '1990_01_01_00_04_00']
#         '1990_01_01_00_04_30',
#         '1990_01_01_00_05_00', 
#         '1990_01_01_00_05_30',
#         '1990_01_01_00_06_00', 
#         '1990_01_01_00_06_30',
#         '1990_01_01_00_07_00',
#         '1990_01_01_00_07_30',
#         '1990_01_01_00_08_00', 
#         '1990_01_01_00_08_30',
#         '1990_01_01_00_09_00']

expno = '92'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_DilatationLaw/'


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
# dx = 10e3
dy = dx
# Nx = 502
# Ny = 502
# Nx = 202

# Ny = 252
# Nx = 102

Ny = 502
Nx = 202

# if (Ny == 252 and Nx == 102):
mask = np.ones((Ny, Nx))
mask[0, :] = 0

N_transect = 51
muphi = 1
log = 1
time = np.arange(1, 14, 1)
# time = np.arange(1, 4, 1)

mu_0 = np.tan(5*np.pi/36)

mu_infty = np.tan(13*np.pi/36)

mu_0 = 0.1

mu_infty = 0.9

#_--------- READING DATA ----------#
datadict = read_data.read_data(expno, dates, outputdir, MuPhi = muphi)
if muphi:
    divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, \
        zeta_tot, eta_tot, uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()
# print(datadict)
# #---------- Analysing Wind Forcing----------#
# analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', muphi)

# # #---------- Analysing Invariants ----------#
# analysis.invariants(datadict, N_transect, time, figdir+expno+'/', muphi)

# #---------- Plotting ----------#
plot.uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, MuPhi = muphi, log = log)

# else: 
        
#         divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, uair_tot, vair_tot = \
#             datadict.values()

# k = -1
# tot_deformation =np.sqrt(divergence_tot[k]**2 + shear_tot[k]**2)



# plot.plot_single(Nx, Ny, dx, tot_deformation, r'Total Deformation (1/day)', figdir+expno+'/totdeformation_{}.png'.format(dates[k]))


# k = 9
# u, v, zeta, eta, p = u_tot[k], v_tot[k], zeta_tot[k], eta_tot[k], p_tot[k]

# print(p)

# sigI, sigII = analysis.compute_invariant_stresses(u, v, zeta, eta, mask, p, Ny, Nx, dx)
# print(sigI.flatten())

# # sigII[sigII>100000] = np.nan
# plt.figure()
# plt.scatter(sigI.ravel(), sigII.ravel(), s= 2, color= 'r')
# plt.plot(np.unique(sigI.ravel()), np.abs(np.tan(5*np.pi/36)*np.unique(sigI.ravel())))
# plt.plot(np.unique(sigI.ravel()), np.abs(np.tan(13*np.pi/36)*np.unique(sigI.ravel())))
# # plt.xscale('log')
# plt.yscale('log')
# # plt.ylim(1e-3, 1e5)
# # plt.xlim(-10, 10)
# plt.savefig('stress.png')