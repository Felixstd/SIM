import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import warnings

warnings.filterwarnings("ignore")


minute = 0
two_hours = 1
ten_minutes = 0

if minute:
        
        start = 1
        start_k = 61
        #--- Time for 1 minute run ---#
        start = datetime(1990, 1, 1, 0, 1, start)
        intervals = [timedelta(seconds=1)] * 58
        
elif ten_minutes:
        #--- Time for 10 minutes run ---#
        start = datetime(1990, 1, 1, 0, 0, 30)
        start_k =1
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)] * int((9*60/30))
        
elif two_hours:
        #--- Time for 2 hours 30 minutes run ---#
        start = datetime(1990, 1, 1, 0, 00, 30)
        start_k =1
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)]*1  + [timedelta(minutes=4)]*0 + [timedelta(minutes=5)] *0

dates = [(start + sum(intervals[:i], timedelta())).strftime('%Y_%m_%d_%H_%M_%S') for i in range(len(intervals)+1)]
print(dates)



expno = '81'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_NewP_NewDilat/'


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
dx = 1e3*10
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
log = 0
time = np.arange(1, 14, 1)
# time = np.arange(1, 4, 1)

mu_0 = np.tan(5*np.pi/36)

mu_infty = np.tan(13*np.pi/36)

mu_0 = 0.1

mu_infty = 0.9

#_--------- READING DATA ----------#
datadict = read_data.read_data(expno, start_k, dates, outputdir, MuPhi = muphi)
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