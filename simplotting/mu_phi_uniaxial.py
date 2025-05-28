import numpy as np
import xarray as xr
from simplotting.data import read_data
from simplotting.data.readnamelist import namelist
from simplotting.utils.TimeUtilities import TimeUtility
from simplotting.plotting import plot
from simplotting.analysis import analysis

import matplotlib.pyplot as plt


import configparser
import os
import warnings
import matplotlib.lines as mlines

warnings.filterwarnings("ignore")


#---- Reading the Configuration file for the analysis ----#

config_exp = configparser.ConfigParser()
config_exp.read('./namelistAnalysis')

#---- Setting the Dates ----#
Dates_Config = TimeUtility(configuration = config_exp['Time'])
dates_pre = list(TimeUtility.dates_analysis(Dates_Config, nmax = Dates_Config.nsteps))
dates = [date.strftime('%Y_%m_%d_%H_%M_%S') for date in dates_pre ]

#---- Reading the Parameters for the analysis
Parameters = namelist(configuration_exp  = config_exp['Experiment'], 
                      configuration_rheo = config_exp['Rheology'], 
                      configuration_fig  = config_exp['Figures'], 
                      configuration_time = config_exp['Time'])

expnos = Parameters.expno



#--- Looping through the different experiments to analyse ---#

u_exp = []
v_exp = []
for i, expno in enumerate(expnos):
    expno = "{:02d}".format(expno)

    if not os.path.isdir(Parameters.figdir+expno):
        os.mkdir(Parameters.figdir+expno)

    print('Analysising Exp: ', expno)

    N_transect = 51
    log = 0
    time = np.arange(1, len(dates), 1)*Parameters.dt

    #_--------- READING DATA ----------#
    if Parameters.read_all:
        datadict = read_data.read_data(expno, 
                                   Parameters.startk, 
                                   dates, 
                                   Parameters.outputdir, 
                                   MuPhi = Parameters.muphi)
            
        if Parameters.muphi:
            
            divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, zeta_tot,eta_tot, \
                uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()
        
        else: 
            
            divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot,sigI_tot, sigII_tot, zeta_tot, eta_tot, uair_tot, vair_tot = \
                datadict.values()
                
        shear = shearI_tot[0]
        print('shear',np.vstack([shear[100:, 0],shear[100:, -2]]).T)
        u = u_tot[0]
        v = v_tot[0]
        
        u_exp.append(u)
        v_exp.append(v)
    

    if Parameters.savevar:
        
        file_A = "./A_time_{}.nc".format(expno)
        file_sigI = "./sigI_time_{}.nc".format(expno)
        file_sigII = "./sigII_time_{}.nc".format(expno)
        
        if os.path.exists(file_A) == False:
            
            read_data.read_data_var('A', 
                            expno, 
                            Parameters, 
                            Parameters.startk, 
                            dates, 
                            Parameters.outputdir, 
                            'A_time_{}.nc'.format(expno))
        
        # if os.path.exists(file_sigI) == False:
    
        #     read_data.read_data_var('sigInorm', 
        #                     expno, 
        #                     Parameters, 
        #                     Parameters.startk, 
        #                     dates, 
        #                     Parameters.outputdir, 
        #                     'sigI_time_{}.nc'.format(expno))

        # if os.path.exists(file_sigII) == False:
        #     read_data.read_data_var('sigIInorm', 
        #                     expno, 
        #                     Parameters, 
        #                     Parameters.startk, 
        #                     dates, 
        #                     Parameters.outputdir, 
        #                     'sigII_time_{}.nc'.format(expno))
    
    
    # #---------- Analysing Wind Forcing----------#
    # analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', muphi)

    # # #---------- Analysing Invariants ----------#
    # analysis.invariants(dates, expno, datadict, Ny, Nx, dx, maskC, figdir, mu_0, mu_infty)

    # #---------- Plotting ----------#
    if Parameters.plotfields:
        plot.uniaxial(dates, 
                  expno, 
                  datadict, 
                  Parameters.dx, 
                  Parameters.figdir, 
                  Parameters.mu0, 
                  Parameters.muinf, 
                  Parameters.microangle, 
                  MuPhi = Parameters.muphi, 
                  log = log)
    

#     mean_mu_time, mean_div_time, mean_shear_time = analysis.mean_values(dates, datadict, dx, dt)
#     plot.plot_mean(mean_mu_time, time, r'$\langle\mu\rangle$', 'mean_mu.png', expno, figdir)
#     plot.plot_mean(mean_div_time, time, r'$\langle\ \dot{\epsilon}_I \rangle$ (day$^{-1}$)', 'mean_div.png', expno, figdir)
#     plot.plot_mean(mean_shear_time, time, r'$\langle\ \dot{\epsilon}_\mathrm{II} \rangle$ (day$^{-1}$)', 'mean_shear.png', expno, figdir)
if Parameters.veltransect:
    print('Analysis of Tracer (u, v) transects')
    analysis.velocity_transects(datadict, dates, Parameters.dx,Parameters.figdir, expno, MuPhi = Parameters.muphi)

if Parameters.tracertransect:
    print('Analysis of Tracer (A, h) transects')
    analysis.tracers_transect(datadict, dates, Parameters.dx,Parameters.figdir, expno, MuPhi = Parameters.muphi)
#     
    # plot.totdef_uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)

if Parameters.plotdilatation:
    plot.dilatation('./', expnos, Parameters)

if Parameters.plotsingle:
    dataset_A = xr.open_dataset('A_time_{:02d}.nc'.format(Parameters.expno[0]))
    A_ice_tot= dataset_A["data"]
    for i, A in enumerate(A_ice_tot):
        plot.plot_single(Parameters.Nx,
                        Parameters.Ny, 
                        Parameters.dx, 
                        1-A, 
                        'A', 
                        Parameters.figdir+'{}/SeaIceConc/A_{}.png'.format(expno, dates[i]))

u = u_exp[0] - u_exp[1]
v = v_exp[0] - v_exp[1]
eps = np.finfo(u.dtype).eps
print(np.float64(np.sqrt(np.sum(v**2))), eps)

# print(np.shape(v), np.shape(u))
# x_u = np.arange(501) + 0.5  # 0.5 to 499.5 (500 values)
# y_u = np.arange(200)            # 0 to 199
# X_u, Y_u = np.meshgrid(x_u, y_u)

# # Grid for v (on north-south faces)
# x_v = np.arange(500)            # 0 to 499
# y_v = np.arange(201) + 0.5 # 0.5 to 199.5 (200 values)
# X_v, Y_v = np.meshgrid(x_v, y_v)
# u = u[:, :-1]
# v = v[:-1, :]
# # Create a plot
# plt.figure(figsize=(12, 6))
# # plt.quiver(X_u[180:, :], Y_u[180:, :], u[180:, :], np.zeros_like(u)[180:, :], color='r', scale=20, label='u (zonal)', alpha=0.6)
# # plt.quiver(X_v[180:, :], Y_v[180:, :], np.zeros_like(v)[180:, :], v[180:, :], color='b', scale=20, label='v (meridional)', alpha=0.6)

# # plt.figure()
# plt.quiver(u[180::2,::5], v[180::2, ::5])
# plt.savefig('u_diff_freedrift.png')


# u  =u_tot[-1]
# print(np.shape(u))
# # u_split = np.hsplit(u, 2)
# y = np.shape(u)[1] // 2  # 250
# u_split_1 = u[:, :y]
# u_split_2 = u[:, y+1:]

# # print(np.shape(u_split))
# plt.figure()
# plt.pcolor(u_split_1 - u_split_2)
# plt.colorbar()
# plt.savefig('error_vp.png')





            # h = h_tot[0]
    # u = u_tot[0]
    # v = v_tot[0]
    # p = p_tot[0]
    # # mu = muI_tot[1]
    # A = A_tot[0]
    # div = divergence_tot[0]
    # shear = shear_tot[0]
    # zeta = zeta_tot[0]
    # eta = eta_tot[0]
        # print(np.shape(shearI))
        # print('I', np.vstack([I[100:, 1],I[100:, -1]]).T)
        # print('shearI', shearI[:, -1])
        # print('shearI', shearI[:, 1])

        # with np.printoptions(threshold=np.inf):
            # print('shearI', np.where(shearI == -999)[0])
        # print(np.where(shearI == -999)[0], shearI[0, -1])
        # print('shearI', np.vstack([shearI[201, 0],shearI[201, 501]]).T)
        # print('mu', np.vstack([mu[100:, 1],mu[100:, -1]]).T)
        # print('zeta', np.vstack([zeta[100:, 1],zeta[100:, -1]]).T)
        # print('eta', np.vstack([eta[100:, 1],eta[100:, -1]]).T)
    # print(np.shape(u))
    # print('h', np.vstack([h[100:, 1],h[100:, -1]]).T)
    # print('A', np.vstack([A[100:, 1],A[100:, -1]]).T)
    # print('u', np.vstack([u[100:, 0],u[100:, -1]]).T)
    # print('v', np.vstack([v[100:, 0],v[100:, -1]]).T)
    # print('p', np.vstack([p[100:, 0],p[100:, -1]]).T)
    # print('div',np.vstack([div[100:, 1],div[100:, -1]]).T)
    # print('shear',np.vstack([shear[100:, 1],shear[100:, -1]]).T)
    # print('shear',np.vstack([shear[100:, 0],shear[100:, -2]]).T)
    # print('zeta', np.vstack([zeta[100:, 1],zeta[100:, -1]]).T)
    # print('eta', np.vstack([eta[100:, 1],eta[100:, -1]]).T)