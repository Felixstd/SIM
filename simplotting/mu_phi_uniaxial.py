import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
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
        start = datetime(1990, 1, 1, 0, 0, 30)
        start_k = 1
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)]*1  + [timedelta(minutes=4)]*1 + [timedelta(minutes=5)]*23 #+ [timedelta(minutes=3)]*1 + [timedelta(minutes=1)]*1

dates = [(start + sum(intervals[:i], timedelta())).strftime('%Y_%m_%d_%H_%M_%S') for i in range(len(intervals)+1)]


expno = '01'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_Runs/'
for i in range(16, 17):
    expno = "{:02d}".format(i)

    if not os.path.isdir(figdir+expno):
        os.mkdir(figdir+expno)

    print('Analysising Exp: ', expno)

    dx = 1e3*10
    dy = dx

    # mu_0 = 0.1
    # if i < 8:
    #     mu_0 = 0.1
    # else:
    #     mu_0 = 0.36
    mu_0 = 0.2
    mu_infty = 0.9
    angle_phi = 10*np.pi/180


    Ny = 502
    Nx = 202

    # if (Ny == 252 and Nx == 102):
    mask = np.ones((Ny, Nx))
    mask[0, :] = 0

    N_transect = 51
    muphi = 1
    log = 0
    time = np.arange(1, 14, 1)



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
    plot.uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)
    
#     plot.totdef_uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)
    
