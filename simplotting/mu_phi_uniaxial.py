import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import warnings
import argparse

warnings.filterwarnings("ignore")
 
parser = argparse.ArgumentParser(description="For analysis",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-t", "--time", type=int, help = "1 for min, 2 for 10 min, 3 for 2hours and 4 for 10 hours simulation")
parser.add_argument("-e", "--expno", type=int, help = "Experiment number ")
parser.add_argument("-mu0", "--mu0", type=float, help = "Lower bound for mu ")
parser.add_argument("-muinf", "--muinfty", type=float, help = "higher bound for mu ")
args = vars(parser.parse_args())

time = args['time']
expno = args['expno']
mu_0 = args['mu0']
mu_infty = args['muinfty']

if time == 1:
        
        start = 1
        start_k = 61
        #--- Time for 1 minute run ---#
        start = datetime(1990, 1, 1, 0, 1, start)
        intervals = [timedelta(seconds=1)] * 58
        
elif time == 2:
        #--- Time for 10 minutes run ---#
        start = datetime(1990, 1, 1, 0, 0, 30)
        start_k =1
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)] * int((9*60/30))
        
elif time == 3:
        # --- Time for 2 hours 30 minutes run ---#
        start = datetime(1990, 1, 1, 0, 5, 00)
        start_k = 1
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)]*0  + [timedelta(minutes=4)]*0 + [timedelta(minutes=5)]*24 #+ [timedelta(minutes=3)]*1 + [timedelta(minutes=1)]*1
       
elif time == 4:
        #--- Time for 2 hours 30 minutes run ---#
        # start = datetime(1990, 1, 1, 8, 0, 00)
        # start_k = 48
        start = datetime(1990, 1, 1, 0, 10 , 00)
        start_k = 1
        intervals = [timedelta(minutes=10)]*6 #+ [timedelta(minutes=3)]*1 + [timedelta(minutes=1)]*1


dates = [(start + sum(intervals[:i], timedelta())).strftime('%Y_%m_%d_%H_%M_%S') for i in range(len(intervals)+1)]

print(dates)

#------- MuPhi -------#
#---- Uniaxial Tests ----#
outputdir = "/storage/fstdenis/output_sim_MuPhi_Runs_Tests_Dilat/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_Runs_Tests_Dilat/'

#---- Shear Tests ----#
outputdir = "/storage/fstdenis/output_sim_MuPhi_Runs_ShearExperiments/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_Runs_ShearExperiments/'


#------- Mu -------#
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/Mu_ShearExperiments/'



for i in range(expno, expno+1):
    expno = "{:02d}".format(i)

    if not os.path.isdir(figdir+expno):
        os.mkdir(figdir+expno)

    print('Analysising Exp: ', expno)

    dx = 1e3*10
    dy = dx

    angle_phi = 20*np.pi/180


    Ny = 502
    Nx = 202

    # if (Ny == 252 and Nx == 102):
    mask = np.ones((Ny, Nx))
    mask[0, :] = 0
    
#     filemask=outputdir+"mask.dat"
#     maskC = np.genfromtxt(filemask, dtype=None)

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
    # analysis.invariants(dates, expno, datadict, Ny, Nx, dx, maskC, figdir, mu_0, mu_infty)

    # #---------- Plotting ----------#
    plot.uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)
    
#     plot.totdef_uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)
    
