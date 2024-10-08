import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore")


# For 1 to 13
dates = ['1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04', 
         '1990_01_01_00_00_05', '1990_01_01_00_00_06', '1990_01_01_00_00_07', '1990_01_01_00_00_08', 
         '1990_01_01_00_00_09', '1990_01_01_00_00_10', '1990_01_01_00_00_11', '1990_01_01_00_00_12', 
         '1990_01_01_00_00_13', '1990_01_01_00_00_14']

expno = '99'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/aos/home/fstdenis/SIM/Experiments/MuPhi_GoodWay_Convergence/'


Hibler_P = 0
MuPhi_P  = 0
Hibler   = 0
compression = 0
phi = 1


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
    


dx = 1e3/2
dy = dx
Nx = 202
Ny = 502
N_transect = 101
time = np.arange(1, 15, 1)

mu_0 = np.tan(5*np.pi/36)

mu_infty = np.tan(13*np.pi/36)

#---------- READING DATA ----------#
datadict = read_data.read_data(expno, dates, outputdir, MuPhi = True)

#---------- Analysing Wind Forcing----------#
analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', True)

#---------- Analysing Invariants ----------#
analysis.invariants(datadict, N_transect, time, figdir+expno+'/', True)

#---------- Plotting ----------#
plot.uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, MuPhi = True)


