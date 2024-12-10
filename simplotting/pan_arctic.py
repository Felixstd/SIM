import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import warnings
import os


expno = '09'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_Runs/'

if not os.path.isdir(figdir+expno):
    os.mkdir(figdir+expno)

#---- Constants ----#
dx = 1e3*10
dy = dx
mu_0 = 0.36
mu_infty = 0.9
angle_phi = 10*np.pi/180

Ny = 438
Nx = 518

muphi = 1

#----- Dates -----#

start = datetime(1990, 1, 8, 0, 00, 00)
start_k = 85
# Time interval (30 seconds initially, then 5-minute steps)
intervals = [timedelta(hours=1, minutes=40)]*0 + [timedelta(hours=2)]*20
dates = [(start + sum(intervals[:i], timedelta())).strftime('%Y_%m_%d_%H_%M_%S') for i in range(len(intervals)+1)]

print(dates)


datadict = read_data.read_data(expno, start_k, dates, outputdir, MuPhi = muphi)

plot.plot_data_panarctic(figdir+expno+'/', datadict, dates, Nx, Ny, 1e3*10, MuPhi = True)