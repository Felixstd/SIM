import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot

# 1990-01-01:00:00:00  starting date
# 1990-01-01:05:30:00  end date
# 1990-01-01:01:00:00  posting date
# 1990-01-01:01:30:00  posting date
# 1990-01-01:02:00:00  posting date
# 1990-01-01:02:30:00  posting date
# 1990-01-01:03:00:00  posting date
# 1990-01-01:03:30:00  posting date
# 1990-01-01:04:00:00  posting date
# 1990-01-01:04:30:00  posting date
# 1990-01-01:05:00:00  posting date
# 1990-01-01:05:05:00  posting date
# 1990-01-01:05:10:00  posting date
# 1990-01-01:05:15:00  posting date
# 1990-01-01:05:20:00  posting date
# 1990-01-01:05:25:00  posting date

# 1                    input namelist
# 0                    restart
# 15                   exp version
# 1990-01-01:00:00:00  starting date
# 1990-01-01:02:30:00  end date
# 1990-01-01:00:15:00  posting date
# 1990-01-01:00:30:00  posting date
# 1990-01-01:00:45:00  posting date
# 1990-01-01:01:00:00  posting date
# 1990-01-01:01:15:00  posting date
# 1990-01-01:01:30:00  posting date
# 1990-01-01:01:45:00  posting date
# 1990-01-01:02:00:00  posting date
# 1990-01-01:02:15:00  posting date
# 1990-01-01:02:20:00  posting date
# 1990-01-01:02:25:00  posting date
# 1990-01-01:02:28:00  posting date
# 1990-01-01:02:29:00  posting date
# stop                 end of post date string

# For 1 to 13
dates = ['1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04', 
         '1990_01_01_00_00_05', '1990_01_01_00_00_06', '1990_01_01_00_00_07', '1990_01_01_00_00_08', 
         '1990_01_01_00_00_09', '1990_01_01_00_00_10', '1990_01_01_00_00_11', '1990_01_01_00_00_12', 
         '1990_01_01_00_00_13', '1990_01_01_00_00_14']



expno = '99'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/aos/home/fstdenis/SIM/Experiments/Sensitivity_Experiment/'


Hibler_P = 0
MuPhi_P  = 1
Hibler   = 0
compression = 0


if Hibler_P:
    figdir = figdir+'Hibler_P/'

if MuPhi_P:
    figdir = figdir+'MuPhi_P/'
    
if Hibler:
    figdir = figdir+'VP/'
    
if compression:
    figdir = '/aos/home/fstdenis/SIM/Experiments/Ice_top_Ice/'
    



dx = 1e3/2
# dt = 5*60

datadict = read_data.read_data(expno, dates, outputdir, MuPhi = True)

plot.uniaxial(dates, expno, datadict, dx, figdir, MuPhi = True)


#Exp 62
# plot.histograms(datadict, 0.4, 0.5, dates, figdir, expno)