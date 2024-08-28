import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot



dates = ['1990_01_01_00_15_00', '1990_01_01_00_30_00', '1990_01_01_00_45_00', '1990_01_01_01_00_00', '1990_01_01_01_15_00','1990_01_01_01_30_00',  \
    '1990_01_01_01_45_00', '1990_01_01_02_00_00', '1990_01_01_02_05_00', '1990_01_01_02_10_00', '1990_01_01_02_15_00', '1990_01_01_02_20_00', '1990_01_01_02_25_00']
expno = '57'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/aos/home/fstdenis/SIM/Experiments/'

dx = 1e3
dt = 5*60

datadict_57 = read_data.read_data(expno, dates, outputdir, MuPhi = True)

plot.uniaxial(dates, expno, datadict_57, dx, figdir, MuPhi = True)