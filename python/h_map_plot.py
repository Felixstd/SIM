import numpy as np
import math
import matplotlib.pyplot as plt
import datetime

exp="09"
dx="40"
date="1990_01_02_00_00"
outputdir = "/aos/home/fstdenis/SIM/output/"
file1= outputdir+'h' + date + "." +exp


h = np.genfromtxt(file1, dtype=None)

filemask="./MASKS/mask"+dx+"_1_0_1.dat"
mask = np.genfromtxt(filemask, dtype=None)

print(h.shape[0])
print(h.shape[1])

hplot = np.zeros(h.shape)
hplot[:,:]=h[:,:]
hplot[mask==0]=np.nan
hplot = np.ma.masked_invalid(hplot)

datetime_format = '%Y_%m_%d_%H_%M'
mydates = datetime.datetime.strptime(date, datetime_format)

fileout="h" + date +"_" + exp +".png"
cmap = plt.get_cmap('BuPu')
cmap.set_bad(color = 'gray', alpha = 1.)
plt.pcolormesh(hplot, cmap=cmap, vmin=0, vmax=4)
plt.colorbar(label = 'thickness (m)')
plt.axis([0, h.shape[1], 0, h.shape[0]])
plt.title(mydates)
plt.savefig(fileout, dpi = 500, bbox_inches = 'tight')

# plt.savefig('test_2.png')
