import numpy as np
import math
import matplotlib.pyplot as plt

exp1="01"
exp2="02"
dx="80"
date="1998_01_01_00_00"
file1="outputSIM/h" + date + "." +exp1
h1 = np.genfromtxt(file1, dtype=None)
file2="outputSIM/h" + date + "." +exp2
h2 = np.genfromtxt(file2, dtype=None)

filemask="../MASKS/mask"+dx+"_1_0_1.dat"
mask = np.genfromtxt(filemask, dtype=None)

print(h1.shape[0])
print(h1.shape[1])

hdiff = np.zeros(h1.shape)
hdiff[:,:]=h2[:,:]-h1[:,:]
hdiff[mask==0]=np.nan
hdiff = np.ma.masked_invalid(hdiff)

fileout="FIGS/hdiff" + date +"_" + exp2 + "_" + exp1 + ".png"
cmap = plt.get_cmap('bwr')
cmap.set_bad(color = 'gray', alpha = 1.)
plt.pcolormesh(hdiff, cmap=cmap, vmin=-1, vmax=1)
plt.colorbar()
plt.axis([0, h1.shape[1], 0, h1.shape[0]])
plt.savefig(fileout)

plt.show()
