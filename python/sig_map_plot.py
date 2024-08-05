import numpy as np
import math
from copy import copy,deepcopy
import mpl_toolkits.basemap as basemap
import matplotlib.pyplot as plt

nx=256
ny=nx
nxstr=str(nx)
exp="32"
Pmax=30
kN2N=1000.0
fileinfo="_nx0"+nxstr+"_sbc11_10_sbc22_10_sbc12_00."+exp
file1="output/sigI"+fileinfo
file2="output/sigII"+fileinfo

if nx == 256:
  dx=20.0
elif nx == 512:
  dx=10.0
else:
  print 'not coded yet'

xa=dx*np.arange(nx+2)
ya=dx*np.arange(ny+2)

sigI = np.genfromtxt(file1, dtype=None)
sigI = np.true_divide(sigI, kN2N)
fileout="FIGS/p_map_exp"+exp+".png"
plt.pcolor(xa,ya,sigI, cmap='BuPu', vmin=0, vmax=Pmax)
plt.colorbar()
plt.axis([dx, nx*dx, dx, ny*dx])
plt.savefig(fileout)
#plt.show()


plt.figure(2)
sigII = np.genfromtxt(file2, dtype=None)
sigII = np.true_divide(sigII, kN2N)
fileout="FIGS/q_map_exp"+exp+".png"
plt.pcolor(xa,ya,sigII, cmap='BuPu', vmin=0, vmax=Pmax/4)
plt.colorbar()
plt.axis([dx, nx*dx, dx, ny*dx])
plt.savefig(fileout)

plt.show()