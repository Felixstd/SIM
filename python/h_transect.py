import numpy as np
import math
import matplotlib.pyplot as plt

exp="30"
xtrans=0 # if = 0 : ytransect
pos=25 #x or y position of transect
date="1990_03_02_00_00"
file1="outputSIM/h" + date + "." +exp
h = np.genfromtxt(file1, dtype=None)

print(h.shape)

if xtrans == 1:
        htrans = np.zeros(h.shape[1])
        htrans=h[pos,:]
else:
        htrans = np.zeros(h.shape[0])
        htrans=h[:,pos]


print(htrans.shape[0])

plt.plot(htrans, markersize=6, color=[1, 0.33, 0.16])
plt.axis([0, htrans.shape[0], 0, 1.5])    
plt.ylabel('Sea ice thickness (m)', fontsize=18)
plt.xlabel('grid cell', fontsize=18)
fileout="FIGS/htrans" + date +"_" + exp +".png"
plt.savefig(fileout)

plt.show()
