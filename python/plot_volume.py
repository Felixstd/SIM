import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from datetime import datetime

exp1="85"
file1="outputSIM/sea_ice_volume." +exp1
exp2="86"
file2="outputSIM/sea_ice_volume." +exp2

lista1 = [] # dates
listvol1 = [] # volume
with open(file1, 'r') as file:
   for row in file:
      a, b = row.split()
      lista1.append(int(a))
      listvol1.append(float(b))

datetp=np.array(lista1)

dtime1 = ["" for i in range(len(datetp))]

for i in range(len(datetp)):
   dstr=str(datetp[i])
   dtime1[i]=datetime.strptime(dstr, "%Y%m%d")  

lista2 = [] # dates                                                  
listvol2 = [] # volume

with open(file2, 'r') as file:
   for row in file:
      a, b = row.split()
      lista2.append(int(a))
      listvol2.append(float(b))

datetp=np.array(lista2)

dtime2 = ["" for i in range(len(datetp))]

for i in range(len(datetp)):
   dstr=str(datetp[i])
   dtime2[i]=datetime.strptime(dstr, "%Y%m%d")

#plt.plot_date(dtime,listb)
plt.plot(dtime1,listvol1, label='upwind')
plt.plot(dtime2,listvol2, label='SL')
plt.legend(loc="upper right")
plt.ylabel('Sea ice volume (10$^3$ km$^3$)', fontsize=18)
plt.locator_params(axis='y', nbins=4) 
#plt.locator_params(axis='x', nbins=6)
plt.xlim(datetime.strptime(str(200201), '%Y%m'), datetime.strptime(str(200301), '%Y%m'))
plt.ylim(0, 30)
fileout="FIGS/volume_" + exp1 + "_" +exp2 +".png"
#plt.axis([0, h.shape[1], 0, h.shape[0]])
plt.savefig(fileout)

plt.show()
