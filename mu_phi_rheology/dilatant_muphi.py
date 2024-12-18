import numpy as np
import matplotlib.pyplot as plt



wspeed = 1
Tramp = 60*60*2

tstep = np.arange(0, 200)
Deltax = 1e4


Tramp = 200*Deltax/3

rampfactor = np.tanh(tstep*Deltax/Tramp)


rampfactor = np.exp(-((abs(tstep-100))*Deltax/Tramp)**6)
vspeed = wspeed*rampfactor

plt.figure()
plt.plot(tstep*Deltax/(1e4), vspeed)
plt.savefig('speed.png')