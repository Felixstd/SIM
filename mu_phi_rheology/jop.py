import numpy as np
import matplotlib.pyplot as plt



I = np.arange(0, 5, 0.01)

mu_s = np.tan(np.rad2deg(20.9))
mu_2 = np.tan(np.rad2deg(32.76))
I_0 = 0.279

mu = mu_s + (mu_2 - mu_s)/(I_0/(I+1))


plt.figure()
plt.plot(I, mu)
plt.savefig('mu.png', dpi = 500)