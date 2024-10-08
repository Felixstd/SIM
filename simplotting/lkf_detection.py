import numpy as np
from simplotting.data import read_data
from simplotting.LKF import lkf
from simplotting.LKF import lkf_hutter
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

dates = ['1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04', 
         '1990_01_01_00_00_05', '1990_01_01_00_00_06', '1990_01_01_00_00_07', '1990_01_01_00_00_08', 
         '1990_01_01_00_00_09', '1990_01_01_00_00_10', '1990_01_01_00_00_11', '1990_01_01_00_00_12', 
         '1990_01_01_00_00_13', '1990_01_01_00_00_14']


expno = '11'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/aos/home/fstdenis/SIM/Experiments/MuPhi_GoodWay_Convergence/'

max_kernel=5
min_kernel=1
dog_thres=0
dis_thres=4
ellp_fac=3
angle_thres=35
eps_thres=0.5
lmin=4
max_ind=500 
use_eps=False
skeleton_kernel=0
k = -4

datadict = read_data.read_data(expno, dates, outputdir, MuPhi = True)

divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, \
    uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()

shear_day = shear_tot[k][:501, 20:180]
div_day   = divergence_tot[k][:501, 20:180]
mu_day = muI_tot[k][:501, 20:180]

#find the range of mu for this day
mu_unique = np.unique(mu_day)
print(np.shape(mu_unique))


eps_tot = np.sqrt(shear_day**2 + div_day**2)

# print(np.shape(eps_tot))
# plt.pcolormesh(eps_tot)
# plt.savefig('text.png')

pairs_lkf, angles_lkf = lkf_hutter.lkf_detection(eps_tot, max_kernel, min_kernel, dog_thres, dis_thres, \
    ellp_fac, angle_thres, eps_thres, lmin, max_ind, use_eps, skeleton_kernel)

angles_lkf, idx_unique = np.unique(angles_lkf, return_index = True)
pairs_lkf = pairs_lkf[idx_unique[0]]

fracture_angles_range = lkf.mohr_fracture_angle(mu_unique)
mu_angles_lkf = lkf.mohr_friction_coeffcient(angles_lkf)

plt.figure()
plt.scatter(mu_unique, fracture_angles_range, color = 'k')
plt.scatter(mu_angles_lkf, angles_lkf, color = 'r', label = 'LKF angle')
plt.xlabel(r'$\mu$')
plt.ylabel(r'$2\theta$')
plt.legend()
plt.grid()
plt.savefig('fracture_angles.png')


print(pairs_lkf,angles_lkf)

        #    plt.plot(seg[1,:]+0.5,seg[0,:]+0.5)
        #     plt.plot(seg2[1,:]+0.5,seg2[0,:]+0.5)