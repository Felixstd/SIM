import numpy as np
from simplotting.data import read_data
from simplotting.LKF import lkf
from simplotting.LKF import lkf_hutter
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import scienceplots

plt.style.use('science')

dates = ['1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04', 
         '1990_01_01_00_00_05', '1990_01_01_00_00_06', '1990_01_01_00_00_07', '1990_01_01_00_00_08', 
         '1990_01_01_00_00_09', '1990_01_01_00_00_10', '1990_01_01_00_00_11', '1990_01_01_00_00_12', 
         '1990_01_01_00_00_13', '1990_01_01_00_00_14']


expno = '01'
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/aos/home/fstdenis/SIM/Experiments/MuPhi_GoodWay_Convergence/'

max_kernel=5
min_kernel=1
dog_thres=14
dis_thres=4
ellp_fac=3
angle_thres=35
eps_thres=0.5
lmin=4
max_ind=500 
use_eps=False
skeleton_kernel=0
k = -4
tests = 0

datadict = read_data.read_data(expno, dates, outputdir, MuPhi = True)

divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, zeta_tot, \
    uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()

shear_day = shear_tot[k][:501, 20:180]
div_day   = divergence_tot[k][:501, 20:180]
mu_day = muI_tot[k][:501, 20:180]
Y = np.arange(0, 501)
X = np.arange(20,180)

#find the range of mu for this day
mu_unique = np.unique(mu_day)
print(np.shape(mu_unique))

eps_tot = np.sqrt(shear_day**2 + div_day**2)

pairs_lkf_non, angles_lkf = lkf_hutter.lkf_detection(eps_tot, max_kernel, min_kernel, dog_thres, dis_thres, \
                ellp_fac, angle_thres, eps_thres, lmin, max_ind, use_eps, skeleton_kernel)

angles_lkf, idx_unique = np.unique(angles_lkf, return_index = True)
pairs_lkf_non = pairs_lkf_non[idx_unique[0]]


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


fig = plt.figure()
ax = plt.axes()
ax.set_aspect('equal', adjustable='box')
pc = plt.pcolormesh(shear_day,cmap = cmocean.cm.ice, norm = colors.Normalize(vmin=1e-4, vmax=1e-2))
fig.colorbar(pc, label = r'$\dot{\epsilon}_\mathrm{II}$ (1/s)')
for seg in pairs_lkf_non:
    # for seg in pairs:
        plt.plot(seg[1,:]+0.5,seg[0,:]+0.5, linewidth = 2, alpha = 0.5)
plt.title(r'$2\theta = {}$'.format(np.round(angles_lkf)))
plt.savefig('lkf_identif_{}.png'.format(dis_thres))


if tests:

    for dis_thres in range(0, 10):

        try : 
            pairs_lkf_non, angles_lkf = lkf_hutter.lkf_detection(eps_tot, max_kernel, min_kernel, dog_thres, dis_thres, \
                ellp_fac, angle_thres, eps_thres, lmin, max_ind, use_eps, skeleton_kernel)
        except:
            pairs_lkf_non = [np.array([[0],
        [0]]), np.array([[ 0],
        [0]])]
            angles_lkf =[0]

        try:
            angles_lkf, idx_unique = np.unique(angles_lkf, return_index = True)
            pairs_lkf_non = pairs_lkf[idx_unique[0]]
        except:
            pairs_lkf = pairs_lkf_non
        # print(pairs_lkf_non)
        # for pair in pairs_lkf_non:
        #     print(pair)
        #     for seg in pair:
        #         print(seg[1, :])
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


        fig = plt.figure()
        ax = plt.axes()
        ax.set_aspect('equal', adjustable='box')
        pc = plt.pcolormesh(shear_day,cmap = cmocean.cm.ice, norm = colors.Normalize(vmin=1e-4, vmax=1e-2))
        fig.colorbar(pc, label = r'$\dot{\epsilon}_\mathrm{II}$ (1/s)')
        try:
            for pairs in pairs_lkf_non:
                for seg in pairs:
                    plt.plot(seg[1,:]+0.5,seg[0,:]+0.5, linewidth = 2, alpha = 0.5)
        except:
            _ = 0
        plt.title(r'$2\theta = {}$'.format(np.round(angles_lkf)))
        plt.savefig('lkf_identif_{}.png'.format(dis_thres))


        print(pairs_lkf,angles_lkf)

        #    plt.plot(seg[1,:]+0.5,seg[0,:]+0.5)
        #     plt.plot(seg2[1,:]+0.5,seg2[0,:]+0.5)