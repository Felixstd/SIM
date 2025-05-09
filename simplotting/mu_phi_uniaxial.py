import numpy as np
from simplotting.data import read_data
from simplotting.plotting import plot
from simplotting.analysis import analysis
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import warnings
import argparse
import matplotlib.lines as mlines

warnings.filterwarnings("ignore")
 
parser = argparse.ArgumentParser(description="For analysis",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-t", "--time", type=int, help = "1 for min, 2 for 10 min, 3 for 2hours and 4 for 10 hours simulation")
parser.add_argument("-e", "--expno", type=int, help = "Experiment number ")
parser.add_argument("-mu0", "--mu0", type=float, help = "Lower bound for mu ")
parser.add_argument("-muinf", "--muinfty", type=float, help = "higher bound for mu ")
args = vars(parser.parse_args())

time = args['time']
expno = args['expno']
mu_0 = args['mu0']
mu_infty = args['muinfty']

if time == 1:
        
        start = 1
        start_k = 61
        #--- Time for 1 minute run ---#
        start = datetime(1990, 1, 1, 0, 1, start)
        intervals = [timedelta(seconds=1)] * 58
        
elif time == 2:
        #--- Time for 10 minutes run ---#
        start = datetime(1990, 1, 1, 0, 0, 30)
        start_k =1
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)] * int((9*60/30))
        
elif time == 3:
        # --- Time for 2 hours 30 minutes run ---#
        start = datetime(1990, 1, 1, 0, 5, 00)
        start_k = 6
        # Time interval (30 seconds initially, then 5-minute steps)
        intervals = [timedelta(seconds=30)]*5  + [timedelta(minutes=4)]*0 + [timedelta(minutes=5)]*24 #+ [timedelta(minutes=3)]*1 + [timedelta(minutes=1)]*1
       
elif time == 4:
        #--- Time for 2 hours 30 minutes run ---#
        start = datetime(1990, 1, 1,23,0, 00)
        start_k =138
        # start = datetime(1990, 1, 1, 0 ,10, 00)
        # start_k = 1
        intervals = [timedelta(minutes=10)]*5#+ [timedelta(minutes=3)]*1 + [timedelta(minutes=1)]*1
        # intervals = [timedelta(minutes=10)]*142 #+ [timedelta(minutes=3)]*1 + [timedelta(minutes=1)]*1

dates = [(start + sum(intervals[:i], timedelta())).strftime('%Y_%m_%d_%H_%M_%S') for i in range(len(intervals)+1)]
print(dates)

#------- MuPhi -------#
#---- Uniaxial Tests ----#
outputdir = "/storage/fstdenis/output_sim_MuPhi_Runs_Tests_Dilat/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_Runs_Tests_Dilat/'

#---- Shear Tests ----#
outputdir = "/storage/fstdenis/output_sim_MuPhi_Runs_ShearExperiments/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuPhi_Runs_ShearExperiments/'


#------- Mu -------#
outputdir = "/storage/fstdenis/output_sim/output_mu_ShearExperiments/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/Mu_ShearExperiments/'

#------- Mu Periodic Boundary Conditions -------#
outputdir = "/storage/fstdenis/output_sim_mu_PeriodicBC/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/Mu_PeriodicBC/'

#------- Mu Periodic Boundary Conditions -------#
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/VP_PeriodicBC/'

#------- Mu Parameter Sensibility -------#
outputdir = "/storage/fstdenis/output_sim/"
figdir = '/storage/fstdenis/Experiments_Results_MuPhi/MuCompressibleParameters/'

markers = ['.', 'x', '^', '*']
fig_ellipse = plt.figure()
ax_ellipse = plt.axes()

ax_ellipse.set_aspect('equal')
ax_ellipse.grid()
expnos = [expno]
# expnos = [1, 2, 3, 4]
# expnos = [5, 6, 7, 8]
# mub = [1, 2, 5, 10] 
mub = [1]
# for i in range(expno, expno+1):


fig = plt.figure()
ax = plt.axes()

for i, expno in enumerate(expnos):
    expno = "{:02d}".format(expno)

    if not os.path.isdir(figdir+expno):
        os.mkdir(figdir+expno)

    print('Analysising Exp: ', expno)

    dx = 1e3*10
    dy = dx
    dt = 10

    angle_phi = 20*np.pi/180

    Ny = 502
    Nx = 202

    # if (Ny == 252 and Nx == 102):
    mask = np.ones((Ny, Nx))
    mask[0, :] = 0
    
    filemask=outputdir+"mask.dat"
    maskC = np.genfromtxt(filemask, dtype=None)
    
    N_transect = 51
    muphi = 1
    log = 0
    time = np.arange(1, len(intervals)+2, 1)*dt

    #_--------- READING DATA ----------#
    datadict = read_data.read_data(expno, start_k, dates, outputdir, MuPhi = muphi)
    if muphi:
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot, sigI_tot, sigII_tot, \
            zeta_tot, eta_tot, uair_tot, vair_tot, muI_tot, phi_tot, I_tot, shearI_tot, Pmax_tot, Peq_tot = datadict.values()
            
        I = I_tot[0]
        shearI = shearI_tot[0]
        mu = muI_tot[0]
        zeta = zeta_tot[0]
        eta = eta_tot[0]
        print('shearI', np.vstack([shearI[:, 1],shearI[:, -1]]).T)
    
    else: 
        
        divergence_tot, shear_tot, h_tot, A_tot, p_tot, u_tot, v_tot,sigI_tot, sigII_tot, zeta_tot, eta_tot, uair_tot, vair_tot = \
            datadict.values()
            
    # datadict_A = read_data.read_data_var('A', expno, start_k, dates, outputdir)
    # A_tot = datadict_A['variable_date']
    # print(A_tot)
    # print(datadict)
    # #---------- Analysing Wind Forcing----------#
    # analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', muphi)

    # # #---------- Analysing Invariants ----------#
    # analysis.invariants(dates, expno, datadict, Ny, Nx, dx, maskC, figdir, mu_0, mu_infty)

    # #---------- Plotting ----------#
    plot.uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)
    
    #---- Dilatation ----#
    # dilatation_time = []
    # time = np.arange(0, len(dates))*10/60
    # for c in range(len(A_tot)) : 
    #     A = A_tot[c] 
    #     dilatation_x =  np.cumsum(1-A, axis = 1)
    #     dilatation_y =  np.cumsum(1-A, axis = 0)
    #     dilat = np.nanmean(1-A)
    #     dilat = np.sum(1-A)/(Nx*Ny)
    #     dilatation_time.append(dilat)
    
    
    # ax.plot(time, dilatation_time, label = r'$\mu_b = {}$'.format(mub[i]))
    
    
    # j = -1
    # sigI = sigI_tot[j]
    # sigII = sigII_tot[j]
    
    # sigI_1 = sigI[0:90, :]
    # sigI_2 = sigI[90:110, :]
    # sigI_3 = sigI[110:, :]
    
    # sigII_1 = sigII[0:90, :]
    # sigII_2 = sigII[90:110, :]
    # sigII_3 = sigII[110:, :]
    
    # ax_ellipse.scatter(sigI_1.flatten(), sigII_1.flatten(), alpha = 0.2,color = 'r', marker = markers[i], s=2)
    # ax_ellipse.scatter(sigI_3.flatten(), sigII_3.flatten(),  alpha = 0.2,color = 'r', marker = markers[i], s=2)
    # ax_ellipse.scatter(sigI_2.flatten(), sigII_2.flatten(),  alpha = 0.2,color = 'g' ,marker = markers[i], s=2)
    
#     mean_mu_time, mean_div_time, mean_shear_time = analysis.mean_values(dates, datadict, dx, dt)
#     plot.plot_mean(mean_mu_time, time, r'$\langle\mu\rangle$', 'mean_mu.png', expno, figdir)
#     plot.plot_mean(mean_div_time, time, r'$\langle\ \dot{\epsilon}_I \rangle$ (day$^{-1}$)', 'mean_div.png', expno, figdir)
#     plot.plot_mean(mean_shear_time, time, r'$\langle\ \dot{\epsilon}_\mathrm{II} \rangle$ (day$^{-1}$)', 'mean_shear.png', expno, figdir)
    
    # analysis.velocity_transects(datadict, dates, dx,figdir, expno, MuPhi = muphi)
#     
    # plot.totdef_uniaxial(dates, expno, datadict, dx, figdir, mu_0, mu_infty, angle_phi, MuPhi = muphi, log = log)

ax.set_ylabel(r'$\overline{1-A}$')
ax.set_xlabel(r'Time (hr)')
ax.grid()
                        
fig.legend(loc='upper center', 
        labelcolor='linecolor',  bbox_to_anchor=(1.05, 0.9))
plt.savefig('dilatation.png')

    
u  =u_tot[-1]
print(np.shape(u))
# u_split = np.hsplit(u, 2)
y = np.shape(u)[1] // 2  # 250
u_split_1 = u[:, :y]
u_split_2 = u[:, y+1:]

# print(np.shape(u_split))
plt.figure()
plt.pcolor(u_split_1 - u_split_2)
plt.colorbar()
plt.savefig('error_vp.png')




sigI = sigI_tot[-1]
sigII = sigII_tot[-1]
print(np.shape(sigI))

sigI_t = np.linspace(-1 ,0, 2000 )
# sigI_t = np.arange(0,1, 0.0001 )
sigII_t = (1/2)*np.sqrt((1/2)**2-(sigI_t+1/2)**2)
colors = ['blue', 'green', 'red']


    
ext = mlines.Line2D([], [], color='r', marker= 'None', linestyle='None',
                      markersize=1, label='Extremities')
cent = mlines.Line2D([], [], color='g', marker= 'None', linestyle='None',
                      markersize=1, label='Center')
zerodeg = mlines.Line2D([], [], color='k', marker= '.', linestyle='None',
                      markersize=1, label='0°')
tendeg = mlines.Line2D([], [], color='k', marker= 'x', linestyle='None',
                      markersize=1, label='10°')
twentydeg = mlines.Line2D([], [], color='k', marker= '^', linestyle='None',
                      markersize=1, label='20°')
thirtydeg = mlines.Line2D([], [], color='k', marker= '*', linestyle='None',
                      markersize=1, label='30°')



                        
fig_ellipse.legend(loc='upper center', 
        handles=[ext, cent, zerodeg, tendeg, twentydeg, thirtydeg], labelcolor='linecolor',  bbox_to_anchor=(1.05, 0.9))
    
    
ax_ellipse.plot(sigI_t, sigII_t, color = 'k')
ax_ellipse.plot(sigI_t, -sigII_t, color = 'k')
ax_ellipse.set_xlabel(r'$\sigma_I/P$')
ax_ellipse.set_ylabel(r'$\sigma_{II}/P$')
fig_ellipse.savefig('sigI_sigII.png', dpi=500)


            # h = h_tot[0]
    # u = u_tot[0]
    # v = v_tot[0]
    # p = p_tot[0]
    # # mu = muI_tot[1]
    # A = A_tot[0]
    # div = divergence_tot[0]
    # shear = shear_tot[0]
    # zeta = zeta_tot[0]
    # eta = eta_tot[0]
        # print(np.shape(shearI))
        # print('I', np.vstack([I[100:, 1],I[100:, -1]]).T)
        # print('shearI', shearI[:, -1])
        # print('shearI', shearI[:, 1])

        # with np.printoptions(threshold=np.inf):
            # print('shearI', np.where(shearI == -999)[0])
        # print(np.where(shearI == -999)[0], shearI[0, -1])
        # print('shearI', np.vstack([shearI[201, 0],shearI[201, 501]]).T)
        # print('mu', np.vstack([mu[100:, 1],mu[100:, -1]]).T)
        # print('zeta', np.vstack([zeta[100:, 1],zeta[100:, -1]]).T)
        # print('eta', np.vstack([eta[100:, 1],eta[100:, -1]]).T)
    # print(np.shape(u))
    # print('h', np.vstack([h[100:, 1],h[100:, -1]]).T)
    # print('A', np.vstack([A[100:, 1],A[100:, -1]]).T)
    # print('u', np.vstack([u[100:, 0],u[100:, -1]]).T)
    # print('v', np.vstack([v[100:, 0],v[100:, -1]]).T)
    # print('p', np.vstack([p[100:, 0],p[100:, -1]]).T)
    # print('div',np.vstack([div[100:, 1],div[100:, -1]]).T)
    # print('shear',np.vstack([shear[100:, 1],shear[100:, -1]]).T)
    # print('shear',np.vstack([shear[100:, 0],shear[100:, -2]]).T)
    # print('zeta', np.vstack([zeta[100:, 1],zeta[100:, -1]]).T)
    # print('eta', np.vstack([eta[100:, 1],eta[100:, -1]]).T)