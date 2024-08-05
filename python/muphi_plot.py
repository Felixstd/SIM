"""
This script is used to plot the output of the SIM model considering the 
mu(I) - Phi(I) rheology for pan-arctic simulations


"""


import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import datetime
import cmocean as cm

import warnings
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")


def phi_theo(I, c_phi) : 
    
    return 1 - c_phi*I


def load_params_muphi(filename, mask):
    """
    Function used to load the different variables and apply the 
    mask considering the right resolution


    Args:
        filename (str): filename of the variable and exp number
        mask (array): bool array for the Arctic mask

    Returns:
        _type_: _description_
    """

    #loading the array
    var = np.loadtxt(filename, dtype=None)

    #new variable for plotting and analysis
    varplot = np.zeros(var.shape)
    
    #applying the mask
    varplot[:,:] = var[:,:]
    varplot[mask == 0] = np.nan #putting nan on land 
    var_nonmasked = np.copy(varplot)
    varplot = np.ma.masked_invalid(varplot)

    return varplot, var_nonmasked


exp = "27"
exp_vp = '24'
dx = "20"
# dates = ["1990_01_03_00_20", "1990_01_03_20_00"]
# dates = ["1994_01_01_00_00", "1994_01_02_00_00", "1994_01_03_00_00"]

# dates = ["1990_01_02_00_00", "1990_01_03_00_00", "1990_01_07_00_00", "1990_01_14_00_00", "1990_01_21_00_00", "1990_01_28_00_00", "1990_01_30_00_00", "1990_02_01_00_00"]



dates = ["1990_01_14_00_00", "1990_02_14_00_00", "1990_03_14_00_00", "1990_04_14_00_00", "1990_05_14_00_00", "1990_06_01_00_00"]
# 1990-01-01:00:00:00  starting date
# 1990-02-01:00:20:00  end date
# 1990-01-02:00:00:00  posting date
# 1990-01-07:00:00:00  posting date
# 1990-01-14:00:00:00  posting date
# 1990-01-21:00:00:00  posting date
# 1990-01-28:00:00:00  posting date
# 1990-01-30:00:00:00  posting date
# 1990-02-01:00:00:00
outputdir = "/aos/home/fstdenis/SIM/output/"

popt_time = []


filemask="./MASKS/mask"+dx+"_1_0_1.dat"
mask = np.genfromtxt(filemask, dtype=None)

vp_comp = False
fit = False

for date in dates:
    fileI = outputdir+'I' + date + "." +exp
    file_phi = outputdir+'phi_I' + date + "." +exp
    file_mu = outputdir+'mu_I' + date + "." +exp
    file_conc = outputdir+'A' + date + "." +exp
    file_h = outputdir+'h'+ date + '.' + exp
    file_u = outputdir+'u'+ date + '.' + exp
    file_p = outputdir+'p'+ date + '.' + exp
    file_peq = outputdir+'Peq'+ date + '.' + exp
    file_shear = outputdir+'e_II'+ date + '.' + exp
    file_eta = outputdir+'eta'+ date + '.' + exp
    


    I_plot, I_nonmask = load_params_muphi(fileI, mask)
    phi_plot, _ = load_params_muphi(file_phi, mask)
    mu_plot, _ = load_params_muphi(file_mu, mask)
    conc_plot, conc_nonmask = load_params_muphi(file_conc, mask)
    
    h_plot,_ = load_params_muphi(file_conc, mask)
    shear_plot,_  = load_params_muphi(file_shear, mask)
    p_plot,_  = load_params_muphi(file_p, mask)
    peq,_  = load_params_muphi(file_peq, mask)


    datetime_format = '%Y_%m_%d_%H_%M'
    mydates = datetime.datetime.strptime(date, datetime_format)

    fileout="../Experiments/27/mu_I" + date +"_" + exp +".png"
    file2="test_A_I" + date +"_" + exp +".png"


    cmap = mpl.cm.get_cmap("inferno").copy()
    cmap.set_bad(color = 'gray', alpha = 1.)

    cmap_ice = cm.cm.ice.copy()
    # cmap_ice = mpl.cm.get_cmap("PuBu_r").copy()
    cmap_ice.set_bad(color = 'gray', alpha = 1.)

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex = True, sharey = True, figsize = (10, 8))


    pc = ax1.pcolormesh(I_plot, cmap=cmap,vmin = 0, vmax = 0.1)
    fig.colorbar(pc, ax = ax1, label = r'$I$')
    ax1.axis([0, mu_plot.shape[1], 0, mu_plot.shape[0]])

    pc = ax2.pcolormesh(mu_plot, cmap=cmap)
    fig.colorbar(pc,ax = ax2, label = r'$\mu(I)$')
    ax2.axis([0, phi_plot.shape[1], 0, phi_plot.shape[0]])

    phi_plot[phi_plot < 0 ] = 0
    pc = ax3.pcolormesh(phi_plot, cmap=cmap_ice,vmin = 0, vmax = 1)
    fig.colorbar(pc,ax = ax3, label = r'$\Phi(I)$')
    ax3.axis([0, mu_plot.shape[1], 0, mu_plot.shape[0]])

    pc = ax4.pcolormesh(conc_plot, cmap=cmap_ice, vmin = 0, vmax = 1)
    fig.colorbar(pc,ax = ax4, label = r'$A$')
    ax4.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
    
    pc = ax5.pcolormesh(shear_plot, cmap=cmap, vmin = 0, vmax = 1e-6)
    fig.colorbar(pc,ax = ax5, label = r'$\dot \epsilon_{II}$ (1/s)')
    ax5.axis([0, shear_plot.shape[1], 0, shear_plot.shape[0]])
    
    pc = ax6.pcolormesh(h_plot, cmap=cmap)
    fig.colorbar(pc,ax = ax6, label = r'$h$ (m)')
    ax6.axis([0, h_plot.shape[1], 0, h_plot.shape[0]])
    
    pc = ax7.pcolormesh(p_plot, cmap=cmap)
    fig.colorbar(pc,ax = ax7, label = r'$P$ (N/m)')
    ax7.axis([0, p_plot.shape[1], 0, p_plot.shape[0]])
    
    pc = ax8.pcolormesh(peq, cmap=cmap)
    fig.colorbar(pc,ax = ax8, label = r'$P_{eq}$ (N/m)')
    ax7.axis([0, peq.shape[1], 0, peq.shape[0]])

    fig.suptitle(mydates)
    plt.savefig(fileout, dpi = 500, bbox_inches = 'tight')

    
    if fit:
        I_nonmask = I_nonmask.flatten()
        conc_nonmask = conc_nonmask.flatten()
        
        idx_mask = ~np.isnan(conc_nonmask)

        masked_I = I_nonmask[idx_mask]
        masked_conc = conc_nonmask[idx_mask]
        
        idx_less1 = np.where(masked_I < 10)
        masked_I_cut = masked_I[idx_less1]
        masked_conc = masked_conc[idx_less1]
        log_I = np.log(masked_I_cut)

        
        # popt, pcov = curve_fit(phi_theo, log_I[~np.isinf(log_I)], masked_conc[~np.isinf(log_I)])
        popt, pcov = curve_fit(phi_theo, masked_I_cut, masked_conc)

        popt_time.append(popt)
        
        
        I = np.arange(1e-3, 10, 1e-5)
        plt.figure()
        plt.scatter(masked_I_cut, masked_conc, marker = 'o', color = 'k')
        plt.plot(I, 1 - 0.0824*I, color = 'r')
        plt.xlabel('I')
        plt.ylabel('A')
        plt.xscale('log')
        plt.title(mydates)
        plt.grid()
        plt.savefig(file2)
        
    
    if vp_comp:
        file_conc_vp = outputdir+'A' + date + "." +exp_vp
        conc_plot_vp, _ = load_params_muphi(file_conc_vp, mask)
    
        
        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize = (15, 6))
        
        pc = ax1.pcolormesh(phi_plot, cmap=cmap_ice, vmin = 0, vmax = 1)
        fig.colorbar(pc,ax = ax1, label = r'$\Phi(I)$')
        ax1.set_title('$\mu(I) - \Phi(I)$')
        ax1.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
        
        pc = ax2.pcolormesh(conc_plot, cmap=cmap_ice, vmin = 0, vmax = 1)
        fig.colorbar(pc,ax = ax2, label = r'$A$')
        ax1.set_title('$\mu(I) - \Phi(I)$')
        ax2.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
        
        pc = ax3.pcolormesh(conc_plot_vp, cmap=cmap_ice, vmin = 0, vmax = 1)
        fig.colorbar(pc,ax = ax3, label = r'$A$')
        ax3.set_title('VP')
        ax3.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
        
        pc = ax4.pcolormesh(conc_plot-phi_plot, cmap=plt.cm.inferno)
        fig.colorbar(pc,ax = ax4, label = r'$A - \Phi(I)$')
        ax4.set_title('$\mu(I) - \Phi(I)$')
        ax4.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
        
        pc = ax5.pcolormesh(conc_plot_vp-phi_plot, cmap=plt.cm.inferno)
        fig.colorbar(pc,ax = ax5, label = r'$A_{vp} - \Phi(I)$')
        ax5.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
        
        pc = ax6.pcolormesh(conc_plot_vp-conc_plot, cmap=plt.cm.inferno)
        fig.colorbar(pc,ax = ax6, label = r'$A_{vp} - A_{\mu}$')
        ax6.axis([0, conc_plot.shape[1], 0, conc_plot.shape[0]])
        
        fig.suptitle(mydates)
        plt.savefig('comp_SIC_phi_A_20km_'+date+exp+'.png', dpi=500, bbox_inches = 'tight')
    
    
    
    
print('Mean C_phi: ', np.mean(popt_time))

