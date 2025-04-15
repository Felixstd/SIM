import numpy as np
import matplotlib.pyplot as plt


def periodic_BC_test(filename, dir):
    
    var = np.loadtxt(dir+filename)

    diff_var_1 = var[:, -3] - var[:, 0]
    diff_var_2 = var[:, 1] -  var[:, -2]
    diff_var_3 = var[:, 2] -  var[:, -1]
    
    # print( var[:, 2])
    
    print(diff_var_1, diff_var_2,diff_var_3)

dir = '/storage/fstdenis/output_sim/'
file_1 = 'var1_01_1990_01_01_00_10_00.15'
file_2 = 'var2_01_1990_01_01_00_10_00.15'
file_3 = 'var2_04_04_1990_01_01_00_05_00.15'
file_4 = 'var1_03_1990_01_01_00_05_00.18'
file_5 ='var1_05_01_1990_01_01_07_00_00.40'
# periodic_BC_test(file_4, dir)
periodic_BC_test(file_5, dir)