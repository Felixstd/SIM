import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
from skimage.morphology import skeletonize
import scipy.ndimage as ndim

def DoG_leads(in_array,max_kern,min_kern):
    """DoG: Difference of Gaussian Filters Combination as implemented in Linow & Dierking, 2017"""
    
    res = np.zeros(in_array.shape)
    c = np.arange(min_kern,max_kern+1)*0.5
    
    # for i in range(0,c.size-1):

    #     gaus1 = nan_gaussian_filter(in_array,c[i],truncate=2)
    #     gaus2 = nan_gaussian_filter(in_array,c[i+1],truncate=2)
    #     res += (gaus1 - gaus2)
    
    gaus1 = nan_gaussian_filter(in_array,c[0],truncate=2)
    gaus2 = nan_gaussian_filter(in_array,c[-1],truncate=2)
    
    res = (gaus1 - gaus2)
    return res

def nan_gaussian_filter(field,kernel,truncate):
    """ Version of scipy.ndimage.gaussian_filter that considers
    NaNs in the input array by setting them to zero and afterwards
    rescale the output array.
    Source https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    
    Input: field  - field to be filtered
           kernel - kernel of gaussian filter

    Output: gaussian_field - filtered field """
    
    field_nonnan = field.copy()
    mask_nan = np.ones(field.shape)

    field_nonnan[np.isnan(field)] = 0
    mask_nan[np.isnan(field)] = 0

    field_nonnan_f = ndim.gaussian_filter(field_nonnan,kernel,truncate=truncate)
    mask_nan_f = ndim.gaussian_filter(mask_nan,kernel,truncate=truncate)

    gaussian_field = field_nonnan_f/mask_nan_f

    #gaussian_field[np.isnan(field) | np.isnan(gaussian_field)] = 0.
    
    return gaussian_field

def hist_eq(array, number_bins=256):
    """ Histogram equalization
    Input:  array and number_bins (range of possible output valus: 0 to number_bins as integers)
    Output: histogram equalized version of array
    """
    # Compute histogram
    bins_center = np.linspace(np.nanmin(array[~np.isnan(array)]),np.nanmax(array[~np.isnan(array)]),number_bins)
    bins = np.append(bins_center-np.diff(bins_center)[0],(bins_center+np.diff(bins_center)[0])[-1])
    hist,bins = np.histogram(array[~np.isnan(array)].flatten(), bins)

    # Distribute bins equally to create lookup table
    new_values = np.floor((number_bins-1)*np.cumsum(hist/float(array[~np.isnan(array)].size)))

    # Compute equalized array with lookuptable
    array_equalized = np.take(new_values,np.digitize(array[~np.isnan(array)].flatten(),bins)-1)

    new_array_equalized = array.flatten()
    new_array_equalized[~np.isnan(new_array_equalized)]=array_equalized
	
    return new_array_equalized.reshape(array.shape)

def detect_lkf(eps_total):
    print(np.shape(eps_total))

    # eps_total = np.array(eps_total)
    
    proc_eps = np.log(eps_total)
    
    proc_eps[~np.isfinite(proc_eps)] = np.NaN
    
    proc_eps = hist_eq(proc_eps)
    lkf_detect = DoG_leads(proc_eps,5,1)
    
    lkf_detect = (lkf_detect > 5).astype('float')
    lkf_detect[~np.isfinite(proc_eps)] = np.NaN
    
    skeleton = skeletonize(lkf_detect)
    skeleton = np.array(skeleton, dtype = np.uint8)
    
    img_dilation = cv.dilate(skeleton, kernel = np.ones((3, 3), np.uint8) , iterations=1) 
    lines = cv.HoughLines(img_dilation, 1, np.pi / 180, 2)
    
    print(lines)

# Function to calculate the angle between two lines in degrees
    # def calculate_angle(line1, line2):
    #     rho1, theta1 = line1[0]
    #     rho2, theta2 = line2[0]

    #     # Convert angles to degrees and calculate the absolute difference
    #     angle1 = theta1 * 180 / np.pi
    #     angle2 = theta2 * 180 / np.pi
        
    #     return abs(angle1 - angle2)

    # # Filter pairs of lines that form a specific angle
    # desired_angle = 90  # You can change this to any desired angle
    # tolerance = 10  # Set an angle tolerance (e.g., +-10 degrees)
    # line_pairs = []

    # for i in range(len(lines)):
    #     for j in range(i + 1, len(lines)):
    #         angle_diff = calculate_angle(lines[i], lines[j])
        
    #     # If the angle difference is close to the desired angle (e.g., 90 degrees for perpendicular)
    #     if abs(angle_diff - desired_angle) < tolerance:
    #         line_pairs.append((lines[i], lines[j]))
            
    # print(line_pairs)

    plt.figure()
    ax = plt.axes()
    plt.pcolormesh(img_dilation)
    ax.set_aspect('equal', adjustable='box')
    # for i in range(0, len(lines)):
    #     rho = lines[i][0][0]
    #     theta = lines[i][0][1]
    #     a = np.cos(theta)
    #     b = np.sin(theta)
    #     x0 = a * rho
    #     y0 = b * rho
    #     pt1 = (int(x0 + 1000*(-b)), int(y0 + 1000*(a)))
    #     pt2 = (int(x0 - 1000*(-b)), int(y0 - 1000*(a)))
    #     print(pt1, pt2)
    #     plt.plot(pt1, pt2)
    #     plt.plot(pt2[0], pt2[1])
    
    plt.savefig('test.png', bbox_inches= 'tight')
    
def mohr_fracture_angle(mu): 
    """
    This function computes the fracture angle based on the 
    Mohr-Coulomb theory
    
    It assumes that the fracture angle is:
    
    2Theta = pi/2 - phi
    
    where phi is the internal angle of friction 
    
    and mu = tan(phi)

    Args:
        mu (array): the unique friction coefficients

    Returns:
        np.rad2deg(two_theta): 2theta in degrees
    """
    
    two_theta = np.pi/2 - np.arctan(mu)
    
    return np.rad2deg(two_theta)

def mohr_friction_coeffcient(angle): 
    """
    This function computes the friction coefficient 
    based on the mohr coulomb theory.
    
    


    Args:
        angle (array): in degrees

    Returns:
        mu: friction coefficient
    """
    
    mu = np.tan(np.pi/2 - np.deg2rad(angle))
    
    return mu

    
    
    
    
    
    
    
    
    


