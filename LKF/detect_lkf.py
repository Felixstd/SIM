import read_data
import cv2 as cv
import numpy as np
import matplotlib.pyplot as plt
import skimage.morphology 
import scipy.ndimage as ndim

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

def cut_neighbours(img):
    """Function that stencils each pixel with neighbouring pixels
    Input:  image (shape: MxN)
    Output: all neighbours (shape: MxNx3x3)
    """
    img = np.ascontiguousarray(img)  # won't make a copy if not needed
    X, Y = img.shape
    x, y = (1,1)
    overlap = 1
    shape = (((X-2*overlap)//x), ((Y-2*overlap)//y), x+2*overlap, y+2*overlap) # number of patches, patch_shape
    strides = img.itemsize*np.array([Y*x, y, Y, 1])
    return np.lib.stride_tricks.as_strided(img, shape=shape, strides=strides)

def nansum_neighbours(img):
    return np.nansum(cut_neighbours(img),axis=(2,3))

def elliptical_distance(seg_I,seg_II,ellp_fac=1,dis_thres=np.inf):
    """ Function to compute the elliptical distance between two
    segments, where the distance within the segment direction is
    weighted by 1 and the direction perpendicular to the direction
    by the factor 3. The weighted distance is computed from both
    segments and averaged. If the first computation already exceeds
    an threhold the second computation is skipped for efficiency.
    
    Input: seg_I - array with start and end coordinates of seg I
                   (rows dimension, column start-end)
           seg_II - array with start and end coordinates of seg II
           ellp_fac - weighting factor for ellipse
           dis_thres - distance threshold to stop computation

    Output: dis - elliptical distance"""

    # Determine basis vectors along seg_I direction
    e1 = (seg_I[:,0]-seg_I[:,1])
    e1 = (e1/np.sqrt(np.sum(e1**2))).reshape((2,1)) # Normalize basis vector

    e2 = np.dot(np.array([[0,-1],[1,0]]),e1)

    # Project connection vetor on basis vectors
    coeff = np.linalg.solve(np.hstack([e1,e2]),(seg_II[:,0] - seg_I[:,0]))

    if coeff[0]<0: coeff[0] = np.inf

    # Compute weighted distance
    d1 = np.sqrt(np.sum(coeff**2 * np.array([1,ellp_fac])))

    if d1 <= dis_thres:
        # Determine basis vectors along seg_II direction
        e1 = (seg_II[:,0]-seg_II[:,1])
        e1 = (e1/np.sqrt(np.sum(e1**2))).reshape((2,1)) # Normalize basis vector

        e2 = np.dot(np.array([[0,-1],[1,0]]),e1)

        # Project connection vetor on basis vectors
        coeff = np.linalg.solve(np.hstack([e1,e2]),(seg_I[:,0] - seg_II[:,0]))

        # Compute weighted distance
        d2 = np.sqrt(np.sum(coeff**2 * np.array([1,ellp_fac])))

        dis = 0.5*(d1+d2)
    else:
        dis = np.NaN

    return dis

def angle_segs(seg_I,seg_II):
    """ Function to compute the angle between two segments.
    
    Input: seg_I - array with start and end coordinates of seg I
                   (rows dimension, column start-end)
           seg_II - array with start and end coordinates of seg II

    Output: angle - angle between segments"""

    # Determine directions of segments
    e1 = (seg_I[:,0]-seg_I[:,1])
    e1 = (e1/np.sqrt(np.sum(e1**2))) # Normalize basis vector

    f1 = (seg_II[:,0]-seg_II[:,1])
    f1 = (f1/np.sqrt(np.sum(f1**2))) # Normalize basis vector

    # Determine angle between both directions
    angle = np.dot(e1,-f1)
    angle = np.arccos(angle)/np.pi*180

    return angle

def detect_segments(lkf_thin,eps_thres=0.1,max_ind=500):
    """ Function to detect segments of LKFs in thinned binary field
    The aim of this function is to split the binary field into 
    multiple smaller segments, and guarantee that all points in a
    segment belong to the same LKF. To do so a threshold for the
    deformation rate is establishes, which cuts to line that might
    belong to different LKFs. Note that also segments belonging
    to one LKF might be detected as multiple single segments in this 
    step.
    
    Input: lkf_thin  - thinned binary field
           eps_thres - deformation difference threshold to break a 
                       segment

    Output: seg_list - list of segments """


    # ------------------ Find starting points -----------------------
    seg = np.rollaxis(np.array(np.where((nansum_neighbours(lkf_thin)<=2) & (lkf_thin[1:-1,1:-1]==1))),1)
    seg = seg.reshape((seg.shape[0],seg.shape[1],1))
    # seg - array of dimension [N,2,M] with N being the number of segments
    #       M being an index for the point

    # Array of LKF points that have not been detected so far
    nodetect = lkf_thin[1:-1,1:-1].copy() 

    # Set all starting points to zero as they are detected already
    nodetect[(nansum_neighbours(lkf_thin)==2) & (nodetect==1)] = 0.
    nodetect_intm = np.zeros((nodetect.shape[0]+2,
                              nodetect.shape[1]+2))
    nodetect_intm[1:-1,1:-1] = nodetect.copy()

    # Initialize list of active segments
    active_detection = np.arange(seg.shape[0])

    # Deactivate segments that contain only one or two points
    deactivate_segs = np.where(nansum_neighbours(nodetect_intm)[seg[:,0].astype('int'),
                                                                seg[:,1].astype('int')].squeeze() != 1)
    if deactivate_segs[0].size > 0:
        active_detection = np.delete(active_detection,
                                     active_detection[deactivate_segs]) # remove from active list
    

    # --------------------- Detection loop --------------------------

    # Loop parameters
    num_nodetect = np.sum(nodetect) # Number of undetected pixels
    ind = 0 # Index of detection iteration
    max_ind = max_ind # Maximum number of iterations

    angle_point_thres = 5 # Number of last point in segment to compute the critical angel to break segments

    while num_nodetect > 0:
        #print ind, num_nodetect
        # Reduce segment array to active indeces
        seg_active = seg[active_detection]
    
        # Scheme of neighbouring cells
        #
        #   1 | 2 | 3
        #  -----------     ----> y 
        #   8 | X | 4      |
        #  -----------     v   
        #   7 | 6 | 5      x 
        #
    
        x = np.empty(seg_active.shape[:1])*np.NaN  
        y = np.empty(seg_active.shape[:1])*np.NaN

        for ix in [-1,0,1]:
            for iy in [-1,0,1]:
                indx = (seg_active[:,0,ind] + ix).astype('int')
                indy = (seg_active[:,1,ind] + iy).astype('int')
                mask = np.all([indx>=0,indx<nodetect.shape[0],
                               indy>=0, indy<nodetect.shape[1]], axis=0)
                nodetect_intm = np.zeros((nodetect.shape[0]+2,
                                          nodetect.shape[1]+2))
                nodetect_intm[1:-1,1:-1] = nodetect.copy()
                x[(nodetect_intm[indx+1,indy+1] == 1) & mask] = ix
                y[(nodetect_intm[indx+1,indy+1] == 1) & mask] = iy
        
        # Deactivate segments that ended
        deactivate_segs_end = np.where(np.all([np.isnan(x),np.isnan(y)],axis=0))[0]
    
        # Compute new points for segments
        seg_append = (seg_active[:,:,ind] + np.rollaxis(np.vstack([x,y]),1))#.astype('int')

        # Filter for sharp turns (degree larger 45deg)
        new_starts = np.empty((0,2))
        deactivate_segs_ang = np.empty((0,))
        if seg.shape[-1]>1:
            # Compute number of valid points per active segment
            num_points = np.sum(np.all(~np.isnan(seg_active),axis=1),axis=-1)
            # Limit points for the computation of the angle to threshold
            num_points[num_points>angle_point_thres] = angle_point_thres
            dx = (seg_active[:,:,-1]-seg_active[np.arange(seg_active.shape[0]),:,-num_points])/np.stack([num_points-1,num_points-1],axis=1) - (seg_append-seg_active[:,:,-1])
            new_starts = seg_append[np.sum(np.abs(dx),axis=1)>1] # high angle -> new starting point
            #print 'Number of segments broken by angel: %i' %np.sum(np.sum(np.abs(dx),axis=1)>1)

            deactivate_segs_ang = np.where(np.sum(np.abs(dx),axis=1)>1)[0]
    
        # Mark pixels as detected
        mask = np.all(~np.isnan(seg_append),axis=1) # masks all NaN entries in seg_append
        nodetect[seg_append[:,0][mask].astype('int'),seg_append[:,1][mask].astype('int')] = 0.
        if new_starts.size>0:
            nodetect[new_starts[:,0].astype('int'),new_starts[:,1].astype('int')] = 0
        nodetect_intm[1:-1,1:-1] = nodetect.copy()
    
        # Deactivate pixels with more than one neighbour and activate neighbours
        num_neighbours = nansum_neighbours(nodetect_intm)

        deactivate_segs_muln = np.where(num_neighbours[seg_append[:,0][mask].astype('int'),
                                                  seg_append[:,1][mask].astype('int')].squeeze() > 1)[0]
        deactivate_segs_muln = np.arange(seg_append.shape[0])[mask][deactivate_segs_muln]
        
        if (deactivate_segs_muln.size > 0):
            # Search for possibles neighbours to activate
            seg_deact = seg_append[deactivate_segs_muln]
            neigh_deactivate = np.vstack([seg_deact +
                                          i*np.hstack([np.ones((seg_deact.shape[0],1)),
                                                       np.zeros((seg_deact.shape[0],1))]) +
                                          j*np.hstack([np.zeros((seg_deact.shape[0],1)),
                                                       np.ones((seg_deact.shape[0],1))])
                                          for i in [-1,0,1] for j in [-1,0,1]])
            # new_starts_deact_ind = np.rollaxis(np.array(np.where((num_neighbours[neigh_deactivate[:,0].astype('int'),
            #                                                                      neigh_deactivate[:,1].astype('int')]<=2) & 
            #                                                      (nodetect[neigh_deactivate[:,0].astype('int'),
            #                                                                neigh_deactivate[:,1].astype('int')]==1))),1).squeeze()
            new_starts_deact_ind = np.rollaxis(np.array(np.where((num_neighbours[neigh_deactivate[:,0].astype('int'),
                                                                                 neigh_deactivate[:,1].astype('int')]<=2) & 
                                                                 (nodetect[neigh_deactivate[:,0].astype('int'),
                                                                           neigh_deactivate[:,1].astype('int')]==1))),1).squeeze()
            if new_starts_deact_ind.size>0:
                new_starts = np.append(new_starts,neigh_deactivate[new_starts_deact_ind].reshape((new_starts_deact_ind.size,2)),axis=0)
                nodetect[new_starts[:,0].astype('int'),new_starts[:,1].astype('int')] = 0
                nodetect_intm[1:-1,1:-1] = nodetect.copy()
                # if ind<5:
                #     print(new_starts.shape)
                #     print(new_starts)
                #     print(deactivate_segs_muln.shape)
                #     print(deactivate_segs_muln)

                

        # Test for segements that are on the same point
        nan_mask_segs = np.all(~np.isnan(seg_append),axis=-1)
        
        ravel_seg_append = np.ravel_multi_index((seg_append[nan_mask_segs,0].astype('int'),
                                                 seg_append[nan_mask_segs,1].astype('int')),
                                                lkf_thin[1:-1,1:-1].shape)
        seg_head_unique, seg_head_counts = np.unique(ravel_seg_append,return_counts=True)
        deactivate_segs_samehead = np.empty((0,))
        seg_head_continue = seg_head_unique[seg_head_counts==1]

        if np.any(seg_head_counts>1):
            deactivate_segs_samehead = np.hstack([np.where(ravel_seg_append==ihead)
                                                  for ihead in seg_head_unique[seg_head_counts>1]]).squeeze()
            new_starts = np.concatenate([new_starts,np.vstack(np.unravel_index(seg_head_unique[seg_head_counts>1],
                                                                               lkf_thin[1:-1,1:-1].shape)).T])
        #print(deactivate_segs_samehead)


        # Remove sharp turns from seg_append (here because search for new starting points
        # needs to run beforehand)
        if seg.shape[-1]>1:
            seg_append[np.sum(np.abs(dx),axis=1)>1,:] = np.NaN # Remove from appending list



        # Plot intermediate results
        if ind<5:#ind%5==0:
            do_plot = False
        else:
            do_plot = False
        
        if do_plot:
            fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(9,5),tight_layout=True)
            ax[0].pcolormesh(num_neighbours.copy())
            for i in range(seg.shape[0]):
                if np.any(active_detection==i):
                    col = 'r'
                else:
                    col = 'g'
                ax[0].plot(seg[i,1,:]+0.5,seg[i,0,:]+0.5,col)
                ax[0].text(seg[i,1,~np.isnan(seg[i,1,:])][-1]+0.5,seg[i,0,~np.isnan(seg[i,1,:])][-1]+0.5,'%i' %i,color='w')
            for i in range(neigh_deactivate.shape[0]):
                ax[0].plot(neigh_deactivate[i,1]+0.5,neigh_deactivate[i,0]+0.5,'m.')
            for i in range(new_starts.shape[0]):
                ax[0].plot(new_starts[i,1]+0.5,new_starts[i,0]+0.5,'c.')
            for i in range(active_detection.size):
                if np.any(deactivate_segs_end.copy()==i):
                    mark = 'x'
                elif np.any(deactivate_segs_ang.copy()==i):
                    mark = 'v'
                elif np.any(deactivate_segs_muln.copy()==i):
                    mark = 's'
                elif np.any(deactivate_segs_samehead.copy()==i):
                    mark = '>'
                else:
                    mark = '.'
                if ~np.isnan(seg_append[i,1]):
                    ax[0].plot(seg_append[i,1]+0.5,seg_append[i,0]+0.5,color='r',marker=mark)
                else:
                    ax[0].plot(seg[active_detection[i],1,-1]+0.5,seg[active_detection[i],0,-1]+0.5,color='r',marker=mark)
            

            #plt.figure()
            ax[1].pcolormesh(nodetect.copy()+lkf_thin[1:-1,1:-1])
            for i in range(seg.shape[0]):
                if np.any(active_detection==i):
                    col = 'r'
                else:
                    col = 'g'
                ax[1].plot(seg[i,1,:]+0.5,seg[i,0,:]+0.5,col)
                ax[1].text(seg[i,1,~np.isnan(seg[i,1,:])][-1]+0.5,seg[i,0,~np.isnan(seg[i,1,:])][-1]+0.5,'%i' %i,color='w')
            for i in range(active_detection.size):
                if np.any(deactivate_segs_end.copy()==i):
                    mark = 'x'
                elif np.any(deactivate_segs_ang.copy()==i):
                    mark = 'v'
                elif np.any(deactivate_segs_muln.copy()==i):
                    mark = 's'
                elif np.any(deactivate_segs_samehead.copy()==i):
                    mark = 'd'
                else:
                    mark = '.'
                if ~np.isnan(seg_append[i,1]):
                    ax[1].plot(seg_append[i,1]+0.5,seg_append[i,0]+0.5,color='r',marker=mark)
                else:
                    ax[1].plot(seg[active_detection[i],1,-1]+0.5,seg[active_detection[i],0,-1]+0.5,color='r',marker=mark)

            ax[0].set_xlim([380,395])
            ax[0].set_ylim([180,197])
            for iax in ax: iax.set_aspect('equal')


        
            
        # Test for multiple times same start
        new_starts_unique, new_starts_counts = np.unique(np.ravel_multi_index((new_starts[:,0].astype('int'),
                                                                               new_starts[:,1].astype('int')),
                                                                              lkf_thin[1:-1,1:-1].shape),
                                                         return_counts=True)
        
        new_starts_unique = np.array([i_seg_start for i_seg_start in new_starts_unique if not np.any(seg_head_unique==i_seg_start)],dtype='int')
        
        # if np.any(new_starts_counts > 1):
        #     # print 'Warning: %i starting points arises maximum %i-times' %(np.sum(new_starts_counts>1),
        #     #                                                               np.max(new_starts_counts))
        #     new_starts = np.vstack(np.unravel_index(new_starts_unique,lkf_thin[1:-1,1:-1].shape)).T
        new_starts = np.vstack(np.unravel_index(new_starts_unique,lkf_thin[1:-1,1:-1].shape)).T
                
        # Append new positions of this detection step
        num_new_starts = new_starts.shape[0]
        # Initialize list of new segment elements
        seg_append_time = np.empty((seg.shape[0],2))*np.NaN
        seg_append_time[active_detection] = seg_append
        seg_append_time = np.append(seg_append_time,new_starts,axis=0)
        seg_old_shape = seg.shape[0]
        # Fill up seg with NaNs for new starts
        seg = np.append(seg,np.empty((num_new_starts,2,seg.shape[-1]))*np.NaN,axis=0)
        # Append seg with new detected pixels
        seg = np.append(seg,seg_append_time.reshape(seg_append_time.shape[0],2,1),axis=-1)

        # Deactivate segments if finished
        active_detection_old = active_detection.copy()
        if np.any([(deactivate_segs_muln.size > 0),
                   (deactivate_segs_ang.size > 0),
                   (deactivate_segs_end.size > 0),
                   (deactivate_segs_samehead.size > 0)]):
            deactivate_segs = np.unique(np.append(deactivate_segs_muln,
                                                  np.append(deactivate_segs_ang,deactivate_segs_end)))
            deactivate_segs = np.unique(np.hstack([deactivate_segs_muln,deactivate_segs_ang,
                                                   deactivate_segs_end,deactivate_segs_samehead]))
            active_detection = np.delete(active_detection,deactivate_segs.astype('int')) # remove from active list

        # Activate new segments that started in this iteration
        active_detection = np.append(active_detection,np.arange(seg_old_shape, seg_old_shape + num_new_starts))

        # Compute number of undetected points and update ind
        num_nodetect = np.sum(nodetect)
        ind += 1
        
        if ind > max_ind:
            break
        if active_detection.size == 0:
            fac_starts = 100
            new_starts = np.append(np.rollaxis(np.array(np.where((nansum_neighbours(nodetect_intm)==3) & (nodetect==1))),1)[::fac_starts,:],
                                   np.rollaxis(np.array(np.where((nansum_neighbours(nodetect_intm)==3) & (nodetect==1))),1)[1::fac_starts,:],axis=0)
            
            # Mark new starts as detected
            if new_starts.size>0:
                nodetect[new_starts[:,0].astype('int'),new_starts[:,1].astype('int')] = 0
            nodetect_intm[1:-1,1:-1] = nodetect.copy()

            # Add new generated end points as well
            new_starts = np.append(new_starts,np.rollaxis(np.array(np.where((nansum_neighbours(nodetect_intm)<=2) & (nodetect==1))),1),axis=0)

            # Mark new starts as detected
            if new_starts.size>0:
                nodetect[new_starts[:,0].astype('int'),new_starts[:,1].astype('int')] = 0
            nodetect_intm[1:-1,1:-1] = nodetect.copy()

            # Test for multiple times same start
            new_starts_unique, new_starts_counts = np.unique(np.ravel_multi_index((new_starts[:,0].astype('int'),
                                                                                   new_starts[:,1].astype('int')),
                                                                                  lkf_thin[1:-1,1:-1].shape),
                                                             return_counts=True)
            if np.any(new_starts_counts > 1):
                print ('Warning: %i starting points arises maximum %i-times' %(np.sum(new_starts_counts>1),
                                                                              np.max(new_starts_counts)))

            # Append new positions of this detection step
            num_new_starts = new_starts.shape[0]
            seg_old_shape = seg.shape[0]
            # Fill up seg with NaNs for new starts
            seg = np.append(seg,np.empty((num_new_starts,2,seg.shape[-1]))*np.NaN,axis=0)
            # Fill in new start values 
            seg[seg_old_shape:,:,-1] = new_starts

            # Activate new segments that started in this iteration
            active_detection = np.append(active_detection,np.arange(seg_old_shape, seg_old_shape + num_new_starts))
            
            
            
            if active_detection.size == 0:
                break

    return seg

def compute_mn_eps(eps,seg):
    eps_mn = np.zeros(len(seg))
    
    for i in range(len(seg)):
        eps_mn[i] = np.mean(eps[1:-1,1:-1][seg[i][0,:],seg[i][1,:]])

    return eps_mn

def seg_reconnection(seg,segs,eps_segs,num_points_segs,dis_thres,angle_thres,eps_thres,ellp_fac=1):
    """ Function that does the reconnection
    
    Input: seg - list of segments
           segs - array with arrays containing coordinates 
                  of segments
                  (number segments, rows dimension, column start-end)
           eps_segs - mean deformation rate of segs
           num_points_segs - array of the  number of points of all 
                            segments
           angle_thres - angle threshold for reconnection
           ellp_fac - weighting factor for ellipse
           dis_thres - distance threshold for reconnection
           eps_thres - threshold difference in deformation rate 

    Output: seg - new list of reconnected segments"""
    

    # 1. Initialize probability matrix
    prob_ma = init_prob_matrix(segs,eps_segs,dis_thres,angle_thres,eps_thres,ellp_fac=ellp_fac)
    
    # 2. Loop over matrix and reconnect within one iteration the pair
    #    of segments that minimizes the probability matrix
    
    #    - Loop parameters
    ind = 0 
    num_pos_reconnect = np.sum(prob_ma[:,:,0]<1)
    max_ind = 500

    # loop (break criteria: no connection possible or max iterations are reached)
    while num_pos_reconnect >= 1:
        
        # 2.a. Find minimum of probability matrix
        ind_connect = np.unravel_index(np.nanargmin(prob_ma[:,:,0]),
                                       prob_ma[:,:,0].shape)
        ori_connect = prob_ma[ind_connect][1:]

        # 2.b. Update segments
        seg, segs, eps_segs, num_points_segs = update_segs(ind_connect,ori_connect,seg,segs,eps_segs,num_points_segs)
        
        # 2.c. Update probability matrix
        prob_ma = update_prob_matrix(prob_ma,ind_connect,segs,eps_segs,dis_thres,angle_thres,eps_thres,ellp_fac=ellp_fac)

        # 2.d. Update loop parameters
        ind += 1
        num_pos_reconnect = np.sum(prob_ma[:,:,0]<1)

        if ind>= max_ind:
            break

    return seg

def compute_prob(seg_I,segs,eps_segI,eps_segs,dis_thres,angle_thres,eps_thres,ellp_fac=1):
    """ Function to compute the probabilty for each segment in segs
    to be a succession of seg_I given the critical parameters for
    distance dis_thres, the angle angle_thres, and the deformation
    rate eps_thres.
    
    Input: seg_I - array coordinates of seg I
                   (rows dimension, column start-end)
           segs - array with arrays containing coordinates of other
                  segments
                  (number segments, rows dimension, column start-end)
           eps_segI - mean deformation rate of seg_I
           eps_seg - mean deformation rate of segs
           ellp_fac - weighting factor for ellipse
           dis_thres - distance threshold to stop computation
           eps_thres - threshold difference in deformation rate 

    Output: prob - probablility metric of segs, in second column
                   orientation information of segI and in third column
                   orientation information of the corresponding segment
                   from segs is stored in case of reconnection"""
    
    # 1. Check for similarity of deformation rates
    p_eps = np.abs(eps_segs-eps_segI)/eps_thres
    p_eps[p_eps > 1] = np.nan

    mask_eps = ~np.isnan(p_eps)
    segs_i = segs[mask_eps]
    
    
    # 2. Find corresponding starting and end points and first instance of 
    #    distance thresholding
    ori_segI, ori_segs, mask_dis_i = find_pos_connect(seg_I,segs_i,dis_thres)
     
    segs_i = segs_i[mask_dis_i]
    ori_segI = ori_segI[mask_dis_i]
    ori_segs = ori_segs[mask_dis_i]


    # 3. Check angle between segments and angle thresholding
    p_ang = np.zeros(p_eps.shape) * np.nan
    mask_ang = np.zeros(segs_i.shape[0]).astype('bool')
    
    for i in range(segs_i.shape[0]):
        # Resort arrays if necessary to have proper orientation with starting 
        # points
        if ori_segI[i] == 1:
            seg_I_i = seg_I[:,::-1].copy()
        else:
            seg_I_i = seg_I.copy()
        if ori_segs[i] == 1:
            seg_II_i = segs_i[i][:,::-1].copy()
        else:
            seg_II_i = segs_i[i].copy()

        # Determine angle
        p_ang[np.arange(p_ang.size)[mask_eps][mask_dis_i][i]] = angle_segs(seg_I_i,seg_II_i)/angle_thres
        mask_ang[i] = (p_ang[mask_eps][mask_dis_i][i]<=1)

    p_ang[p_ang>1] = np.nan
    
    segs_i = segs_i[mask_ang]
    ori_segI = ori_segI[mask_ang]
    ori_segs = ori_segs[mask_ang]

        

    # 4. Compute elliptical distance and final distance thresholding
    p_dis = np.zeros(p_eps.shape) * np.nan
    
    for i in range(segs_i.shape[0]):
        # Resort arrays if necessary to have proper orientation with starting 
        # points
        if ori_segI[i] == 1:
            seg_I_i = seg_I[:,::-1].copy()
        else:
            seg_I_i = seg_I.copy()
        if ori_segs[i] == 1:
            seg_II_i = segs_i[i][:,::-1].copy()
        else:
            seg_II_i = segs_i[i].copy()

        # Determine distance
        p_dis[np.arange(p_dis.size)[mask_eps][mask_dis_i][mask_ang][i]] = elliptical_distance(seg_I_i,seg_II_i,ellp_fac=ellp_fac,dis_thres=dis_thres)/dis_thres

    p_dis[p_dis>1] = np.nan


    # 5. Compute joint probability as sum of all three components

    prob = np.sqrt(p_eps**2 + p_ang**2 + p_dis**2)


    # 6. Save orientation of the corresponding connection partners
    ori_segI_all = np.zeros(prob.size) * np.nan
    ori_segs_all = np.zeros(prob.size) * np.nan
    ori_segI_all[np.arange(p_dis.size)[mask_eps][mask_dis_i][mask_ang]] = ori_segI
    ori_segs_all[np.arange(p_dis.size)[mask_eps][mask_dis_i][mask_ang]] = ori_segs

    return np.rollaxis(np.stack([prob,ori_segI_all,ori_segs_all]),1)

def find_pos_connect(seg_I,segs,dis_thres):
    """ Function to determine the possible connection segments
    and to compute both corresponding starting point. The latter
    information is given as arrays of the orientation where 1 means
    that the current orientation has the starting point in the
    first column and -1 that the starting point is in the second
    column. These orientation can be used by indexing to flip the
    the array to the right order [:,::i] where i=1,-1
    
    Input: seg_I - array coordinates of seg I
                   (rows dimension, column start-end)
           segs - array with arrays containing coordinates of other
                  segments
                  (number segments, rows dimension, column start-end)

    Output: ori_segI - required orientation of seg_I
            ori_segs - required orientation of segs 
            mask - mask of all segments in segs that fulfill distance
                   criteria"""

    # Compute displacement from starting and end points in segs from start in seg_I
    disp_start = np.rollaxis(np.rollaxis(segs,-1,start=1)-seg_I[:,0],-1,start=1)

    # Compute displacement from starting and end points in segs from end in seg_I
    disp_end   = np.rollaxis(np.rollaxis(segs,-1,start=1)-seg_I[:,1],-1,start=1)

    # # Filter for larger displacements than dis_thres
    # mask = np.all([np.all(np.any(np.abs(disp_start)>dis_thres,axis=1),axis=1),
    #                np.all(np.any(np.abs(disp_end  )>dis_thres,axis=1),axis=1)],axis=0)
    # disp_start[mask,:,:] = np.NaN
    # disp_end[mask,:,:]   = np.NaN

    # Compute distance only for filtered displacements
    dis = np.hstack([np.sqrt(np.sum(disp_start**2,axis=1)).reshape((segs.shape[0],1,2)),
                     np.sqrt(np.sum(disp_end**2,  axis=1)).reshape((segs.shape[0],1,2))])

    # Give out combination of orientation of segments with starting point being first column
    ori_segI = np.argmin(dis,axis=1)
    ori_segs = np.argmin(np.min(dis,axis=1),axis=1)
    ori_segI = ori_segI[np.arange(ori_segs.size),ori_segs]

    # Filter for larger displacements than dis_thres
    mask = np.all(np.abs(dis)>dis_thres,axis=(1,2))

    return ori_segI, ori_segs, ~mask

def init_prob_matrix(segs,eps_segs,dis_thres,angle_thres,eps_thres,ellp_fac=1):
    """ Function to initialize the probability matrix given the 
    probability of all possible combinations of segments to belong
    to the same deformation feature. The probabilty matrix is a
    upper triangular matrix with empty diagonal.
    
    Input: segs - array with arrays containing coordinates of
                  segments
                  (number segments, rows dimension, column start-end)
           eps_seg - mean deformation rate of segs
           ellp_fac - weighting factor for ellipse
           dis_thres - distance threshold to stop computation
           eps_thres - threshold difference in deformation rate 

    Output: prob_ma - probablility matrics of segs"""

    
    # 1. Initialize empty probability matrix
    num_segs = segs.shape[0]
    prob_ma = np.zeros((num_segs,num_segs,3)) * np.nan

    # 2. Loop over all segments an fill
    for i_s in range(num_segs-1):
        prob_ma[i_s,i_s+1:,:] = compute_prob(segs[i_s],segs[i_s+1:],eps_segs[i_s],eps_segs[i_s+1:],dis_thres,angle_thres,eps_thres,ellp_fac=ellp_fac)
        
    return prob_ma

def update_segs(ind_connect,ori_connect,seg,segs,eps_segs,num_points_segs):
    """ Function to update the list of segment seg, array of start
    and end points segs, and the array of mean deformation rates.

    Input: ind_connect - index that were connected in this step
           ori_connect - orientation of segmeents that are reconnected
                         ( 0 means orientation as given in segs is
                         right, if 1 segment needs to be reversed)
           seg - list of segments
           segs - array with arrays containing coordinates 
                  of segments
                  (number segments, rows dimension, column start-end)
           eps_segs - mean deformation rate of segs
           num_points_segs - array of the  number of points of all 
                            segments

    Output: seg - updated list of segments
            segs - updated array with arrays containing coordinates 
                      of segments
                      (number segments, rows dimension, column start-end)
            eps_segs - updated mean deformation rate of segs_up
            num_points_segs - updated array of the  number of points
                                of all segments"""
    
    # 1. Update list of segments seg
    #    - Update smaller index element
    seg[ind_connect[0]] = np.append(seg[ind_connect[0]][:,::int(2*ori_connect[0]-1)],
                                    seg[ind_connect[1]][:,::int(-2*ori_connect[1]+1)],
                                    axis=1)
    #    - Remove larger index element
    seg.pop(ind_connect[1]);
    

    # 2. Update array of end and starting points
    #    - Update smaller index element
    segs[ind_connect[0]] = np.stack([segs[ind_connect[0]][:,::int(2*ori_connect[0]-1)][:,0],
                                     segs[ind_connect[1]][:,::int(-2*ori_connect[1]+1)][:,-1]]).T
    #    - Remove larger index element
    segs = np.delete(segs,(ind_connect[1]),axis=0)


    # 3. Update array of mean deformation rates
    #    - Update smaller index element
    eps_segs[ind_connect[0]] = (((eps_segs[ind_connect[0]]*num_points_segs[ind_connect[0]]) +
                                 (eps_segs[ind_connect[1]]*num_points_segs[ind_connect[1]]))/
                                (num_points_segs[ind_connect[0]]+num_points_segs[ind_connect[1]]))
    #    - Remove larger index element
    eps_segs = np.delete(eps_segs,(ind_connect[1]),axis=0)


    # 4. Update array of number of points of all segments
    #    - Update smaller index element
    num_points_segs[ind_connect[0]] = (num_points_segs[ind_connect[0]] + 
                                       num_points_segs[ind_connect[1]])
    #    - Remove larger index element
    num_points_segs = np.delete(num_points_segs,(ind_connect[1]),axis=0)


    return seg, segs, eps_segs, num_points_segs


def update_prob_matrix(prob_ma,ind_connect,segs_up,eps_segs_up,dis_thres,angle_thres,eps_thres,ellp_fac=1):
    """ Function to update the probability matrix given the 
    probability of all possible combinations of segments to belong
    to the same deformation feature. Only the rows and columns
    corresponding to indeces ind_connect are updated as the others
    remain unchanged. The column and row corresponding to the larger 
    index ind_connect[1] are removed from the matrix and the others
    are recalculated. As the orientation and deformation rate of the
    newly reconnected segment besides its length might have changed 
    a new computation of the entire row is required instead of only
    updating all non NaN values.
    
    Input: prob_ma - probablility matrics of segs
           ind_connect - index that were connected in this step
           segs_up - updated array with arrays containing coordinates 
                     of segments
                     (number segments, rows dimension, column start-end)
           eps_seg_up - updated mean deformation rate of segs_up
           ellp_fac - weighting factor for ellipse
           dis_thres - distance threshold to stop computation
           eps_thres - threshold difference in deformation rate 

    Output: prob_ma_up - probablility matrics of segs"""

    # 1. Remove column and row corresponding to the larger index
    prob_ma = np.delete(np.delete(prob_ma,ind_connect[1],axis=0),
                        ind_connect[1],axis=1)

    # 2. Reevaluate the probabilty in the row for the lower index
    i_s = ind_connect[0]
    prob_ma[i_s,i_s+1:,:] = compute_prob(segs_up[i_s],segs_up[i_s+1:],eps_segs_up[i_s],eps_segs_up[i_s+1:],dis_thres,angle_thres,eps_thres,ellp_fac=ellp_fac)
    prob_ma[:i_s,i_s,:] = compute_prob(segs_up[i_s],segs_up[:i_s],eps_segs_up[i_s],eps_segs_up[:i_s],dis_thres,angle_thres,eps_thres,ellp_fac=ellp_fac)[:,[0,2,1]]

    return prob_ma

def filter_segs_lmin(seg,lmin):
    """ Function to filter all segements in seg where the distance
    between start and end point is below threshold lmin"""
    return [i for i in seg if np.sqrt(np.sum((i[:,0]-i[:,-1])**2))>=lmin]

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
                           
                           
dates = ['1990_01_01_00_00_01', '1990_01_01_00_00_02', '1990_01_01_00_00_03', '1990_01_01_00_00_04', 
         '1990_01_01_00_00_05', '1990_01_01_00_00_06', '1990_01_01_00_00_07', '1990_01_01_00_00_08', 
         '1990_01_01_00_00_09', '1990_01_01_00_00_10', '1990_01_01_00_00_11', '1990_01_01_00_00_12', 
         '1990_01_01_00_00_13', '1990_01_01_00_00_14']


expno = '01'
outputdir = "/storage/fstdenis/output_sim/"

datadict = read_data.read_data(expno, dates, outputdir, MuPhi = True)


divergence_tot, shear_tot, h_tot, A_tot, p_tot, sigI_tot, sigII_tot, muI_tot, phi_tot, I_tot = \
        datadict.values()
        
shear_day =shear_tot[-1]
div_day   = divergence_tot[-1]

eps_tot = np.sqrt(shear_day**2 + div_day**2)


proc_eps = np.log(eps_tot)
proc_eps[~np.isfinite(proc_eps)] = np.NaN

proc_eps = hist_eq(proc_eps)

lkf_detect = DoG_leads(proc_eps,5,1)
### Filter for DoG>0
lkf_detect = (lkf_detect > 15).astype('float')
lkf_detect[~np.isfinite(proc_eps)] = np.NaN

lkf_thin =  skimage.morphology.skeletonize(lkf_detect).astype('float')


seg_f = detect_segments(lkf_thin,max_ind=500) # Returns matrix fill up with NaNs
## Convert matrix to list with arrays containing indexes of points
seg = [seg_f[i][:,~np.any(np.isnan(seg_f[i]),axis=0)].astype('int')
        for i in range(seg_f.shape[0])]
# ## Apply inter junction connection
# seg = connect_inter_junctions(seg,lkf_thin)
## Filter segments that are only points
seg = [i for i in seg if i.size>2]

eps_mn = compute_mn_eps(np.log10(eps_tot),seg)
num_points_segs = np.array([i.size/2. for i in seg])
## Initialize array containing start and end point of segments
segs = np.array([np.stack([i[:,0],i[:,-1]]).T for i in seg])
seg = seg_reconnection(seg,segs,eps_mn,num_points_segs,dis_thres,
                           angle_thres,eps_thres,ellp_fac=ellp_fac)

    # Filter too short segments
seg = filter_segs_lmin(seg,lmin)

    # Convert to indexes of the original input image
segs = [segi+1 for segi in seg]
print(seg)
for seg in segs:
        plt.plot(seg[1,:]+0.5,seg[0,:]+0.5,'k')
plt.savefig('./test_lkf.png', dpi = 500)

# shear_day *= 255.0/shear_day.max() 
# shear_day = np.uint8(shear_day)


# print(type(shear_day))

# # dst = cv.cvtColor(shear_day, cv.COLOR_BGR2GRAY)
# dst = cv.Canny(shear_day, 0, 50, None, 3)
# lines = cv.HoughLines(dst, 1, np.pi / 180, 150, None, 0, 0)

# if lines is not None:
#     for i in range(0, len(lines)):
#         rho = lines[i][0][0]
#         theta = lines[i][0][1]
#         a = np.cos(theta)
#         b = np.sin(theta)
#         x0 = a * rho
#         y0 = b * rho
#         pt1 = (int(x0 + 1000*(-b)), int(y0 + 1000*(a)))
#         pt2 = (int(x0 - 1000*(-b)), int(y0 - 1000*(a)))

#         cv.line(dst, pt1, pt2, (0,0,255), 3, cv.LINE_AA)
# cv.imwrite('test_lkf.png',shear_day) 
# # plt.figure()
# # plt.imshow(lines)
# # plt.savefig('./test_lkf.png', dpi = 500)