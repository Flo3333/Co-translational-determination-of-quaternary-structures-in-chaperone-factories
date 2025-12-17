import numpy as np
import bigfish.stack as stack
import bigfish.segmentation as seg
import bigfish.plot as plot
from pbwrap.segmentation import random_walker_segmentation
from scipy.ndimage import binary_dilation


def create_3D_mask(image_stack: np.ndarray, log_kernel_size= 2, threshold= 40, path_out = None) :
    """
    Convert stack (3D) of Beta-Cadhenin (EGFP) to 3D mask.
    Segmentation is perfomed plane by plane, any croping should be done prior to this operation

    Returns
    -------
        res_mask : np.ndarray,
            3D mask array.
    """

    stack.check_array(image_stack, ndim= 3)
    cad_mask = random_walker_segmentation(stack.maximum_projection(image_stack), 25,90, beta= 0.8)
    cad_mask_cleaned = seg.clean_segmentation(cad_mask, 10000)
    del cad_mask


    image_stack_log = stack.log_filter(image_stack, sigma= log_kernel_size)
    image_stack_log[:, ~cad_mask_cleaned] = 0

    res_mask = np.array([seg.thresholding(image_log, threshold) for image_log in image_stack_log]).astype(bool)

    if path_out != None :
        for slice_idx in range(len(res_mask)) :
            plot.plot_images([image_stack[slice_idx], cad_mask_cleaned, res_mask[slice_idx]], contrast= True, show= False, path_output= path_out + "segmentation_slice_{0}".format(slice_idx +1))

    masks = {
        "cad_mask" : res_mask,
        "analysis_area_mask" : cad_mask_cleaned
    }

    return masks


def create_tiff_detection_check(cy3: np.ndarray, spots, path_output, dot_size = 2):
    """
    Creates 3D tiff image with cy3 channel and spot detected.
    """

    if cy3.ndim == 3 : 
        channel = stack.maximum_projection(cy3)

    if len(spots) == 0 : return 0

    z,y,x = list(zip(*spots))
    spots_mask = np.zeros_like(channel)
    spots_mask[y, x] = 1
    
    #enlarging dots
    if dot_size > 1 : spots_mask = binary_dilation(spots_mask, iterations= dot_size-1)

    spots_mask = stack.rescale(np.array(spots_mask, dtype = channel.dtype))
    
    im = np.zeros([2] + list(channel.shape))
    im[0,:,:] = channel
    im[1,:,:] = spots_mask

    stack.rescale(channel, channel_to_stretch= 0)
    stack.save_image(im, path_output, extension= 'tif')



def output_spot_tiffvisual(channel, spots, path_output, dot_size = 3, rescale = True):
    """Outputs a tiff image with one channel being {channel} and the other a mask containing dots where sports are located.
    
    Parameters
    ----------
        channel : np.ndarray
            3D monochannel image
        spots :  
        path_output : str
        dot_size : int
            in pixels
    """
    stack.check_parameter(channel = (np.ndarray),spots= (list, np.ndarray), path_output = (str), dot_size = (int))
    stack.check_array(channel, ndim= [2,3])
    if channel.ndim == 3 : 
        channel = stack.maximum_projection(channel)
    if len(spots[0]) == 3 : 
        new_spots = []
        for i in range(0,len(spots)) : new_spots += [[spots[i][1], spots[i][2]]] 
        spots = new_spots

    

    spots_mask = np.zeros_like(channel)
    for spot in new_spots :
        spots_mask[spot[0], spot[1]] = 1

    
    #enlarging dots
    if dot_size > 1 : spots_mask = binary_dilation(spots_mask, iterations= dot_size-1)


    spots_mask = stack.rescale(np.array(spots_mask, dtype = channel.dtype))
    
    im = np.zeros([2] + list(channel.shape))
    im[0,:,:] = channel
    im[1,:,:] = spots_mask

    if rescale : channel = stack.rescale(channel, channel_to_stretch= 0)
    stack.save_image(im, path_output, extension= 'tif')



def count_nucleus(nucleus_label: np.ndarray, analysis_area) :

    label = nucleus_label.copy()
    label[~analysis_area] = 0
    count = len(np.unique(label)) - 1 #removing 0 region in count

    return count