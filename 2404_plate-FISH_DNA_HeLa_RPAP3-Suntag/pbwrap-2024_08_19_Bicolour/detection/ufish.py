"""
Submodule related to spot detection using UFish package.
"""

import numpy as np
from ufish.api import UFish
from ..utils import _create_counting_kernel
import bigfish.plot as plot

def init_ufish_model() -> UFish:
    """
    Will initiate a UFish model for inference with Pytorch.
    """

    ufish = UFish(cuda=True)
    ufish.init_model(
        model_type= 'ufish',
        channel_numbers = [16,16,16]
    )

    return ufish


def convert_radius_pxnumber(spot_radius: tuple, voxel_size:tuple) :
    kernel = _create_counting_kernel(spot_radius, voxel_size)
    area = kernel.sum()

    return area

def infer_spots(
        model : UFish,
        image : np.ndarray,
        spot_radius : tuple,
        voxel_size : tuple,
        axes : str = None,
        method : str = 'cc_center',
) -> np.ndarray :
    """
    Detect spot using ufish package (U-net neural network), it needs to have setup the GPU first.

    Parameters
    ----------
        model : UFISH class,
            model for spot inference can be obtained with by calling `init_ufish_model` function.
        image: np.ndarray
            Image to infer spot to. Should 2D or 3D if you have multichannel 2D image please call the function separetly for each channel.
        spot_radius : tuple(z,y,x) or tuple(y,x)
            expected radius of a spot in nanometers
        voxel_size : tuple(z,y,x) or tuple(y,x)
           size of a voxel
        axes : str
            axes for image dimension, example 'zyx' if none ufish will try to guess.

    Returns
    -------
        spots : np.ndarray
            array of dim 2 : shape = (spot number, dimensions number); dtype = int.
    """

    if image.ndim not in [2,3] : 
        raise ValueError("Only 2d or 3d images supported.")
    
    if image.ndim != len(voxel_size) :
        raise ValueError("Image dimension and voxel size length must match.")


    cc_size_threshold = convert_radius_pxnumber(spot_radius, voxel_size)

    spots,im = model.predict(
        image,
        axes=axes,
        spots_calling_method=method,
        cc_size_threshold=20
    )

    if image.ndim == 2 : 
        spots = np.array(list(zip(spots['axis-0'], spots['axis-1'])), dtype=int)
    if image.ndim == 3 : 
        spots = np.array(list(zip(spots['axis-0'], spots['axis-1'], spots['axis-2'])), dtype=int)

    return spots