"""
Submodule used to create measures for the spots colocalisations.
"""
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt


def reconstruct_boolean_signal(image_shape, spot_list: list)-> np.ndarray[bool] :
    signal = np.zeros(image_shape, dtype= bool)
    Z, Y, X = list(zip(*spot_list))
    signal[Z,Y,X] = True

    return signal




def spots_colocalisation(image_shape, spot_list1:list, spot_list2:list, distance: float, voxel_size) :
    """
    Return number of spots from spot_list1 located closer(large) than distance to at least one spot of spot_list2.

    Parameters
    ----------
        image_shape : tuple
        spot_list1 : list
        spot_list2 : list
        distance : nanometer
            distance in nanometer.
        voxel_size : (z,y,x) tuple
    """

    if len(image_shape) != 3 : raise ValueError("Image shape length should be 3 : (Z,Y,X)")

    signal1 = reconstruct_boolean_signal(image_shape, spot_list1)
    signal2 = reconstruct_boolean_signal(image_shape, spot_list2)
    mask = np.logical_not(signal2)
    distance_map = distance_transform_edt(mask, sampling= voxel_size)
    Z,Y,X = zip(*signal1)

    count = (distance_map[Z,Y,X] <= distance).sum()
    return count

def cluster_localisation(clusters_dataframe: pd.DataFrame) :
    return clusters_dataframe["spot_number"].sum()

# Testing

# image_shape = [10,10,10]
# X1 = list(np.random.randint(0,10,20))
# Y1 = list(np.random.randint(0,10,20))
# Z1 = list(np.random.randint(0,10,20))

# X2 = list(np.random.randint(0,10,20))
# Y2 = list(np.random.randint(0,10,20))
# Z2 = list(np.random.randint(0,10,20))

# X1 += [2]
# X2 += [2]
# Y1 += [2]
# Y2 += [2]
# Z1 += [2]
# Z2 += [2]

# spots_list1 = list(zip(Z1,Y1,X1))
# spots_list2 = list(zip(Z2,Y2,X2))

# print("image shape : ", image_shape)
# print("spots_list1 : ", spots_list1)
# print("spots_list2 : ", spots_list2)

# count = spots_colocalisation(image_shape, spots_list1, spots_list2)
# print("count : ", count)