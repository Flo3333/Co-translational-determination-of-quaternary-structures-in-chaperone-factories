import os
import numpy as np
import pandas as pd
from czifile import imread
from .measures import cluster_localisation, spots_colocalisation

def create_Input(path_input) :
    """
    Create Input data frame for Soha colloc analysis.
    """
    dirlist = os.listdir(path_input)
    dirlist_copy = dirlist.copy()
    for file in dirlist_copy :
        if not file.endswith(".czi") : dirlist.remove(file)
    filename = dirlist
    ids = np.arange(len(dirlist))

    dataframe = pd.DataFrame({
        "id" : ids,
        "filename" : filename,
        "path" : [path_input] * len(dirlist)
    })

    return dataframe


def compute_fov(id= np.NaN, rna_name= np.NaN, image_shape= np.NaN, cell_number= np.NaN, cy3_spots= np.NaN, egfp_spots= np.NaN, clusters_dataframe= np.NaN) :

    if isinstance(cy3_spots, (list,tuple, np.ndarray)) :
        cy3_count = len(cy3_spots)
    else : cy3_count = np.NaN
    if isinstance(egfp_spots, (list,tuple, np.ndarray)) :
        egfp_count = len(egfp_spots)
    else : egfp_count = np.NaN

    if not isinstance(clusters_dataframe, pd.DataFrame) :
        cluster_spot_count = np.NaN
        cluster_number = np.NaN
    else :
        cluster_spot_count = cluster_localisation(clusters_dataframe)
        cluster_number = len(clusters_dataframe)
    
    if isinstance(cy3_spots, (list,tuple,np.ndarray)) and isinstance(egfp_spots, (list,tuple,np.ndarray)) and isinstance(image_shape, (list,tuple,np.ndarray)) :
        count = spots_colocalisation(image_shape, cy3_spots, egfp_spots)
    else : count = np.NaN

    fraction_colocalising = count / cy3_count
    assert all(fraction_colocalising <= 1)
    fraction_in_cluster = cluster_spot_count / cy3_count
    assert all(fraction_in_cluster <= 1)

    df = pd.DataFrame({
        "id" : [id]
        ,"rna_name" : [rna_name]
        ,"cell_number" : [cell_number]
        ,"egfp_spots" : [egfp_count]
        ,"cluster_number" : [cluster_number]
        ,"cy3_spots" : [cy3_count]
        ,"colocalisation_count" : [count]
        ,"fraction_colocalising" : [fraction_colocalising]
        ,"cluster_spot_count" : [cluster_spot_count]
        ,"fraction_in_cluster" : [fraction_in_cluster]
    })
    
    return df


def get_image_as_gen(Acquisition: pd.DataFrame, channel = None, Query_index = None, crop= None) :
    if type(Query_index) != type(None) : index = Query_index.copy()
    else : index = Acquisition.sort_index().index
    
    if type(crop) != type(None) :
        if not isinstance(crop, (list,tuple)) : raise TypeError("Wrong type for crop should be list or tuple it is {0}".format(type(crop)))
        if not len(crop) == 2 : raise ValueError("Expected 2 elements in crop argument : [start,stop]")

        min_z,max_z = crop

    for idx in index :
        im = imread(Acquisition.at[idx, "path"] + Acquisition.at[idx, "filename"])
        if channel != None : im = im[channel]
        assert im.shape[3] == 1
        im = im.reshape(im.shape[:3])
        if type(crop) != type(None) : 
            if min_z == None : min_z = 0
            if max_z == None : max_z = im.shape[0]

            im = im[min_z:max_z]

        yield im