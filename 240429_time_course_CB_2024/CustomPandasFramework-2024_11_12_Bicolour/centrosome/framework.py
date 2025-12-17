import os, re
import numpy as np
import pandas as pd
import pbwrap.quantification as quant
from CustomPandasFramework.operations import add_data
from pbwrap.utils import open_image as read_image
from ..utils import check_parameter
from ..integrity.Errors import InputError
from scipy.ndimage import distance_transform_edt

def make_Input(path_in:str, 
               channel,
               regex :str,
               groups : list,
               idx_keys : 'tuple[str]'
               ) :
    """
    Creates Input DataFrame from files placed in the input folder.
    """
    
    check_parameter(path_in = str)
    if not path_in.endswith('/') : path_in += '/'
    if type(groups) != list : groups = list(groups)
    if type(idx_keys) != list : idx_keys = list(idx_keys)
    groups = ['filename'] + groups

    filenames_list = os.listdir(path_in)
    data = []

    for file in filenames_list :
        match_re = re.findall(regex, file)

        if len(match_re) == 0 :
            print("{0} remove from pipeline (not matching regex)".format(match_re))
            continue
        elif len(match_re) > 0 : 
            caught_groups = list(match_re[0])
            if len(groups) -1 != len(caught_groups) : raise ValueError("Regex matched a file but the number of groups in input is inconsistent with regex capturing group number.\n groups : {0}, match_re : {1}".format(len(groups)-1, len(caught_groups)))
            data += [[file] + caught_groups]
    
    Input = pd.DataFrame(columns=groups, data=data)
    file_number = len(Input)
    assert len(Input) == len(Input.drop_duplicates()), "Duplicates found in Input."
    
    if 'channel' in Input.columns : 
        channel_less_Input = Input.drop(['channel'],axis=1)
    else :
        channel_less_Input = Input
    
    channel_less_Input = channel_less_Input.drop(['filename'],axis=1).drop_duplicates()
    merge_columns = channel_less_Input.columns
    channel_less_Input = channel_less_Input.reset_index(drop=False, names='acquisition_id')
    Input = pd.merge(
        channel_less_Input, #This adds id
        Input, #This adds channel & dupplicates rows
        on=list(merge_columns)
    )
    Input = Input.set_index(idx_keys,verify_integrity=False, drop=False).sort_index() #Expected not unique index

    assert len(Input) == file_number, "Acquisition number is not consistent after merge."

    return Input

def open_image(Input: pd.DataFrame,input_path:str, channel, idx, output_list= False, crop_z=None) :

    check_parameter(Input = pd.DataFrame, channel = (str, type(None)), input_path=str)
    if not input_path.endswith('/') : input_path += '/'

    view = Input.loc[idx]
    files = view.loc['filename']
    if isinstance(files, str) : files = [files]
    if type(channel) != type(None) : view = view.loc[view['channel'] == channel]
    if output_list :
        if type(crop_z) != tuple :
            return [read_image(input_path + filename) for filename in files]
        elif len(crop_z) == 2 :
            return [read_image(input_path + filename)[crop_z[0] : crop_z[1]] for filename in files]

    else : 
        if type(crop_z) != tuple :
            return (read_image(input_path + filename) for filename in files)
        elif len(crop_z) == 2 :
            return (read_image(input_path + filename)[crop_z[0] : crop_z[1]] for filename in files)


def add_cell(Cell_df, acquisition_id, cell, dapi_stack, voxel_size = (300,103,103)):
    """
    Returns DataFrame with expected Cell datashape containing all cell level features. Features are computed using bigFish built in functions.
    

    Parameters
    ----------
        acquisition_id : int
            Unique identifier for current acquisition.
        cell : dict
            Dictionary computed from bigFish.multistack.extract_cell
        dapi : raw signal from dapi channel : ndim = 3
    
    Returns
    -------
        new_Cell : pd.Dataframe
    """
    
    check_parameter(Cell_df = pd.DataFrame)
    with np.errstate(divide='ignore', invalid= 'ignore') :
        Cell_result = quant.hub_compute_Cell(
            acquisition_id= acquisition_id, 
            cell= cell, 
            dapi_stack= dapi_stack, 
            voxel_size=voxel_size, 
            pipeline_name= 'centrosome'
            )
        
        Cell = add_data(Cell_df, Cell_result)

    return Cell

def add_fov(Acquisition: pd.DataFrame, Input: pd.DataFrame, acquisition_id, rna_threshold, cell_computed, cell_discarded):
    
    if 'channel' in Input.columns : 
        acquisition = Input.drop('channel',axis=1)
        acquisition = acquisition.loc[acquisition_id].iloc[0] # We use iat because there can be more than 1 acquisition if there are several channels
    else :
        acquisition = Input.loc[acquisition_id]
    acquisition['id'] = acquisition_id
    acquisition['rna_threshold'] = rna_threshold
    acquisition['number_cell_computed'] = cell_computed
    acquisition['number_cell_discarded'] = cell_discarded
    acquisition = pd.DataFrame(
        columns= acquisition.index,
        data= [tuple(acquisition.values)],
        index= [acquisition_id],
    )

    acquisition =  pd.concat([Acquisition, acquisition], axis= 0)


    return acquisition

def empty_frame() :
    return pd.DataFrame()

def add_Spots(
            Spots_df : pd.DataFrame,
            Cell_df : pd.DataFrame,
            spots_fov: pd.DataFrame, 
            signal: np.ndarray, 
            nucleus_mask :np.ndarray, 
            cell_label: np.ndarray, 
            acquisition_id : int,
            centrosome_coords: np.ndarray = None, 
            other_mask: np.ndarray = None,
            voxel_size = None,
            ) :
    
    if nucleus_mask.dtype != bool : raise ValueError("mask are expected to be boolean arrays with 0 as background.")
    if type(other_mask) != type(None) : 
        if other_mask.dtype != bool : raise ValueError("mask are expected to be boolean arrays with 0 as background.")
    if (type(other_mask) != type(None) or type(centrosome_coords) != type(None)) and type(voxel_size) == type(None) : raise ValueError("voxel size is required if centrosome coords or other mask is passed.")
    


    Spots = pd.DataFrame(columns=['acquisition_id','spot_id','Z','Y','X','intensity', 'cell_label','in_nucleus','distance_centrosome','in_other_mask','distance_other_mask','cluster_id'])
    Spots['in_nucleus'] = Spots['in_nucleus'].astype(bool)
    Spots['in_other_mask'] = Spots['in_other_mask'].astype(bool)

    Spots['acquisition_id'] = [acquisition_id] * len(spots_fov)

    Z,Y,X = spots_fov['z'],spots_fov['y'],spots_fov['x'],
    intensity = signal[Z,Y,X]
    label = cell_label[Y,X]
    in_nucleus = nucleus_mask[Y,X]
    if type(other_mask) != type(None) : 
        if other_mask.ndim == 2 : in_other_mask = other_mask[Y,X]
        else : in_other_mask = other_mask[Z,Y,X]
    Spots['Z'] = Z
    Spots['Y'] = Y
    Spots['X'] = X
    Spots['intensity'] = intensity
    Spots['cell_label'] = label
    Spots['in_nucleus'] = in_nucleus
    Spots['cluster_id'] = spots_fov['cluster_id']

    #distance maps
    if type(centrosome_coords) != type(None) :
        if len(centrosome_coords) >0 :
            centrosome_array = np.ones_like(signal, dtype=bool)
            Z_centrosome, Y_centrosome, X_centrosome = zip(*centrosome_coords)
            centrosome_array[Z_centrosome,Y_centrosome,X_centrosome] = 0
            distance_map_centrosome = distance_transform_edt(centrosome_array, sampling= voxel_size)
            distance_centrosome = distance_map_centrosome[Z,Y,X]
            Spots['distance_centrosome'] = distance_centrosome
        else :
            Spots['distance_centrosome'] = np.NaN

    if type(other_mask) != type(None) :
        in_other_mask = other_mask[Z,Y,X]
        distance_map_other_mask = distance_transform_edt(~other_mask, sampling=voxel_size)
        distance_other_mask = distance_map_other_mask[Z,Y,X]
        Spots['in_other_mask'] = in_other_mask
        Spots['distance_other_mask'] = distance_other_mask

    #Getting cell_label
    Spots = pd.merge(
        Spots,
        Cell_df.loc[:,['id','acquisition_id','label']].rename(columns={'id' : 'cell_id', 'label' : 'cell_label'}),
        how= 'inner', #this will drop spots that are in discarded cells.
        on= ['acquisition_id','cell_label']
    )

    #Adding to previous results
    Spots = pd.concat([
        Spots_df,
        Spots
    ],axis=0).drop('spot_id', axis=1).reset_index(drop=True).reset_index(drop=False, names='spot_id')

    return Spots