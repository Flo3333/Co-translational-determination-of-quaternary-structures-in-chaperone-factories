import os, re
import numpy as np
import pandas as pd
import pbwrap.quantification as quant
from CustomPandasFramework.operations import add_data
from bigfish.stack import read_image
from ..utils import check_parameter
from ..integrity.Errors import InputError

def make_Input(path_in:str, channel,
               regex_treated= "(\w*)-([-\w]*)_([\w]*)_([\w\s-]*)f(\d{2})-([\w\s]*)-.*\.tiff",
               regex_untreated = "(\w*)-([-\w]*)_([\w\s-]*)f(\d{2})-([\w\s]*)-.*\.tiff"
               ) :
    """
    Creates Input DataFrame from files placed in the input folder.
    index : ['rna_name', 'treatment', 'fov']
    Keys : ['id', 'filename', 'cell_line', 'channel']
    """
    
    check_parameter(path_in = str)
    if not path_in.endswith('/') : path_in += '/'

    filenames_list = os.listdir(path_in)
    well_list = []
    cell_line_list = []
    treatment_list = []
    rna_list = []
    fov_number_list = []
    channel_list = []
    filenames = filenames_list.copy()

    for file in filenames :
        match_re = re.findall(regex_treated, file)

        if len(match_re) == 0 :
            match_re = re.findall(regex_untreated, file)
            groups = match_re[0]
            well_list.append(groups[0])
            cell_line_list.append(groups[1])
            treatment_list.append("untreated")
            rna_list.append(groups[2])
            fov_number_list.append(int(groups[3]))
            channel_list.append(groups[4])
        elif len(match_re) > 0 : 
            groups = match_re[0]
            well_list.append(groups[0])
            cell_line_list.append(groups[1])
            treatment_list.append(groups[2])
            rna_list.append(groups[3])
            fov_number_list.append(int(groups[4]))
            channel_list.append(groups[5])

        else : filenames_list.remove(file)
    
    Input = pd.DataFrame({
        "well" : well_list,
        "filename" : filenames_list,
        "cell_line" : cell_line_list,
        "treatment" : treatment_list,
        "rna_name" : rna_list,
        "fov" : fov_number_list,
        "channel" : channel_list
    })

    #ids
    group = Input.groupby(['rna_name', 'treatment', 'well', 'fov'])['channel'].count()
    index = group.index
    id_number = len(group)
    id_df = pd.Series(np.arange(id_number), index= index).rename('id')
    Input = pd.merge(id_df, Input, on= ['rna_name', 'treatment', 'well', 'fov'])


    channel_filter = Input.query('channel in {0}'.format(channel)).index
    Input = Input.loc[channel_filter,:].sort_values(['rna_name', 'treatment', 'fov','channel'])
    
    #Integrity checks

    # Coherence filenumbers and channel numbers
    truth = [len(Input[Input['channel'] == chan] == len(Input) / len(channel)) for chan in channel]
    if not all(truth) : raise InputError
    
    # Index is unique
    assert len(Input.query('channel == "{0}"'.format(channel[0]))) == len(Input['id'].unique()), "id is not unique when making Input df."

    # Index refers to unique experiment
    Input = Input.set_index(['rna_name', 'treatment', 'well', 'fov', 'channel'], verify_integrity=True)


    return Input

def open_image(Input: pd.DataFrame,input_path:str, channel, rna_name= None, output_list= False, crop_z=None) :

    check_parameter(Input = pd.DataFrame, channel = str, input_path=str)
    if type(rna_name) == type(None) : rna_name = list(Input['rna_name'].unique())
    if not input_path.endswith('/') : input_path += '/'

    view = Input.query('rna_name == "{0}" and channel == "{1}"'.format(rna_name, channel)).sort_index()
    if output_list :
        if type(crop_z) != tuple :
            return [read_image(input_path + filename) for filename in view['filename']]
        elif len(crop_z) == 2 :
            return [read_image(input_path + filename)[crop_z[0] : crop_z[1]] for filename in view['filename']]

    else : 
        if type(crop_z) != tuple :
            return (read_image(input_path + filename) for filename in view['filename'])
        elif len(crop_z) == 2 :
            return (read_image(input_path + filename)[crop_z[0] : crop_z[1]] for filename in view['filename'])


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
        Cell_result = quant.hub_compute_Cell(acquisition_id= acquisition_id, cell= cell, dapi_stack= dapi_stack, voxel_size=voxel_size, pipeline_name= 'centrosome')
        
        Cell = add_data(Cell_df, Cell_result)

    return Cell

def add_fov(Acquisition: pd.DataFrame, acquisition_id, cell_line, rna, well, fov, treatement, rna_threshold, cell_computed, cell_discarded):
    df = pd.DataFrame({
        'id' : [acquisition_id],
        'cell_line' : [cell_line],
        'rna' : [rna],
        'well' : [well],
        'fov' : [fov],
        'treatment' : [treatement],
        'rna_thresholrd' : [rna_threshold],
        'number_cell_computed' : [cell_computed],
        'number_cell_discarded' : [cell_discarded]
                        })
    return pd.concat([Acquisition, df], axis= 0)

def empty_frame() :
    return pd.DataFrame()