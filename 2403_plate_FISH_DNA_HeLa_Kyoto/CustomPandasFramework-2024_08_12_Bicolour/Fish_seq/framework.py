import os
import re
import pandas as pd
import numpy as np
import pbwrap.quantification.measures as measure
import CustomPandasFramework.utils as utils
import bigfish.stack as stack
from czifile import imread


def create_Input(input_path, regex, columns, multichannel) :
    """
    regex with as many capturing groups as columns names should be passed
    """

    if not multichannel :
        if 'channel' not in columns : raise KeyError("For set of monochannel images 'channel' is expected in columns.")

    #Removing from input files not tif or not czi
    filelist = os.listdir(input_path)
    non_filter_list:list[str] = filelist.copy()
    for file in non_filter_list :
        if not file.endswith('czi') and not file.endswith('tiff') : filelist.remove(file)
    del non_filter_list
    
    data = []
    for filename in filelist :
        try: 
            groups_caught = list(re.findall(regex, filename)[0])
            if len(groups_caught) != len(columns) : raise KeyError("regex capturings groups number and columns number doesn't match\nCaught : len {0} / {1}\nColumns : len {2} : {3}".format(len(groups_caught), groups_caught,  len(columns), columns))
            groups_caught = [filename] + groups_caught

            data.append(groups_caught)
        except IndexError :
            continue

    data_columns = ['filename'] + columns
    Input = pd.DataFrame(columns= data_columns, data= data)

    if not multichannel : #For set of monochannel images we group together inputs from same field of view
        for channel in Input['channel'].unique() :
            Input['filename'] = Input['filename'].str.replace(channel, '')
        Input = Input.groupby('filename').first().reset_index().drop(["channel"], axis= 1)

    Input = Input.reset_index(drop=True).reset_index(drop= False).rename(columns={'index' : 'id'})

    if 'rna1' not in Input.columns :
        raise ValueError("Please rename main rna target column to 'rna1'.")

    if 'rna2' not in Input.columns : 
        print("Warning : no rna2 metadata found : colocalisation will not yield any results. Only rna1 cluster analysis will be performed")
        Input['rna2'] = np.NaN

    return Input

def create_Acquisiton(Input: pd.DataFrame, ) :
    """
    """
    Acquisition = Input.copy().rename(columns={"id" : "acquisition_id"})
    Acquisition['Runtime_date'] = utils.get_datetime()
    Acquisition['colocalisation_distance'] = np.NaN

    #Thresholds
    Acquisition['rna1_threshold'] = np.NaN
    Acquisition['rna2_threshold'] = np.NaN
    Acquisition['suntag_threshold'] = np.NaN

    #boolean values
    Acquisition['rna1_deconvolution_sucess'] = False
    Acquisition['rna2_deconvolution_sucess'] = False
    Acquisition['suntag_deconvolution_sucess'] = False

    return Acquisition


def get_image_as_gen(Input: pd.DataFrame, path_in:str, index_list, multichannel, channels = None, _crop_pixel = 0) :
    
    """
    Parameters
    ----------
        Input : pd.DataFrame from create Input
        path_in : str
            path to image files
        index_list : list[int]
            list of index from Input making a group (from make_analysis_group)
        channel : int or list of int or None
            channel to yield in the generator. Default None : takes all channels
    """

    if isinstance(channels, int) and multichannel:
        channels = [channels]
    elif isinstance(channels, list) and multichannel:
        for elmt in channels : 
            if not isinstance(elmt, int) : raise TypeError('Non int element was found in channels argument for multichannel image')
    elif isinstance(channels, type(None)) and multichannel:
        pass
    elif isinstance(channels, str) and not multichannel :
        pass

    elif isinstance(channels, list) and not multichannel :
        for elmt in channels : 
            if not isinstance(elmt, str) : raise TypeError('Non str element was found in channels argument for set of monochannel images')

    elif multichannel :
        raise TypeError('For multichannel image, channels argument should be of type list[int] or int. It is {0}'.format(type(channels)))
    
    else :
        raise TypeError('For set of monochannel images, argument should be str. It is {0}'.format(type(channels)))


    df = Input.loc[index_list,:]
    if not path_in.endswith('/') : path_in += '/' 
    for index in df.index :
        filename:str = df.at[index, 'filename']
        if not multichannel :
            filename = filename.replace('--', '-{0}-'.format(channels))
            if filename.endswith('.czi') :
                im = imread(path_in + filename)
            else :
                im = stack.read_image(path_in + filename)
            yield im
        
        else :
            if filename.endswith('.czi') :
                im = imread(path_in + filename)
            else :
                im = stack.read_image(path_in + filename)
            if _crop_pixel != 0 : im = im[:,:,_crop_pixel:im.shape[2] - _crop_pixel, _crop_pixel:]
            if type(channels) == type(None) : channels = [channel for channel in range(0,len(im))]
            if len(im.shape) > 4 : im= im.reshape(im.shape[:4])

            if len(channels) == 1 :
                yield im[channels[0]]
            else : yield im

def make_analysis_group(Input:pd.DataFrame, group_keys= ['rna1', 'rna2']) :
    return Input.groupby(group_keys)['id'].apply(list)


def tag_id(Dataframe: pd.DataFrame, tag, id_list:list = ['id']) :
    """
    if an id is NA then the new id is NA as well.
    """

    df = Dataframe.copy()

    for id_name in id_list :
        if id_name not in df.columns : raise KeyError('{0} column from id_list wasnt found in dataframe.'.format(id_name))
        na_idx = df.query("{0}.isna()".format(id_name)).index
        df[id_name] = list(zip(tag, df[id_name]))
        df.loc[na_idx,id_name] = np.NaN
    
    return df


def read_feather(path, tag:str, **read_kargs) :
    df = pd.read_feather(path, **read_kargs)
    df = tag_id(df, tag)
    return df

def load_DataFrames(input:str, DataFrame_name:str) :

    if not input.endswith('/') : input += '/'
    dirlist = os.listdir(input)

    fullpath_list = []
    
    for dir in dirlist : 
        if not os.path.isdir(input + dir) : raise NotADirectoryError("{0} was found in input and is not a directory".format(dir))
        table_list = os.listdir(input+dir)
        if DataFrame_name in table_list : fullpath_list.append((table_list.index(input+dir+DataFrame_name), dir))

    df = pd.concat([read_feather(path= table_fullpath, tag= table_dir) for table_fullpath, table_dir in fullpath_list], axis= 1).set_index('id')

    return df