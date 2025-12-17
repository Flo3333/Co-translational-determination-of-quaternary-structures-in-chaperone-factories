import pandas as pd
import numpy as np
import os
from czifile import imread



def create_Input(data_path:str) :
    """
    Creates an Acquisition DataFrame using a recursive file search in a folder which architecture is expected to be as follows.

    Expected file architecture :
    ----------------------------
        --> Given path /
            --> Protein_name ex : (b-catenin smFISH)
                --> Untreated / Puromycin
                    --> files.czi

    Returns
    -------
        Input : pd.DataFrame
            Pandas Dataframes with keys : (id, filename, path, rna_name, WithPuromycin)
            
            )        
    """
    if not isinstance(data_path, str) : raise TypeError("Incorrect type for data_path arg, expected str. It is {0}".format(type(data_path)))
    if not data_path.endswith('/') : data_path += '/'

    main_dir_list = os.listdir(data_path)

    ids_list = []
    filenames = []
    protein_names = []
    path_list= []
    puromycin_boollist = []
    max_id = 0

    #Recursive reading
    for protein in main_dir_list :
        if not os.path.isdir(data_path + protein) : continue
        for folder in os.listdir(data_path + protein) :
            if folder not in ["Puromycin", "Untreated"] : continue
            filelist = os.listdir(data_path + protein + '/' + folder)
            file_number = len(filelist)

            ids = np.arange(max_id,max_id + file_number)
            max_id = ids.max()
            ids_list.extend(list(ids))
            protein_names.extend([protein] * file_number)
            puromycin_boollist.extend([folder == "Puromycin"] * file_number)
            path_list.extend([data_path + protein + '/' + folder + '/']* file_number)
            filenames.extend(filelist)


    Input = pd.DataFrame(columns= ["id", "filename", "path", "rna_name", "WithPuromycin"], data= zip(*[ids_list, filenames, path_list, protein_names, puromycin_boollist ]))
    return Input


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