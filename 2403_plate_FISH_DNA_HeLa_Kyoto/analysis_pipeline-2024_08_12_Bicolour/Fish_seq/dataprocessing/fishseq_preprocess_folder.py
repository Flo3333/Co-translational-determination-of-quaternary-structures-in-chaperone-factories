"""
Script to assert input data is conformed data is consistent with scripts.

Expected

- Run Folder (lvl1)
|
|
|
--- Cy3_Z-stacks (lvl2)
|
--|-- Location-01 # One folder per fov (lvl3) 
--|-- . . .
--|-- Location-0.
|
|
|
-----|-- img000.*.tif # even numbers for RNA
-----|-- img001.*.tif # odd numbers for washes

|
--- Dapi_Z-stacks (lvl2)
|
--|-- Location-01 # One folder per fov (lvl3)
--|-- . . .
--|-- Location-0.
----|
----|-- img.tif # 1 image for DAPI

"""

import os
import numpy as np
from skimage import io
from tqdm import tqdm

def _lvl1(RUN_PATH:str) :
    """
    returns True if ok else raise FileNotFoundError
    """

    dirlist = os.listdir(RUN_PATH)
    if not "Cy3_Z-stacks" in dirlist : raise FileNotFoundError("'Cy3_Z-stacks' folder not found in run folder.")
    if not "DAPI_Z-stacks" in dirlist : raise FileNotFoundError("'Dapi_Z-stacks' folder not found in run folder.")

    return True

def _lvl2(RUN_PATH:str) :
    """
    returns locations list of ok else raise ValueError
    """
    Cy3_dirlist = os.listdir(RUN_PATH + "/Cy3_Z-stacks")
    Dapi_dirlist = os.listdir(RUN_PATH + "/DAPI_Z-stacks")

    locations_Cy3 = []
    locations_Dapi = []

    for file in Cy3_dirlist :
        if 'Location-' in file : 
            locations_Cy3.append(file)

    for file in Dapi_dirlist :
        if 'Location-' in file : 
            locations_Dapi.append(file)

    locations_Cy3.sort()
    locations_Dapi.sort()

    if locations_Cy3 != locations_Dapi :
        raise ValueError("Missmatch between locations found in Cy3 folder and in Dapi folder")
    else :
        return locations_Cy3

def _lvl3(RUN_PATH, locations) :

    #Dapi
    for location in locations :
        dirlist = os.listdir(RUN_PATH + '/DAPI_Z-stacks/' + location)
        if len(dirlist) == 0 : raise FileNotFoundError("Dapi acquisition not found for location : {0}".format(location))
        elif len(dirlist) > 1 : raise FileNotFoundError("More than 1 dapi stack for location : {0}".format(location))
    
    #Cy3
    file_number = []
    file_dict = {}
    for location in locations :
        dirlist = os.listdir(RUN_PATH + '/Cy3_Z-stacks/' + location)
        file_dict[location] = dirlist.copy()
        file_dict[location].sort()
        for file in dirlist : 
            if not file.endswith(".ome.tif") or file.startswith("._") : file_dict[location].remove(file)
        file_number.append(len(file_dict[location]))
    assert len(np.unique(file_number)) == 1, "Different file numbers found for Cy3 Z-stacks amongst locations : {0}".format(np.unique(file_number))

    return file_dict

def assert_run_folder_integrity(run_path) :
    _lvl1(run_path)
    locations_list = _lvl2(run_path)
    file_dict = _lvl3(run_path, locations_list)

    return file_dict

def _assert_constant_shape(images:'list[np.ndarray]') :

    shapes = [image.shape for image in images]
    print(shapes)
    assert len(np.unique(shapes)) == 1, "images shape are not uniform in same location : {0}".format(np.unique(shapes))

    return shapes[0]

def create_multichannel_tif(location_folder_path:str, location:str, file_dict:dict, filename = 'img000_000_000000_0000000000.ome.tif') :
    
    dapi =  io.imread(location_folder_path.replace("Cy3_Z-stacks", "DAPI_Z-stacks" + "/{0}/{1}".format(location, filename)))
    cy3 = io.imread(location_folder_path + "/{0}/{1}".format(location, filename))

    tif = np.append(cy3, dapi.reshape((1,) + dapi.shape), axis=0)

    return tif
    



RUN_PATH = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2024-03-05 - R2TP_Pool-1_Run1/"


#main
file_dict = assert_run_folder_integrity(RUN_PATH)
print("Folder intergrity : OK")

# Extract run information
locations = list(file_dict.keys())
location_number = len(locations)
cycle_number = len(file_dict[locations[0]])

print("{0} locations found.\nEach location contains {1} cycles.".format(location_number, cycle_number))


print("Preparing files for pipeline input.")

location_folder_path = RUN_PATH + "/Cy3_Z-stacks/"
for location in tqdm(locations) :
    tif = create_multichannel_tif(location_folder_path, location, file_dict)
    print(tif.shape)
