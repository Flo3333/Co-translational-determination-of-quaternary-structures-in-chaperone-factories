"""
This script aims at reading the input folder and preparing data folders and locations for next scripts.
"""

from FishSeq_parameters import RUN_PATH, FOLDER_KEYS, MAP_FILENAME, cycle_regex

import CustomPandasFramework.Fish_seq.folder_preprocessing as prepro
import pandas as pd
import os
import re

#Reading input folder.
file_dict = prepro.assert_run_folder_integrity(
    run_path=RUN_PATH,
    fish_folder=FOLDER_KEYS.get('fish'),
    nucleus_folder=FOLDER_KEYS.get('nucleus')
    )
location_list = list(file_dict.keys())
location_list.sort()
print("{0} locations found.".format(len(location_list)))

#Init pandas DF
COLUMNS = [
    "acquisition_id",
    "location",
    "cycle",
    "full_path",
    "dapi_full_path",
    ]
Acquisition = pd.DataFrame(columns=COLUMNS)

file_index = 0
for location_index, location in enumerate(location_list) :
    
    #Get dapi_path
    dapi_full_path = RUN_PATH + "/{0}/{1}/".format(FOLDER_KEYS.get('nucleus'), location)
    assert len(os.listdir(dapi_full_path)) == 1
    dapi_full_path += os.listdir(dapi_full_path)[0]
    assert os.path.isfile(dapi_full_path)
    
    #Get fish_path
    fish_path = RUN_PATH + "/{0}/{1}/".format(FOLDER_KEYS.get('fish'), location)
    fish_path_list = os.listdir(fish_path)
    fish_path_list.sort()
    for file in fish_path_list :
        cycle = int(re.findall(cycle_regex, file)[0])
        fish_full_path = fish_path + file
        assert os.path.isfile(dapi_full_path)
        Acquisition.loc[file_index] = [file_index, location, cycle, fish_full_path, dapi_full_path]
        file_index += 1

#Cycle mapping
map = pd.read_excel(RUN_PATH + '/' + MAP_FILENAME)
assert all(Acquisition['cycle'].isin(map.iloc[:,0])), "Some cycle are not found in map"
map_cycle_key = map.columns[0]
Acquisition = pd.merge(
    left=Acquisition,
    right=map,
    left_on='cycle',
    right_on=map_cycle_key
).sort_values('acquisition_id').reset_index(drop=True)



#Output
save_path = RUN_PATH + '/result_tables/'
os.makedirs(save_path, exist_ok=True)
Acquisition.to_excel(save_path + '/Acquisition.xlsx')
Acquisition.to_feather(save_path + '/Acquisition.feather')
print("Done")