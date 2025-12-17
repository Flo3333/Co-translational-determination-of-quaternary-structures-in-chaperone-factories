"""
Submodule to prepare colocalisation for FishSeq data architecture.
"""

import pandas as pd
import numpy as np
import bigfish.multistack as multistack
import pbwrap.quantification.measures as measure
from itertools import product
from pbwrap.utils import _keep_spots_in_mask
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from scipy.ndimage import distance_transform_edt

def get_coloc_objects_names(
        detection_id_list : list,
        population = ('all', 'clustered', 'free')
) :
    """
    Return object names for coloc, ex : '16_all', '16_free', '16_clustered', '17_all'.... 
    """
    
    res = ['{0}|{1}'.format(detection_id, pop) for detection_id, pop in product(detection_id_list, population)]
    return res

def get_coloc_object_dict(
        Spots:pd.DataFrame,
        detection_id_list : list,
        population = ('all', 'clustered', 'free')

) :
        """
        Returns {detection_id}_{population} : corresponding spots dict

        PARAMETERS
        ##########
                Spots : pd.DataFrame
                detection_id_list : iterable
                population : population str value to look in Spots['population']; 'all' is a special value getting all spots regardeless of population.
        """
        
        if not 'population' in Spots.columns : raise KeyError('"population" column not found in spot please add a "population" column corresponding to values passed in this call population argument.')
        
        res = {}

        for detection_id, pop in product(detection_id_list, population) :
              
                if pop == 'all' :
                    spots = np.array(list(zip(
                          Spots.loc[Spots['detection_id'] == detection_id]['z'],
                          Spots.loc[Spots['detection_id'] == detection_id]['y'],
                          Spots.loc[Spots['detection_id'] == detection_id]['x'],

                    )), dtype=int)

                else :
                    spots = np.array(list(zip(
                          Spots.loc[(Spots['detection_id'] == detection_id) & (Spots['population'] == pop)]['z'],
                          Spots.loc[(Spots['detection_id'] == detection_id) & (Spots['population'] == pop)]['y'],
                          Spots.loc[(Spots['detection_id'] == detection_id) & (Spots['population'] == pop)]['x'],
                          
                    )), dtype=int)
                
                if len(spots) != 0 :
                        res['{0}|{1}'.format(detection_id,pop)] = spots
                else :
                        res['{0}|{1}'.format(detection_id,pop)] = np.empty(shape=(0,3),dtype=int)
                       
        return res



def _cell_colocalisation(
        cell,
        object1:str,
        object2:str,
        voxel_size,
        coloc_distance,
        z_shape,
        ) :
        
        coloc_dict = {}
        min_y, min_x, max_y, max_x = cell['bbox']
        shape = (z_shape, max_y-min_y, max_x-min_x)
        spots1 = cell[object1]
        spots2 = cell[object2]

        if len(spots1) == 0 or len(spots2) == 0 :
              coloc_dict['count'] = np.NaN

        coloc_dict['count'] = measure.spots_colocalisation(
              image_shape= shape,
              spot_list1=spots1,
              spot_list2=spots2,
              distance=coloc_distance,
              voxel_size=voxel_size
        )

        object1, pop1 = object1.split('|',maxsplit=1)
        object2, pop2 = object2.split('|',maxsplit=1)

        coloc_dict['object1'] = object1
        coloc_dict['object2'] = object2
        coloc_dict['population1'] = pop1
        coloc_dict['population2'] = pop2

        return coloc_dict

def _multi_thread_cell_coloc(
              cell,
              voxel_size,
              coloc_distance,
              z_shape,
              population,
              detection_number: int,
) :
        keys_list = list()
        Colocalisation = pd.DataFrame()
        min_y,min_x,max_y,max_x = cell['bbox']
        map_shape = (detection_number,z_shape, max_y-min_y, max_x-min_x)
        distance_map = np.ones(shape=map_shape, dtype=bool)
        map_index = 0
        index=0
        spots2_number = []

        for object_key in cell.keys() :
               
               if not population in object_key : 
                        continue
               elif len(cell[object_key]) == 0 : 
                        map_index +=1
                        keys_list.append(object_key)
                        spots2_number.append(0)
                        continue
               else :
                        keys_list.append(object_key)
                        spots2_number.append(len(cell[object_key]))
                        spots2 = np.array(cell[object_key])
                        spots2 = _keep_spots_in_mask(spots2, cell['cell_mask'])
                        z,y,x = zip(*spots2)
                        distance_map[map_index][z,y,x] = 0
                        distance_map[map_index] = distance_transform_edt(distance_map[map_index], sampling= voxel_size)
                        map_index +=1

        distance_map = distance_transform_edt(distance_map, sampling=(coloc_distance*10,) + voxel_size)
        distance_map = distance_map <= coloc_distance

        for object1 in cell.keys() :
                if '|' not in object1 : 
                       continue
                
                spots1 = np.array(cell[object1])
                if len(cell[object1]) > 0 :
                        spots1 = _keep_spots_in_mask(spots1, mask=cell['cell_mask'])
                object1, pop1 = object1.split('|',maxsplit=1)
                for map_index, object2 in enumerate(keys_list) :
                        coloc_dict = {}
                        truth_map = distance_map[map_index]

                        if len(spots1) == 0 :
                                coloc_dict['count'] = np.NaN
                        else :
                                z,y,x = zip(*spots1)
                                coloc_dict['count'] = (truth_map[z,y,x]).sum()

                        object2, pop2 = object2.split('|',maxsplit=1)

                        coloc_dict['detection_id1'] = object1
                        coloc_dict['detection_id2'] = object2
                        coloc_dict['population1'] = pop1
                        coloc_dict['population2'] = pop2
                        coloc_dict['label'] = cell['cell_id']
                        coloc_dict['spot1_number'] = len(spots1)
                        coloc_dict['spot2_number'] = spots2_number[map_index]
                
                        Colocalisation = pd.concat([
                              Colocalisation,
                              pd.DataFrame(coloc_dict, index = [index])
                        ],axis=0)
                        index +=1
               

        return Colocalisation

def compute_colocalisation(
            Spots,
            detection_id_list,
            population,
            voxel_size,
            coloc_distance,
            cell_label,
            location,
            z_shape,
            max_workers = 4,
            ) :
      
        detection_number = len(detection_id_list)
        object_dict = get_coloc_object_dict(
                Spots,
                detection_id_list=detection_id_list,
                population=population,
        )
        
        cell_extraction = multistack.extract_cell(
                cell_label=cell_label,
                ndim=3,
                rna_coord=None,
                others_coord= object_dict
            )
      
        cell_number = len(cell_extraction)
        Colocalisation = pd.DataFrame()
        for pop in population :
                with ThreadPoolExecutor(max_workers= max_workers) as executor :
                    cell_quantification_result = list(tqdm(executor.map(
                        _multi_thread_cell_coloc,
                        cell_extraction,
                        [voxel_size]*cell_number,
                        [coloc_distance]*cell_number,
                        [z_shape]*cell_number,
                        [pop] * cell_number,
                        [detection_number]*cell_number
                    ),total=cell_number, desc="colocalisation quantification ({0})".format(pop)))
                Colocalisation = pd.concat([Colocalisation] + cell_quantification_result, axis=0).reset_index(drop=True)
                Colocalisation['location'] = location

        return Colocalisation