import pandas as pd
import numpy as np
import itertools
import pbwrap.quantification.measures as measure
import pbwrap.quantification.hub as cell_quant
import bigfish.classification as features
import bigfish.multistack as multistack
from pbwrap.utils import get_centroids_list

def spots_quantif(rna1_free_spots: np.ndarray,
                  rna2_free_spots: np.ndarray,
                  suntag_free_spots: np.ndarray,
                  rna1_clustered_spots: np.ndarray,
                  rna2_clustered_spots: np.ndarray,
                  suntag_clustered_spots: np.ndarray,
                  cell_label : np.ndarray,
                  fov_cell : pd.DataFrame,
                  rna1_signal : np.ndarray,
                  rna2_signal : np.ndarray,
                  suntag_signal : np.ndarray,
                  ) :
    
    #Keep only spots in quantified cells.
    if len(rna1_free_spots) > 0 :
        Z,Y,X = zip(*rna1_free_spots)
        rna1_free_spots = rna1_free_spots[
            np.isin(cell_label[Y,X], fov_cell["cell_id"])
        ]
    if len(rna2_free_spots) > 0 :
        Z,Y,X = zip(*rna2_free_spots)
        rna2_free_spots = rna2_free_spots[
            np.isin(cell_label[Y,X], fov_cell["cell_id"])
        ]
    if len(suntag_free_spots) > 0 :
        Z,Y,X = zip(*suntag_free_spots)
        suntag_free_spots = suntag_free_spots[
            np.isin(cell_label[Y,X], fov_cell["cell_id"])
        ]

    if len(rna1_clustered_spots) > 0 :
        Z,Y,X = zip(*rna1_clustered_spots)
        rna1_clustered_spots = rna1_clustered_spots[
            np.isin(cell_label[Y,X], fov_cell["cell_id"])
        ]
    if len(rna2_clustered_spots) > 0 :
        Z,Y,X = zip(*rna2_clustered_spots)
        rna2_clustered_spots = rna2_clustered_spots[
            np.isin(cell_label[Y,X], fov_cell["cell_id"])
        ]
    if len(suntag_clustered_spots) > 0 :
        Z,Y,X = zip(*suntag_clustered_spots)
        suntag_clustered_spots = suntag_clustered_spots[
            np.isin(cell_label[Y,X], fov_cell["cell_id"])
        ]

    total_spot_number = len(rna1_free_spots) + len(rna2_free_spots) + len(suntag_free_spots) + len(rna1_clustered_spots) + len(rna2_clustered_spots) + len(suntag_clustered_spots)
    
    if len(rna1_free_spots) > 0 : 
        rna1_free_Z, rna1_free_Y, rna1_free_X = zip(*rna1_free_spots)
    else :
        rna1_free_Z, rna1_free_Y, rna1_free_X = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)

    if len(rna2_free_spots) > 0 : 
        rna2_free_Z, rna2_free_Y, rna2_free_X = zip(*rna2_free_spots)
    else :
        rna2_free_Z, rna2_free_Y, rna2_free_X = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)

    if len(suntag_free_spots) > 0 : 
        suntag_free_Z, suntag_free_Y, suntag_free_X = zip(*suntag_free_spots)
    else :
        suntag_free_Z, suntag_free_Y, suntag_free_X = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)

    if len(rna1_clustered_spots) > 0 : 
        rna1_clustered_Z, rna1_clustered_Y, rna1_clustered_X = zip(*rna1_clustered_spots)
    else :
        rna1_clustered_Z, rna1_clustered_Y, rna1_clustered_X = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)

    if len(rna2_clustered_spots) > 0 : 
        rna2_clustered_Z, rna2_clustered_Y, rna2_clustered_X = zip(*rna2_clustered_spots)
    else :
        rna2_clustered_Z, rna2_clustered_Y, rna2_clustered_X = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)

    if len(suntag_clustered_spots) > 0 : 
        suntag_clustered_Z, suntag_clustered_Y, suntag_clustered_X = zip(*suntag_clustered_spots)
    else :
        suntag_clustered_Z, suntag_clustered_Y, suntag_clustered_X = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)


    data = {
        "spots_id" : np.arange(total_spot_number),
        "cell_id" : list(cell_label[rna1_free_Y, rna1_free_X]) + list(cell_label[rna2_free_Y, rna2_free_X]) + list(cell_label[suntag_free_Y, suntag_free_X]) + list(cell_label[rna1_clustered_Y, rna1_clustered_X]) + list(cell_label[rna2_clustered_Y, rna2_clustered_X]) + list(cell_label[suntag_clustered_Y, suntag_clustered_X]),
        "coordinates" : list(rna1_free_spots) + list(rna2_free_spots) + list(suntag_free_spots) + list(rna1_clustered_spots) + list(rna2_clustered_spots) + list(suntag_clustered_spots),
        "spot_type" : ['rna1'] * len(rna1_free_spots) + ['rna2'] * len(rna2_free_spots) + ['suntag'] * len(suntag_free_spots) + ['rna1'] * len(rna1_clustered_spots) + ['rna2'] * len(rna2_clustered_spots) + ['suntag'] * len(suntag_clustered_spots),
        "population" : ['free'] * len(rna1_free_spots) + ['free'] * len(rna2_free_spots) + ['free'] * len(suntag_free_spots) + ['clustered'] * len(rna1_clustered_spots) + ['clustered'] * len(rna2_clustered_spots) + ['clustered'] * len(suntag_clustered_spots),
        "intensity" : np.concatenate([
            rna1_signal[rna1_free_Z, rna1_free_Y, rna1_free_X], 
            rna2_signal[rna2_free_Z, rna2_free_Y, rna2_free_X], 
            suntag_signal[suntag_free_Z, suntag_free_Y, suntag_free_X],
            rna1_signal[rna1_clustered_Z, rna1_clustered_Y, rna1_clustered_X], 
            rna2_signal[rna2_clustered_Z, rna2_clustered_Y, rna2_clustered_X], 
            suntag_signal[suntag_clustered_Z, suntag_clustered_Y, suntag_clustered_X]

            ])
    }

    spots_df = pd.DataFrame(data)
    spots_df["coordinates"] = spots_df['coordinates'].apply(tuple)

    return spots_df     

def cell_quantif(
                cell_label,
                nucleus_label,
                voxel_size,
                nucleus_signal,
                colocalisation_distance,
                
                rna1_spots, 
                rna2_spots, 
                suntag_spots, 
                rna1_cluster, 
                rna2_cluster, 
                suntag_cluster,
                rna1_clustered_spots,
                rna2_clustered_spots,
                suntag_clustered_spots,
                rna1_free_spots,
                rna2_free_spots,
                suntag_free_spots,

                rna1_signal,
                rna2_signal,
                suntag_signal,
                
                ):

    colloc_obj = {
            "rna1" : rna1_spots,
            "rna2" : rna2_spots,
            "suntag" : suntag_spots,
            "rna1_clustered" : rna1_clustered_spots,
            "rna2_clustered" : rna2_clustered_spots,
            "suntag_clustered" : suntag_clustered_spots,
            "rna1_free" : rna1_free_spots,
            "rna2_free" : rna2_free_spots,
            "suntag_free" : suntag_free_spots,
        }
    assert nucleus_signal.ndim == 3, "Nucleus signal given should be the 3D stack monochannel"
    z_shape = nucleus_signal.shape[0]

    coloc_df = pd.DataFrame(columns= _get_coloc_features_names(colloc_obj)) #empty
    del colloc_obj['rna1'] #To avoid repetition with rna1 that is passed as regulard rna_coord in extract_cell 
    rna1_df = pd.DataFrame(columns= cell_quant._get_minimal_columns_names(nucleus_quant=True)) #empty
    rna2_df = pd.DataFrame(columns= cell_quant._get_minimal_columns_names()) #empty
    suntag_df = pd.DataFrame(columns= cell_quant._get_minimal_columns_names()) #empty

    #During extract_cell a cell id is the value of the cell label, so we can use this id for merging data at the end.
    #Quantif rna1 + coloc
    for key, value in colloc_obj.copy().items() :
        if len(value) == 0 : del colloc_obj[key]

    colloc_obj['foci'] = rna1_cluster

    rna1_result = multistack.extract_cell(
        cell_label= cell_label,
        ndim= 3,
        others_image= {"smfish" : rna1_signal},
        nuc_label = nucleus_label,
        rna_coord= rna1_spots if len(rna1_spots) > 0 else np.array([0,0,0]), #else in case no spot detected so that coloc loop still works
        others_coord = colloc_obj
    )

    if len(rna1_spots) > 0 :
        for cell in rna1_result : 
            rna1_quantif = cell_quant._main_cell_quant(
                cell=cell,
                voxel_size=voxel_size,
                dapi_stack=nucleus_signal, #Need to be given only for first cell_quant

            )
            rna1_quantif['cluster_number'] = len(cell['foci'])

            rna1_df = pd.concat([
                rna1_df,
                rna1_quantif
            ],
            axis=0
            )
    
    for cell in rna1_result : 
        coloc_quantif = _cell_coloc_quantification(
                cell,
                colocalisation_distance=colocalisation_distance,
                voxel_size=voxel_size,
                z_shape= z_shape
            )

        coloc_df = pd.concat([
            coloc_df,
            coloc_quantif
        ],
        axis=0
        )
    
    #Quantif rna2
    if len(rna2_spots) > 0 :
        rna2_result = multistack.extract_cell(
            cell_label=cell_label,
            others_image= {"smfish" : rna2_signal},
            ndim=3,
            nuc_label=nucleus_label,
            rna_coord= rna2_spots,
            others_coord= {'foci' : rna2_cluster}
        )

        for cell in rna2_result :
            rna2_quantif = cell_quant._main_cell_quant(
                cell=cell,
                voxel_size=voxel_size,

            )
            rna2_quantif['cluster_number'] = len(cell['foci'])
            rna2_df = pd.concat([
                rna2_df,
                rna2_quantif
            ],
            axis=0
            )


    #Quantif suntag
    if len(suntag_spots) > 0 :
        suntag_result = multistack.extract_cell(
            cell_label=cell_label,
            ndim=3,
            others_image= {"smfish" : suntag_signal},
            nuc_label=nucleus_label,
            rna_coord= suntag_spots,
            others_coord= {'foci' : suntag_cluster}
        )

        for cell in suntag_result :
            suntag_quantif = cell_quant._main_cell_quant(
                cell=cell,
                voxel_size=voxel_size,

            )
            suntag_quantif['cluster_number'] = len(cell['foci'])
            suntag_df = pd.concat([
                suntag_df,
                suntag_quantif
            ],
            axis=0
            )
    else :
        suntag_df = pd.DataFrame(columns=cell_quant._get_minimal_columns_names())


    #Merge result
    rna1_df['cell_id'] = rna1_df['label']
    rna2_df['cell_id'] = rna2_df['label']
    suntag_df['cell_id'] = suntag_df['label']

    rna1_df = rna1_df.set_index('cell_id', verify_integrity=True)
    rna2_df = rna2_df.set_index('cell_id', verify_integrity=True)
    suntag_df = suntag_df.set_index('cell_id', verify_integrity=True)
    coloc_df = coloc_df.set_index('cell_id', verify_integrity=True)    

    rna1_df = rna1_df.add_prefix('rna1_')
    rna2_df = rna2_df.add_prefix('rna2_')
    suntag_df = suntag_df.add_prefix('suntag_')

    for result in [rna2_df, suntag_df, coloc_df] :
        rna1_df = pd.merge(
            left=rna1_df,
            right=result,
            how='outer',
            left_index=True,
            right_index=True
        )
    
    rna1_df = rna1_df.reset_index(drop=False)
    rna1_df['label'] = rna1_df['cell_id']

    return rna1_df


def _keep_spots_in_mask(spots,mask:np.ndarray) :
    spot_zip = zip(*spots)
    spots_indexer = tuple(list(spot_zip))
    if mask.ndim == 2 and len(spots[0]) == 3 :
        spots_indexer = spots_indexer[1:]

    spots = spots[mask[spots_indexer]]

    return spots

def _cell_coloc_quantification(cell:dict, colocalisation_distance, voxel_size, z_shape) :
    
    #Empty spots list are handled in measure.spots_colocalisation and will return nan.
    
    cell_id = cell['cell_id']
    cell_mask = cell['cell_mask']
    objects = {
        'rna1' : cell.get('rna_coord'),
        'rna2' : cell.get('rna2'),
        'suntag' : cell.get('suntag'),
        'rna1_clustered' : cell.get('rna1_clustered'),
        'rna2_clustered' : cell.get('rna2_clustered'),
        'suntag_clustered' : cell.get('suntag_clustered'),
        'rna1_free' : cell.get('rna1_free'),
        'rna2_free' : cell.get('rna2_free'),
        'suntag_free' : cell.get('suntag_free'),

    }

    #Cleaning spots that are in bounding box but not in segmentation
    for key in objects.keys() :
        spots = objects[key]
        if type(spots) != type(None) : 
            if len(spots) == 0 : continue
            spots = _keep_spots_in_mask(spots, cell_mask)
            objects[key] = spots

    min_y, min_x, max_y, max_x = cell['bbox']
    shape = (z_shape, max_y-min_y, max_x-min_x)

    rna1_number = 0 if type(objects['rna1']) == type(None) else len(objects['rna1'])
    rna2_number = 0 if type(objects['rna2']) == type(None) else len(objects['rna2'])
    suntag_number = 0 if type(objects['suntag']) == type(None) else len(objects['suntag'])
    rna1_clustered_number = 0 if type(objects['rna1_clustered']) == type(None) else len(objects['rna1_clustered'])
    rna2_clustered_number = 0 if type(objects['rna2_clustered']) == type(None) else len(objects['rna2_clustered'])
    suntag_clustered_number = 0 if type(objects['suntag_clustered']) == type(None) else len(objects['suntag_clustered'])
    rna1_free_number = 0 if type(objects['rna1_free']) == type(None) else len(objects['rna1_free'])
    rna2_free_number = 0 if type(objects['rna2_free']) == type(None) else len(objects['rna2_free'])
    suntag_free_number = 0 if type(objects['suntag_free']) == type(None) else len(objects['suntag_free'])


    measures = {
        "cell_id" : [cell['cell_id']],
        "rna1_number" : [rna1_number],
        "rna2_number" : [rna2_number],
        "suntag_number" : [suntag_number],
        "rna1_clustered_number" : [rna1_clustered_number],
        "rna2_clustered_number" : [rna2_clustered_number],
        "suntag_clustered_number" : [suntag_clustered_number],
        "rna1_free_number" : [rna1_free_number],
        "rna2_free_number" : [rna2_free_number],
        "suntag_free_number" : [suntag_free_number],
        }

    for object1, object2 in list(itertools.combinations(objects.keys(), 2)) : # give all the combination of keys of length 2 in a list; all elements are a tuple such as ('rna1', 'rna2')
        
        if object2 in object1 or object1 in object2 : continue
        elif 'rna1' in object1 and 'rna1' in object2 : continue
        elif 'rna2' in object1 and 'rna2' in object2 : continue
        elif 'suntag' in object1 and 'suntag' in object2 : continue
        elif object1 == 'cell_id' or object2 == 'cell_id': continue

        coloc_value_forth = measure.spots_colocalisation(
            shape,
            objects[object1], 
            objects[object2], 
            distance=colocalisation_distance, 
            voxel_size=voxel_size
            )
        
        coloc_value_back = measure.spots_colocalisation(
            shape,
            objects[object2], 
            objects[object1], 
            distance=colocalisation_distance, 
            voxel_size=voxel_size
            )
        
        measures.update({
            '{0}_{1}_colocalisation_count'.format(object1, object2) : [coloc_value_forth],
            '{1}_{0}_colocalisation_count'.format(object1, object2) : [coloc_value_back]
        })
    
    return pd.DataFrame(measures)



def _get_coloc_features_names(objects: list) :

    features_names = ["{0}_number".format(key) for key in objects]
    features_names.append("cell_id")

    for object1, object2 in list(itertools.combinations(objects, 2)) : # give all the combination of keys of length 2 in a list; all elements are a tuple such as ('rna1', 'rna2')
        
        if 'rna1' in object1 and 'rna1' in object2 : continue
        elif 'rna2' in object1 and 'rna2' in object2 : continue
        elif 'suntag' in object1 and 'suntag' in object2 : continue
        elif object1 == 'cell_id' or object2 == 'cell_id': continue

        features_names.append('{0}_{1}_colocalisation_count'.format(object1, object2))

    return features_names


def compute_Colocalisation(Colocalisation: pd.DataFrame, AcquisitionId: int, spots1: list, spots2: list, suntag_spots:list, rna1_clusterd_spots: pd.DataFrame, rna2_clusterd_spots: pd.DataFrame, image_shape, colocalisation_distance, voxel_size) :

    """
    Computes and add new colocalisation measurements to Colocalisation dataframe.

    Takes as parameters 2 lists of spots plus a third (optional) one.

    Parameters
    ----------
    AcquisitionId : int
    id : int
    spots1 : list[np.array]
        list of spots which we are trying to see colocalize around protein of interest. In other words those are the spots that are counted.
    spots2 : list[np.array]
        list of spots which are gathering spots around them. In other words, these spots are NOT counted, they are used as the condition to count spots from 'spots1' argument.
    image_shape : tuple(z;y;x)
    colocalisation_distance : int, float
        distance in nanometer up to which spots are considered as colocalising.
    voxel_size : tuple(z;y;x)
        scale factors px <-> nm.
    
    Returns
    -------
    Colocalisation_df : pd.DataFrame
        Result table.
        Keys : 'id', 'AcquisitionId', 'rna1_number',  'rna2_number', 'rna1_rna2_colocalisation_count', 'rna2_rna1_colocalisation_count'
            /!\ 'id' is equivalent to 'AcquisitionId', since there is (1,1) relation between those table I decided to have 2 keys with those names for potential data manipulation later on.  
    """

    if not isinstance(Colocalisation, pd.DataFrame) : raise TypeError("Colocalisation should be a DataFrame it is a {0}".format(type(Colocalisation)))
    if not isinstance(AcquisitionId, (int, np.int64, np.int0, np.int8, np.int16, np.int32)) : raise TypeError("AcquisitionId should be a int it is a {0}".format(type(AcquisitionId)))
    if not isinstance(spots1, (list, np.ndarray)) : raise TypeError("spots1 should be a list it is a {0}".format(type(spots1)))
    if not isinstance(spots2, (list, np.ndarray)) : raise TypeError("spots2 should be a list it is a {0}".format(type(spots2)))
    if not isinstance(suntag_spots, (list, np.ndarray)) : raise TypeError("suntag should be a list it is a {0}".format(type(suntag_spots)))
    if not isinstance(rna1_clusterd_spots, pd.DataFrame) : raise TypeError("rna1_clusterd_spots should be a DataFrame it is a {0}".format(type(rna1_clusterd_spots)))
    if not isinstance(rna2_clusterd_spots, pd.DataFrame) : raise TypeError("rna2_clusterd_spots should be a DataFrame it is a {0}".format(type(rna2_clusterd_spots)))
    if not isinstance(colocalisation_distance, (float, int)) : raise TypeError("colocalisation_distance should be a int/float it is a {0}".format(type(colocalisation_distance)))
    if not isinstance(image_shape, (tuple, list)) : raise TypeError("image_shape should be a tuple/list it is a {0}".format(type(image_shape)))
    if not isinstance(voxel_size, (tuple, list)) : raise TypeError("voxel_size should be a tuple/list it is a {0}".format(type(voxel_size)))

    if len(voxel_size) != 3 : raise ValueError("voxel_size should be given in 3D : (z,y,x)")
    
    rna1_number = len(spots1)
    if rna1_number == 0 : rna1_number = np.NaN
    rna2_number = len(spots2)
    if rna2_number == 0 : rna2_number = np.NaN
    suntag_number = len(suntag_spots)
    if suntag_number == 0 : suntag_number = np.NaN
    rna1_rna2_colocalisation_count = measure.spots_colocalisation(image_shape, spots1, spots2, distance=colocalisation_distance, voxel_size=voxel_size)
    rna2_rna1_colocalisation_count = measure.spots_colocalisation(image_shape, spots2, spots1, distance=colocalisation_distance, voxel_size=voxel_size)
    rna1_number_1000nm_around_rna2 = measure.spots_colocalisation(image_shape, spots1, spots2, distance=1000, voxel_size=voxel_size)
    rna1_number_2000nm_around_rna2 = measure.spots_colocalisation(image_shape, spots1, spots2, distance=2000, voxel_size=voxel_size)
    rna1_suntag_colocalisation_count = measure.spots_colocalisation(image_shape, spots1, suntag_spots, distance=colocalisation_distance, voxel_size=voxel_size)
    rna2_suntag_colocalisation_count = measure.spots_colocalisation(image_shape, spots2, suntag_spots, distance=colocalisation_distance, voxel_size=voxel_size)
    rna1_cluster_number = rna1_clusterd_spots['cluster_id'].count()
    rna2_cluster_number = rna2_clusterd_spots['cluster_id'].count()

    Colocalisation_df = pd.DataFrame({
        'id' : [AcquisitionId],
        'AcquisitionId' : [AcquisitionId],
        'rna1_number' : [rna1_number],
        'rna2_number' : [rna2_number],
        'suntag_number' : [suntag_number],
        'rna1_cluster_number' : [rna1_cluster_number],
        'rna2_cluster_number' : [rna2_cluster_number],
        'rna1_rna2_colocalisation_count' : [rna1_rna2_colocalisation_count],
        'rna2_rna1_colocalisation_count' : [rna2_rna1_colocalisation_count],
        'rna1_number_1000nm_around_rna2' : [rna1_number_1000nm_around_rna2],
        'rna1_number_2000nm_around_rna2' : [rna1_number_2000nm_around_rna2],
        'rna1_suntag_colocalisation_count' : [rna1_suntag_colocalisation_count],
        'rna2_suntag_colocalisation_count' : [rna2_suntag_colocalisation_count],

    })
    Colocalisation_df = pd.concat([Colocalisation, Colocalisation_df], axis= 0)
    return Colocalisation_df


def compute_ClusteredSpots(ClusteredSpots,
                            AcquisitionId,
                            rna_clustered_spots_dataframe : pd.DataFrame,
                            sun_tag_spots,
                            radius_nm,
                            image_shape,
                            voxel_size,
                            rna_name:str,
                            cell_mask=None
                            ) :
    """
    Add quantification to ClusteredSpots dataframe.


    Parameters
    ----------
        ClusteredSpots : pd.DataFrame
        AcquisitionId : int
        rna_clustered_spots_dataframe : pd.DataFrame
        suntag_spots : np.ndarray(shape=(spots_num,3))
        radius_nm : float
        image_shape : np.ndarray, tuple
        voxel_size : tuple
            (z,y,x)
        rna_name: str
        cell_mask : np.ndarray

    """
    
    if not isinstance(rna_clustered_spots_dataframe, pd.DataFrame) : raise TypeError("rna2_clustered_spots_dataframe should be pd.DataFrame it is a {0}".format(type(rna_clustered_spots_dataframe)))

    spots_list = get_centroids_list(rna_clustered_spots_dataframe)
    random_spots = _random_spot_list(image_shape, spot_number= len(spots_list), mask=cell_mask)

    colocalising_suntag_number_list = measure.spots_multicolocalisation(spots_list, sun_tag_spots, radius_nm=radius_nm, image_shape=image_shape, voxel_size=voxel_size)
    clostest_suntag_distance_list = measure.closest_spot_distance(spots_list, sun_tag_spots, image_shape=image_shape, voxel_size=voxel_size)
    clostest_suntag_random_distance_list = measure.closest_spot_distance(random_spots, sun_tag_spots, image_shape=image_shape, voxel_size=voxel_size)

    median_suntag_distance = np.median(clostest_suntag_random_distance_list)
    mean_suntag_distance = np.mean(clostest_suntag_random_distance_list)
    del clostest_suntag_random_distance_list, random_spots

    df = rna_clustered_spots_dataframe.copy()
    df["AcquisitionId"] = [AcquisitionId] * len(rna_clustered_spots_dataframe)
    df["colocalising_suntag_number"] = np.NaN
    df["colocalising_suntag_number"] = colocalising_suntag_number_list
    df["distance_closest_suntag"] = clostest_suntag_distance_list
    df["median_suntag_distance"] = median_suntag_distance
    df["mean_suntag_distance"] = mean_suntag_distance
    df['rna'] = rna_name

    df = pd.concat([ClusteredSpots, df], axis= 0)
    return df


def _random_spot_list_no_mask(shape, spot_number, off_set_list=None) :

    if type(off_set_list) == type(None) : off_set_list = (0,) * len(shape)
    if len(shape) != len(off_set_list) : raise ValueError("shape and off_set_list must have the same length.")

    random_generator = np.random.default_rng()
    Z,Y,X = [random_generator.integers(
        low= 0 + off_set,
        high= axis + off_set,
        size= spot_number
        ) for axis, off_set in zip (shape, off_set_list)]
    
    return list(zip(Z,Y,X))


def _random_spot_list_from_mask(shape, spot_number, mask: np.ndarray) :
    if mask.dtype != bool : 
        mask = mask.astype(bool)
        print("Warning mask dtype was not bool. Mask as been cast to bool.")
    
    coordinates_grid = np.indices(shape)
    if len(shape) == 3 :
        Z, Y , X = coordinates_grid[0][mask], coordinates_grid[1][mask], coordinates_grid[2][mask]
        spot_list = np.array(list(zip(Z,Y,X)), dtype= int)
    elif len(shape) == 2:
        Y, X = coordinates_grid[0][mask], coordinates_grid[1][mask]
        spot_list = np.array(list(zip(Y,X)), dtype= int)
    else : raise ValueError("Shape should either be 3d or 2d. It is {0}d".format(len(shape)))
    del coordinates_grid
    
    if len(spot_list) == 0 : return(np.empty(shape=(0,0), dtype= int))
    random_generator = np.random.default_rng()
    random_spot_index = random_generator.integers(0, len(spot_list), size= spot_number)

    return spot_list[random_spot_index]


def _random_spot_list(shape, spot_number, off_set_list=None, mask=None) :

    if spot_number == 0 :
        return np.empty(shape=(0,0), dtype= int)

    if type(mask) != type(None) and type(off_set_list) != type(None) :
        raise ValueError("off_set_list and mask parameters cannot be used simultaneously please either one to None.")
    
    elif type(mask) != type(None) and type(off_set_list) == type(None) :
        return _random_spot_list_from_mask(shape, spot_number,mask)
    
    else :
        return _random_spot_list_no_mask(shape, spot_number, off_set_list)