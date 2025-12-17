"""
Module to perform multithreaded spot detection with all post processing.
"""
import pandas as pd
import numpy as np

from .wrapper import full_detection, cell_quantification
from numpy import NaN

def multi_thread_full_detection(
        image,
        voxel_size,
        threshold,
        spot_radius,
        alpha,
        beta,
        gamma,
        artifact_size,
        cluster_radius,
        min_spot_per_cluster,
        filename:str,
        detection_id,
) :

    try :
        res = full_detection(image,voxel_size,threshold,spot_radius,alpha,beta,gamma,artifact_size,cluster_radius,min_spot_per_cluster,filename,detection_id,)
    
    except Exception as error :
        print("Thread {0} : Error occured\n{1}".format(detection_id, str(error)))
        keys = ['spots','spots_post_decomp','clustered_spots_dataframe','clusters_dataframe','clusters','clustered_spots','free_spots','threshold','voxel_size','spot_radius','alpha','beta','gamma','artifact_size','cluster_radius','min_spot_per_cluster',]
        res = dict.fromkeys(keys, NaN)
        res['detection_id'] = detection_id
        print("Thread {0} : Returning NaN values".format(detection_id))

    finally :
        return res




def build_Spots_and_Cluster_df(detection_result : dict) :
    """
    Made to build Spots pandas dataframe from result of multi_threaded call to `multi_thread_full_detection`.

    Parameter
    ---------
        detection_result
    """

    SPOTS_COLUMNS = [
        'detection_id',
        'spots',
        'spots_post_decomp',
        'clustered_spots_dataframe',
        'clusters_dataframe',
        'clusters',
        'clustered_spots',
        'free_spots',   
    ]
    
    detection_res = pd.DataFrame(columns=SPOTS_COLUMNS, data=detection_result)
    detection_res = detection_res.set_index('detection_id', verify_integrity=True, drop=False)

    Spots = pd.DataFrame()
    Clusters = pd.DataFrame()

    for detection_id in detection_res.index :
        spots: pd.DataFrame = detection_res.at[detection_id, 'clustered_spots_dataframe']
        clusters: pd.DataFrame = detection_res.at[detection_id, 'clusters_dataframe']

        #names
        spots = spots.rename(columns={'id' : 'spot_id'})
        clusters = clusters.rename(columns={'id' : 'cluster_id'})
        
        #is_clustered
        spots['cluster_id'] = spots['cluster_id'].replace(-1,np.NaN)
        clusters['cluster_id'] = clusters['cluster_id'].replace(-1,np.NaN)
        spots['population'] = spots['cluster_id'].isna()
        spots['population'] = spots['population'].replace({True : 'free', False : 'clustered'})

        #Detection_id
        spots['detection_id'] = detection_id
        clusters['detection_id'] = detection_id

        Spots = pd.concat([
            Spots,
            spots,
        ],axis=0)

        Clusters = pd.concat([
            Clusters,
            clusters,
        ],axis=0)


    Spots = Spots.reset_index(drop=True)
    Clusters = Clusters.reset_index(drop=True)

    return Spots, Clusters


def multi_thread_cell_quantification(
        acquisition_id,
        detection_id,
        spots,
        clusters,
        voxel_size,
        cell_label,
        nucleus_label,
        fish_signal_3D,
        dapi_signal_3D,
        ndim=3,
) :
    print("Thread is starting detection {0}".format(detection_id))

    try :
        res = cell_quantification(
        acquisition_id,
        detection_id,
        spots,
        clusters,
        voxel_size,
        cell_label,
        nucleus_label,
        fish_signal_3D,
        dapi_signal_3D,
        ndim=3
        )
    
    except Exception as error :
        print("Thread {0} : Error occured\n{1}".format(detection_id, str(error)))
        res = pd.DataFrame()
        print("Thread {0} : Returning empty df".format(detection_id))

    finally :
        print("Thread {0} : finishing".format(detection_id))
        return res