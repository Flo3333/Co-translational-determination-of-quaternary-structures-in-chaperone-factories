
import numpy as np
import pandas as pd
import pbwrap.plot as plot
import bigfish.detection as bigfish
import bigfish.multistack as multistack
import pbwrap.quantification.hub as cell_quant
from .clusters import cluster_detection, remove_artifact, get_centroids_array
from .bigfish_wrapper import cluster_deconvolution
from pbwrap.utils import _keep_spots_in_mask
from tqdm import tqdm


def full_detection (
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
    """
    Performs full detection process : thread-compatible.

    1. detect spot (bigfish.detection)
    2. cluster_deconvolution (pbwrap.detection)
    3. remove_artificats (pbwrap.detection) (skipped if artificat_size is None)
    4. cluster_detection
    5. plot detection visual (skipped if filename is None)

    RETURN
    ------
    res : dict
        {
        'id'
        'spots'
        'spots_post_decomp'
        'clustered_spots_dataframe',
        'clusters_dataframe'
        'clusters'
        'clustered_spots'
        'free_spots'
        }


    """
    spots, threshold = bigfish.detect_spots(
            images= image,
            threshold=threshold,
            return_threshold = True,
            voxel_size=voxel_size,
            spot_radius=spot_radius,
        )

    spots_post_decomp = cluster_deconvolution(
        image=image,
        spots=spots,
        spot_radius=spot_radius,
        voxel_size=voxel_size,
        alpha=alpha,
        beta= beta,
        sigma=gamma,
        timer=0 #s
    )

    if type(artifact_size) != type(None) :
        spots_post_decomp = remove_artifact(
            spots_post_decomp, 
            artifact_radius=artifact_size, 
            voxel_size=voxel_size
            )

    spots_cluster_results= cluster_detection(
        spots_post_decomp,
        voxel_size= voxel_size, 
        nb_min_spots=min_spot_per_cluster,
        radius=cluster_radius, 
        keys_to_compute=['clustered_spots_dataframe', 'clusters_dataframe', 'clusters']
        )
        
    clustered_spots_dataframe, clusters_dataframe, clusters = spots_cluster_results['clustered_spots_dataframe'], spots_cluster_results['clusters_dataframe'], spots_cluster_results['clusters']
    clustered_spots = get_centroids_array(clustered_spots_dataframe.query("not cluster_id.isna()"))
    free_spots = get_centroids_array(clustered_spots_dataframe.query("cluster_id.isna()"))
    if type(filename) != type(None) :
        plot.output_spot_tiffvisual(np.max(image, axis=0), [spots, spots_post_decomp, clustered_spots, free_spots], filename, dot_size= 2)
        filename = filename.split('/')[-1]

    clustered_spots_dataframe['intensity'] = image[
        list(clustered_spots_dataframe['z']),
        list(clustered_spots_dataframe['y']),
        list(clustered_spots_dataframe['x']),
    ]

    res = {
        'detection_id' : detection_id,
        'spots' : spots,
        'spots_post_decomp' : spots_post_decomp,
        'clustered_spots_dataframe' : clustered_spots_dataframe,
        'clusters_dataframe' : clusters_dataframe,
        'clusters' : clusters,
        'clustered_spots' : clustered_spots,
        'free_spots' : free_spots,
    }
    

    return res


def cell_quantification(
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
        talk=False,
        
) :
    
    if len(spots) == 0 :
        spots = np.empty(shape=(0,3), dtype=int)
    if len(clusters) == 0 :
        clusters = np.empty(shape=(0,5),dtype=int)

    cell_extraction = multistack.extract_cell(
            cell_label=cell_label,
            others_image= {"smfish" : fish_signal_3D},
            ndim=ndim,
            nuc_label=nucleus_label,
            rna_coord= spots,
            others_coord= {'foci' : clusters}
        )
    
    cell_df = pd.DataFrame()
    for cell in cell_extraction :
        new_cell = cell_quant._main_cell_quant(
            cell=cell,
            voxel_size=voxel_size,
            dapi_stack=dapi_signal_3D,
        )

        new_cell['acquisition_id'] = acquisition_id
        new_cell['detection_id'] = detection_id
        new_cell['cluster_number'] = len(cell['foci'])
        cell_spots = cell['rna_coord']
        if len(cell_spots) > 0 : cell_spots = _keep_spots_in_mask(cell_spots, cell['cell_mask'])
        new_cell['rna_number'] = len(cell_spots)

        cell_df = pd.concat([
            cell_df,
            new_cell
        ],axis=0)

    if talk : print("detection {0} done.".format(detection_id))
    return cell_df