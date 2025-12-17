
"""
Spot quantification for Sequential Fish data

This script use results from FishSeq_pipeline_segmentation.py that must be run before
"""

import os
import numpy as np
import dataprocessing.fishseq_preprocess_folder as prepro
import bigfish.stack as stack
import bigfish.plot as plot
import bigfish.detection as bigfish
import pbwrap.detection as detection
from pbwrap.utils import open_image
from CustomPandasFramework.utils import get_datetime
from tqdm import tqdm

#########
## USER PARAMETERS
#########

RUN_PATH = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2024-03-05 - R2TP_Pool-1_Run1/"
OUTPUT_PATH = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/output"
DATASET_NAME = "R2TP_Pool-1_Run1"
from input_config.FishSeq_input_config import THRESHOLD_DICT

voxel_size = (300,103,103)
spot_size = (150,100,100)

alpha = 0.5
beta = 1
sigma = 0.3

cluster_size = 400 #nm
min_spot_per_cluster = 4

artifact_radius = 1400

# Main Script

#Loading data
file_dict = prepro.assert_run_folder_integrity(run_path=RUN_PATH)
location_list = list(file_dict.keys())
threshold_dict = THRESHOLD_DICT[DATASET_NAME]

#preparing folders
save_path = OUTPUT_PATH + "{0}_{1}/".format(get_datetime(), DATASET_NAME)
result_path = save_path + "result_tables/"
detection_path = save_path + "detection/"
os.makedirs(save_path)
os.makedirs(result_path)
os.makedirs(detection_path)

#Main loop
for location_id, location in enumerate(location_list) :

    filename = file_dict[location][0]
    image_path = RUN_PATH + "/Cy3_Z-stacks/{0}/{1}".format(location, filename) 
    multichannel_stack = open_image(image_path)# This open 4D multichannel image (all the images are loaded in one call)
    image_number = len(multichannel_stack)
    id_list = range(image_number)

    #Detection
    for detection_id in tqdm(id_list) :
        im_stack = multichannel_stack[detection_id]
        threshold = threshold_dict.get(detection_id)
        
        spots = bigfish.detect_spots(
            images= im_stack,
            threshold=threshold,
            return_threshold=False,
            voxel_size=voxel_size,
            spot_radius=spot_size,
        )

        spots_post_decomp = detection.cluster_deconvolution(
            image=im_stack,
            spots=spots,
            spot_radius=spot_size,
            alpha=alpha,
            beta= beta,
            sigma=sigma,
            timer=300 #s
        )

        spots_post_decomp = detection.remove_artifact(
            spots_post_decomp, 
            artifact_radius=artifact_radius, 
            voxel_size=voxel_size
            )

        spots_cluster_results= detection.cluster_detection(
            spots_post_decomp,
            voxel_size= voxel_size, 
            nb_min_spots=min_spot_per_cluster,
            radius=cluster_size, 
            keys_to_compute=['clustered_spots_dataframe', 'clusters_dataframe', 'clusters']
            )
        
        clustered_spots_dataframe, clusters_dataframe, clusters = spots_cluster_results['clustered_spots_dataframe'], spots_cluster_results['clusters_dataframe'], spots_cluster_results['clusters']
        clustered_spots = detection.get_centroids_array(clustered_spots_dataframe.query("not cluster_id.isna()"))
        free_spots = detection.get_centroids_array(clustered_spots_dataframe.query("cluster_id.isna()"))
        plot.output_spot_tiffvisual(np.max(im_stack, axis=0), [spots, spots_post_decomp, clustered_spots, free_spots], detection_path + "{0}_{1}_detection.tif".format(filename, detection_id), dot_size= 2)
        
        
        
    


