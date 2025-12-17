import os
import numpy as np
import pandas as pd
import time, signal

import CustomPandasFramework.log as logs
from CustomPandasFramework.computer_interface import get_datetime
import CustomPandasFramework.Fish_seq.framework as framework
import CustomPandasFramework.Fish_seq.quantification as quantification
import CustomPandasFramework.Fish_seq.run_integrity as integrity

import bigfish.stack as stack
import bigfish.plot as plot
from bigfish.multistack import match_nuc_cell

import pbwrap.segmentation as segmentation
import pbwrap.plot as plt
import pbwrap.detection as detection
from pbwrap.integrity import detectiontimeout_handler

DATASET = "240429_time_course_CB" # For a new dataset add config first in input_config.py
MULTICHANNEL = True

path_in = "/media/floricslimani/SSD4To/SSD_floricslimani/Fish_seq/input/"
path_out = "/media/floricslimani/SSD4To/SSD_floricslimani/Fish_seq/output/"

#### PARAMETERS ####
#channel from 0 to 3
rna_spot_radius = (300,100,100)
suntag_spot_radius = (300,100,100)
voxel_size = (240, 103, 103) #nm

#Detection
#Auto
threshold_penalty_rna1 = 1.2
threshold_penalty_rna2 = 1.2
threshold_penalty_suntag = 0.4

#Manual overall setting taking priority if != None
THRESHOLD_RNA1 = None
THRESHOLD_RNA2 = None
THRESHOLD_SUNTAG = None

#Artifact_correction
REMOVE_ARTIFACTS = True
ARTIFACT_RADIUS = 1400


alpha = 0.5
beta = 1
sigma = 3
_crop_pixel = 0
colocalisation_distance = 310 # ~ 3 pixels
rna_cluster_radius = 500
min_spot_number_per_cluster_rna1 = 5
min_spot_number_per_cluster_rna2 = 5
min_spot_number_per_cluster_suntag = 5
sample_size = 8

#Segmentation
nucleus_model = 'nuclei'
nucleus_size = 180
cytoplasm_model = 'cyto3'
cytoplasm_size = 250

do_detection_visuals = True # Set True if you want to check the dense region decomposition or the cluster detection else colocalisation_visuals are enough.
do_colocalisation_visuals = True

#CONSTANT
EMPTY_SPOTS = np.empty(shape=(0,0), dtype=int) 


#########
#Pipeline
clock = time.process_time()
from input_config.BC_input_config import REGEX, CHANNELS, COLUMNS, GROUP_KEYS, THRESHOLDS
regex = REGEX[DATASET]
columns = COLUMNS[DATASET]
channels = CHANNELS[DATASET]
group_keys = GROUP_KEYS[DATASET]
thresholds: dict = THRESHOLDS.setdefault(DATASET, dict()) #Above overall thresholds takes priority

Input = framework.create_Input(path_in, regex=regex, columns=columns, multichannel = MULTICHANNEL)
analysis_group = framework.make_analysis_group(Input=Input, group_keys= group_keys)
date = get_datetime()
path_out += DATASET + '_' + date + '/'
os.makedirs(path_out, exist_ok= True)
signal.signal(signal.SIGALRM, detectiontimeout_handler) #Initiating timeout handling

# Log
parameters_log = logs.parameter_log("parameters.txt", path= path_out)
parameters_log.add_parameters(DATASET, MULTICHANNEL, rna_spot_radius, suntag_spot_radius, voxel_size, THRESHOLD_RNA1, THRESHOLD_RNA2, THRESHOLD_SUNTAG, threshold_penalty_rna1, threshold_penalty_rna2, threshold_penalty_suntag,  alpha,  beta, _crop_pixel, colocalisation_distance, colocalisation_distance, rna_cluster_radius, min_spot_number_per_cluster_rna1, min_spot_number_per_cluster_rna2, min_spot_number_per_cluster_suntag, nucleus_model, nucleus_size, cytoplasm_model, cytoplasm_size)
parameters_log.write()
Acquisition = framework.create_Acquisiton(Input)
print(Acquisition.value_counts(subset=group_keys))
Cell = pd.DataFrame(columns=['cell_id'])
Spots = pd.DataFrame(columns=['cell_id','spots_id'])
for group_index in analysis_group.index :

    rna1_threshold, rna2_threshold, suntag_threshold = thresholds.setdefault(group_index, (None,None,None))

    group = analysis_group.loc[group_index]
    print('\033[33m' + "\nNext group : {0}".format(group_index) + '\033[00m')
    output_path = path_out + "/{0}/".format(group_index)
    os.makedirs(output_path, exist_ok= True)

    if 'rna2' in Input.columns :
        has_rna2 = not any(
            Input.loc[group]['rna2'].isna()
        )

    else : has_rna2 = False

    if 'suntag' in Input.columns :
        has_suntag = not any(
            Input.loc[group]['suntag'].isna()
        )
    else :
        has_suntag = False


    if type(THRESHOLD_RNA1) != type(None):
        rna1_threshold = THRESHOLD_RNA1
    elif type(rna1_threshold) != type(None) :
        pass 
    else : 
        print("computing threshold for rna1...")
        rna1_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, multichannel=MULTICHANNEL, index_list= group, channels= channels['rna1'], _crop_pixel= _crop_pixel)
        rna1_threshold = detection.compute_auto_threshold(rna1_gen, voxel_size=voxel_size, spot_radius=rna_spot_radius, im_number=sample_size) * threshold_penalty_rna1
    print('thresholrd rna1 : ', rna1_threshold)

    if has_rna2 :
        if type(THRESHOLD_RNA2) != type(None) : 
            rna2_threshold = THRESHOLD_RNA2
        elif type(rna2_threshold) != type(None) :
            pass
        else : 
            print("computing threshold for rna2...")
            rna2_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, multichannel=MULTICHANNEL, index_list= group, channels= channels['rna2'], _crop_pixel= _crop_pixel)
            rna2_threshold = detection.compute_auto_threshold(rna2_gen, voxel_size=voxel_size, spot_radius=rna_spot_radius, im_number=sample_size) * threshold_penalty_rna2
        print('thresholrd rna2 : ', rna2_threshold)

    if has_suntag :
        print("computing threshold for suntag...")
        if type(THRESHOLD_SUNTAG) != type(None) : 
            suntag_threshold = THRESHOLD_SUNTAG
        elif type(suntag_threshold) != type(None) :
            pass
        else : 
            suntag_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, multichannel=MULTICHANNEL, index_list= group, channels= channels['suntag'], _crop_pixel= _crop_pixel)
            suntag_threshold = detection.compute_auto_threshold(suntag_gen, voxel_size=voxel_size, spot_radius=rna_spot_radius, im_number=sample_size) * threshold_penalty_suntag
        print('thresholrd suntag : ', suntag_threshold)

    if MULTICHANNEL : 
        im_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, multichannel=MULTICHANNEL, index_list= group, _crop_pixel= _crop_pixel)
    else : 
        dapi_gen = framework.get_image_as_gen(Input=Input, path_in=path_in, multichannel=MULTICHANNEL, index_list=group, channels= channels["dapi"], _crop_pixel=_crop_pixel)
        rna1_gen = framework.get_image_as_gen(Input=Input, path_in=path_in, multichannel=MULTICHANNEL, index_list=group, channels= channels["rna1"], _crop_pixel=_crop_pixel)
        if has_rna2 : rna2_gen = framework.get_image_as_gen(Input=Input, path_in=path_in, multichannel=MULTICHANNEL, index_list=group, channels= channels["rna2"], _crop_pixel=_crop_pixel)
        if has_suntag : suntag_gen = framework.get_image_as_gen(Input=Input, path_in=path_in, multichannel=MULTICHANNEL, index_list=group, channels= channels["suntag"], _crop_pixel=_crop_pixel)

    for input_id in group :
        
        filename = Input.at[input_id, 'filename'].replace('czi','')
        rna1, rna2 = Input.at[input_id, 'rna1'], Input.at[input_id, 'rna2']
        index = Input.at[input_id, 'id']
        bool_rna1_decomp = True
        bool_rna2_decomp = True

        print("Computing fov {0}".format(filename))

        if MULTICHANNEL :
            im = next(im_gen)
            dapi_chan = CHANNELS.get(DATASET)['dapi']
            rna1_chan = CHANNELS.get(DATASET)['rna1']
            if has_rna2 : rna2_chan = CHANNELS.get(DATASET)['rna2']
            if has_suntag : suntag_chan = CHANNELS.get(DATASET)['suntag']

            dapi = im[dapi_chan]
            rna1 = im[rna1_chan]
            if has_rna2 : rna2 = im[rna2_chan]
            if has_suntag : suntag = im[suntag_chan]
        else :
            dapi = next(dapi_gen)
            rna1 = next(rna1_gen)
            if has_rna2 : rna2 = next(rna2_gen)
            if has_suntag : suntag = next(suntag_gen)

        #detection
        rna1_spots = detection.detect_spots(rna1, threshold= rna1_threshold, ndim= 3, voxel_size= voxel_size, spot_radius= rna_spot_radius)
        if has_rna2 : rna2_spots = detection.detect_spots(rna2, threshold= rna2_threshold, ndim= 3, voxel_size= voxel_size, spot_radius= rna_spot_radius)
        if has_suntag : suntag_spots = detection.detect_spots(suntag, threshold= suntag_threshold, ndim= 3, voxel_size= voxel_size, spot_radius= rna_spot_radius)

        #image projection
        nucleus = stack.mean_projection(dapi)
        cyto = stack.mean_projection(rna1)
        shape = nucleus.shape
        rna1_max = stack.maximum_projection(rna1)
        if has_rna2 : rna2_max = stack.maximum_projection(rna2)
        if has_suntag : suntag_max = stack.maximum_projection(suntag)

        #Segmentation
        print(" Segmenting...")
        nucleus_label = segmentation.Nucleus_segmentation(nucleus, diameter= nucleus_size, use_gpu= True, model_type= nucleus_model)
        cell_label = segmentation.Cytoplasm_segmentation(cyto, nucleus, diameter= cytoplasm_size, use_gpu= True, model_type= cytoplasm_model)
        nucleus_label, cell_label = match_nuc_cell(nuc_label=nucleus_label, cell_label= cell_label, single_nuc=True, cell_alone=False)

        print(" Spots clusters deconvolution...")
        rna1_spots_postdecomp = detection.cluster_deconvolution(rna1, rna1_spots, rna_spot_radius, voxel_size=voxel_size, alpha= alpha, beta= beta, sigma=sigma, timer= 60)
        
        if has_rna2 : 
            rna2_spots_postdecomp = detection.cluster_deconvolution(rna2, rna2_spots, rna_spot_radius, voxel_size=voxel_size, alpha= alpha, beta= beta, sigma= sigma, timer= 60)
        else :
            rna2_spots_postdecomp = EMPTY_SPOTS

        if has_suntag : 
            suntag_spots_postdecomp = detection.cluster_deconvolution(suntag, suntag_spots, suntag_spot_radius, voxel_size=voxel_size, alpha= alpha, beta= beta, sigma= sigma, timer= 60)
        else :
            suntag_spots_postdecomp = EMPTY_SPOTS

        if REMOVE_ARTIFACTS:
            print(" Removing artifacts...")
            rna1_spots_postdecomp = detection.remove_artifact(rna1_spots_postdecomp, artifact_radius=ARTIFACT_RADIUS, voxel_size=voxel_size)
            if has_rna2 : rna2_spots_postdecomp = detection.remove_artifact(rna2_spots_postdecomp, artifact_radius=ARTIFACT_RADIUS, voxel_size=voxel_size)
            if has_suntag : suntag_spots_postdecomp = detection.remove_artifact(suntag_spots_postdecomp, artifact_radius=ARTIFACT_RADIUS, voxel_size=voxel_size)


        print(" Clusters computation...")
        rna1_cluster_results= detection.cluster_detection(rna1_spots_postdecomp, voxel_size= voxel_size, nb_min_spots=min_spot_number_per_cluster_rna1,radius=rna_cluster_radius, keys_to_compute=['clustered_spots_dataframe', 'clusters_dataframe', 'clusters'])
        rna1_clustered_spots_dataframe, rna1_clusters_dataframe= rna1_cluster_results['clustered_spots_dataframe'], rna1_cluster_results['clusters_dataframe']
        
        rna2_cluster_results= detection.cluster_detection(rna2_spots_postdecomp, voxel_size= voxel_size, nb_min_spots=min_spot_number_per_cluster_rna2,radius=rna_cluster_radius, keys_to_compute=['clustered_spots_dataframe', 'clusters_dataframe','clusters'])
        rna2_clustered_spots_dataframe, rna2_clusters_dataframe= rna2_cluster_results['clustered_spots_dataframe'], rna2_cluster_results['clusters_dataframe']
        
        suntag_cluster_results= detection.cluster_detection(suntag_spots_postdecomp, voxel_size= voxel_size, nb_min_spots=min_spot_number_per_cluster_suntag,radius=rna_cluster_radius, keys_to_compute=['clustered_spots_dataframe', 'clusters_dataframe','clusters'])
        suntag_clustered_spots_dataframe, suntag_clusters_dataframe= suntag_cluster_results['clustered_spots_dataframe'], suntag_cluster_results['clusters_dataframe']

        print(" Plotting check visuals...")
        rna1_cluster_spots = detection.get_centroids_array(rna1_clustered_spots_dataframe.query("not cluster_id.isna()"))
        rna1_uncluster_spots = detection.get_centroids_array(rna1_clustered_spots_dataframe.query("cluster_id.isna()"))
        rna2_cluster_spots = detection.get_centroids_array(rna2_clustered_spots_dataframe.query("not cluster_id.isna()"))
        rna2_uncluster_spots = detection.get_centroids_array(rna2_clustered_spots_dataframe.query("cluster_id.isna()"))
        suntag_cluster_spots = detection.get_centroids_array(suntag_clustered_spots_dataframe.query("not cluster_id.isna()"))
        suntag_uncluster_spots = detection.get_centroids_array(suntag_clustered_spots_dataframe.query("cluster_id.isna()"))
        rna1_cluster = detection.get_centroids_array(rna1_clusters_dataframe)
        rna2_cluster = detection.get_centroids_array(rna2_clusters_dataframe)
        suntag_cluster = detection.get_centroids_array(suntag_clusters_dataframe)

        if do_detection_visuals :
            plot.plot_segmentation_boundary(cyto, cell_label, nucleus_label, boundary_size=3, rescale=True, path_output= output_path + "{0}_segmentation".format(filename), show=False)
            plt.output_spot_tiffvisual(rna1_max, [rna1_spots, rna1_spots_postdecomp, rna1_cluster_spots, rna1_uncluster_spots], output_path + "{0}_rna1_detection.tif".format(filename), dot_size= 2)
            if has_rna2 : plt.output_spot_tiffvisual(rna2_max, [rna2_spots, rna2_spots_postdecomp , rna2_cluster_spots, rna2_uncluster_spots], output_path + "{0}_rna2_detection.tif".format(filename), dot_size= 2)
            if has_suntag : plt.output_spot_tiffvisual(suntag_max, [suntag_spots, suntag_spots_postdecomp , suntag_cluster_spots, suntag_uncluster_spots], output_path + "{0}_suntag_detection.tif".format(filename), dot_size= 2)

        if len(rna1_spots_postdecomp) == 0 : 
            print("rna1 detection yielded 0 spots; computing next fov...")
            continue
        elif (len(rna2_spots_postdecomp) == 0 and has_rna2) : 
            print("rna2 detection yielded 0 spots; computing next fov...")
            continue
        elif(len(suntag_spots_postdecomp)==0 and has_suntag) :
            print("suntag detection yielded 0 spots; computing next fov...")
            continue

        fov_cell = quantification.cell_quantif( #cell_id is cell label value
            cell_label=cell_label,
            nucleus_label= nucleus_label,
            voxel_size=voxel_size,
            nucleus_signal= dapi,
            colocalisation_distance=colocalisation_distance,
            rna1_spots= rna1_spots_postdecomp,
            rna2_spots= rna2_spots_postdecomp,
            suntag_spots = suntag_spots_postdecomp,
            rna1_cluster= rna1_cluster_results['clusters'],
            rna2_cluster= rna2_cluster_results['clusters'],
            suntag_cluster= suntag_cluster_results['clusters'],
            rna1_clustered_spots= rna1_cluster_spots,
            rna2_clustered_spots= rna2_cluster_spots,
            suntag_clustered_spots= suntag_cluster_spots,
            rna1_free_spots= rna1_uncluster_spots,
            rna2_free_spots= rna2_uncluster_spots,
            suntag_free_spots= suntag_uncluster_spots,
            rna1_signal = rna1_max,
            rna2_signal = rna2_max if has_rna2 else np.empty(shape=rna1.shape),
            suntag_signal = suntag_max if has_suntag else np.empty(shape=rna1.shape)
        )

        spots_cell = quantification.spots_quantif(
            rna1_clustered_spots_dataframe,
            rna2_clustered_spots_dataframe,
            suntag_clustered_spots_dataframe,
            cell_label= cell_label,
            nucleus_label=nucleus_label,
            fov_cell = fov_cell,
            rna1_signal = rna1,
            rna2_signal = rna2 if has_rna2 else np.empty(shape=rna1.shape),
            suntag_signal = suntag if has_suntag else np.empty(shape=rna1.shape)
        )

        fov_cell["acquisition_id"] = input_id
        spots_cell["acquisition_id"] = input_id
        max_cell_id = Cell["cell_id"].max() if len(Cell) > 0 else 0
        fov_cell["cell_id"] +=  max_cell_id #increatment cell_id to stay unique across all fovs
        spots_cell["cell_id"] += max_cell_id #increatment cell_id to stay unique across all fovs
        
        Cell = pd.concat(
            [Cell, fov_cell],
            axis = 0
        )

        Spots = pd.concat(
            [Spots, spots_cell],
            axis = 0
        )

        #Acquisition
        Acquisition.at[input_id, "rna1_threshold"] = rna1_threshold
        Acquisition.at[input_id, "colocalisation_distance"] = colocalisation_distance
        if has_rna2 : Acquisition.at[input_id, "rna2_threshold"] = rna2_threshold
        if has_suntag : Acquisition.at[input_id, "suntag_threshold"] = suntag_threshold

        Acquisition.at[input_id, "rna1_deconvolution_sucess"] = len(rna1_spots_postdecomp) > 0
        if has_rna2 : Acquisition.at[input_id, "rna2_deconvolution_sucess"] = len(rna2_spots_postdecomp) > 0
        if has_suntag : Acquisition.at[input_id, "suntag_deconvolution_sucess"] = len(suntag_spots_postdecomp) > 0

    print('\033[32m' + "End of acquisition group : saving results." + '\033[00m')
    output_path = path_out + "/result_tables/"
    os.makedirs(output_path, exist_ok= True)
    Acquisition.sort_values('acquisition_id').reset_index(drop=True).to_feather(output_path + "Acquisition.feather")
    Acquisition.sort_values('acquisition_id').reset_index(drop=True).to_excel(output_path + "Acquisition.xlsx")
    Cell.sort_values('cell_id').reset_index(drop=True).to_feather(output_path + "Cell.feather")
    Cell.sort_values('cell_id').reset_index(drop=True).to_excel(output_path + "Cell.xlsx")
    Spots.sort_values('spots_id').reset_index(drop=True).to_feather(output_path + "Spots.feather")

print('\033[32m' + "End of pipeline : integrity checks."+ '\033[00m')
integrity_log = integrity.full_integrity_check(Acquisition, Cell, Spots)
integrity.write_integrity_log(path_out, integrity_log)

print("Done.")
print('pipelinetime : ', time.process_time() - clock)
