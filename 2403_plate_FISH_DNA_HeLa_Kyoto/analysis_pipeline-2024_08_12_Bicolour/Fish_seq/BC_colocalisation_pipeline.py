import os
import time, signal
import CustomPandasFramework.FishBicolor.framework as framework
import CustomPandasFramework.FishBicolor.quantification as quantification
import CustomPandasFramework.FishBicolor.DataFrames as df
import bigfish.stack as stack
import pbwrap.segmentation as segmentation
import pbwrap.plot as plt
import pbwrap.detection as detection
import pbwrap.preprocessing.preprocessing as preprocessing
import CustomPandasFramework.log as logs
from pbwrap.integrity import detectiontimeout_handler
from CustomPandasFramework.computer_interface import get_datetime

path_in = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_input/"
path_out = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/"

#### PARAMETERS ####
#channel from 0 to 3
suntag_channel = 1
rna1_channel = 2
rna2_channel = 3
cell_diameter = 150 #pixels
spot_rna = 0.75
spot_suntag = 0.75
rna_spot_radius = (int(round(350*spot_rna)),int(round(150*spot_rna)),int(round(150*spot_rna)))
suntag_spot_radius = (int(round(350*spot_suntag)),int(round(150*spot_suntag)),int(round(150*spot_suntag)))
voxel_size = (300, 103, 103) #nm
threshold_penalty_rna1 = 2
threshold_penalty_rna2 = 1.5
threshold_penalty_suntag = 2
alpha = 0.4
alpha_suntag = 0.8
beta_suntag = 4
beta = 1.2
_crop_pixel = 10
colocalisation_distance = 310 # ~ 3 pixels
rna_cluster_radius = 750
min_spot_number_per_cluster = 3

suntag_artifact_radius = 103 * 10

do_detection_visuals = False # Set True if you want to check the dense region decomposition or the cluster detection else colocalisation_visuals are enough.
do_colocalisation_visuals = True

#### #####
#Pipeline
clock = time.process_time()
Input = framework.create_Input(path_in)
print(Input)

analysis_group = framework.make_analysis_group(Input=Input)
date = get_datetime()
path_out += date + '/'
os.makedirs(path_out, exist_ok= True)
signal.signal(signal.SIGALRM, detectiontimeout_handler) #Initiating timeout handling

# Log
parameters_log = logs.parameter_log("parameters.txt", path= path_out)
parameters_log.add_parameters(cell_diameter, spot_rna, spot_suntag, rna_spot_radius, suntag_spot_radius, suntag_spot_radius, voxel_size, threshold_penalty_rna1, threshold_penalty_rna2, threshold_penalty_suntag, alpha, alpha_suntag, beta, beta_suntag,_crop_pixel, colocalisation_distance, colocalisation_distance, rna_cluster_radius, min_spot_number_per_cluster)
parameters_log.write()
Acquisition = framework.create_Acquisiton(Input)
Colocalisation = df.create_EmptyColocalisation()
ClusteredSpots = df.create_EmptyClusteredSpots()

for group_index in analysis_group.index :
    group = analysis_group.loc[group_index]
    rna1, rna2 = group_index
    print("\nNext group : ({0};{1})".format(rna1, rna2))
    output_path = path_out + "/{0}_{1}/".format(rna1,rna2)
    os.makedirs(output_path, exist_ok= True)

    print("detecting spots for rna1...")
    rna1_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, channels= rna1_channel, _crop_pixel= _crop_pixel)
    rna1_threshold = detection.compute_auto_threshold(rna1_gen, voxel_size=voxel_size, spot_radius=rna_spot_radius, im_number=8) * threshold_penalty_rna1
    rna1_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, channels= rna1_channel, _crop_pixel= _crop_pixel)
    rna1_spots_list = detection.detect_spots(rna1_gen, threshold=rna1_threshold, ndim= 3, voxel_size= voxel_size, spot_radius= rna_spot_radius)
    rna1_spots_gen = iter(rna1_spots_list)
    del rna1_gen, rna1_spots_list

    print("detecting spots for rna2...")
    rna2_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, channels= rna2_channel, _crop_pixel= _crop_pixel)
    rna2_threshold = detection.compute_auto_threshold(rna2_gen, voxel_size=voxel_size, spot_radius=rna_spot_radius, im_number=8) * threshold_penalty_rna2
    rna2_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, channels= rna2_channel, _crop_pixel= _crop_pixel)
    rna2_spots_list = detection.detect_spots(rna2_gen, threshold= rna2_threshold, ndim= 3, voxel_size= voxel_size, spot_radius= rna_spot_radius)
    rna2_spots_gen = iter(rna2_spots_list)
    del rna2_gen, rna2_spots_list

    print("detecting spots for suntag...")
    suntag_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, channels= suntag_channel, _crop_pixel= _crop_pixel)
    suntag_threshold = detection.compute_auto_threshold(suntag_gen, voxel_size=voxel_size, spot_radius=rna_spot_radius, im_number=8) * threshold_penalty_suntag
    suntag_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, channels= suntag_channel, _crop_pixel= _crop_pixel)
    suntag_spots_list = detection.detect_spots(suntag_gen, threshold=suntag_threshold, ndim= 3, voxel_size= voxel_size, spot_radius= suntag_spot_radius)
    suntag_spots_gen = iter(suntag_spots_list)
    del suntag_gen, suntag_spots_list

    im_gen = framework.get_image_as_gen(Input=Input, path_in= path_in, index_list= group, _crop_pixel= _crop_pixel)

    for index in group :
        filename = Input.at[index, 'filename'].replace('czi','')
        bool_rna1_decomp = True
        bool_rna2_decomp = True
        bool_suntag_decomp = True

        print("Computing fov {0}".format(filename))
        im = next(im_gen)
        shape = im.shape[1:]
        im[suntag_channel] = preprocessing.remove_mean_gaussian_background(im[suntag_channel], sigma= 1, voxel_size= voxel_size)
        rna1_spots = next(rna1_spots_gen)
        rna2_spots = next(rna2_spots_gen)
        suntag_spots = next(suntag_spots_gen)
        rna1_max = stack.maximum_projection(im[2])
        rna2_max = stack.maximum_projection(im[3])
        suntag_max = stack.maximum_projection(im[1])

        print(" Spots clusters deconvolution...")
        rna1_spots_postdecomp = detection.cluster_deconvolution(im[2], rna1_spots, rna_spot_radius, voxel_size=voxel_size, alpha= alpha, beta= beta, timer= 60)
        rna2_spots_postdecomp = detection.cluster_deconvolution(im[3], rna2_spots, rna_spot_radius, voxel_size=voxel_size, alpha= alpha, beta= beta, timer= 60)
        suntag_spots_postdecomp = detection.cluster_deconvolution(im[1], suntag_spots, rna_spot_radius, voxel_size=voxel_size, alpha= alpha_suntag, beta= beta_suntag, timer= 60)
        suntag_spots_without_artifacts = detection.remove_artifact(suntag_spots_postdecomp, artifact_radius=suntag_artifact_radius, voxel_size=voxel_size, spot_density= 2)
        suntag_spots_without_artifacts = detection.remove_artifact(suntag_spots_without_artifacts, artifact_radius=25*106, voxel_size=voxel_size, spot_density= 2)
        print(" {0} suntag artifact spots were removed.".format(len(suntag_spots_postdecomp) - len(suntag_spots_without_artifacts)))

        print(" Clusters computation...")
        rna1_clustered_spots_dataframe= detection.cluster_detection(rna1_spots_postdecomp, voxel_size= voxel_size, nb_min_spots=min_spot_number_per_cluster,radius=rna_cluster_radius, keys_to_compute=['clustered_spots_dataframe'])['clustered_spots_dataframe']
        rna2_clustered_spots_dataframe= detection.cluster_detection(rna2_spots_postdecomp, voxel_size= voxel_size, nb_min_spots=min_spot_number_per_cluster,radius=rna_cluster_radius, keys_to_compute=['clustered_spots_dataframe'])['clustered_spots_dataframe']

        print(" Plotting check visuals...")
        rna1_cluster_spots = detection.get_centroids_list(rna1_clustered_spots_dataframe.query("not cluster_id.isna()"))
        rna1_uncluster_spots = detection.get_centroids_list(rna1_clustered_spots_dataframe.query("cluster_id.isna()"))
        rna2_cluster_spots = detection.get_centroids_list(rna2_clustered_spots_dataframe.query("not cluster_id.isna()"))
        rna2_uncluster_spots = detection.get_centroids_list(rna2_clustered_spots_dataframe.query("cluster_id.isna()"))
        if do_detection_visuals :
            plt.output_spot_tiffvisual(rna1_max, [rna1_spots, rna1_spots_postdecomp, rna1_cluster_spots, rna1_uncluster_spots], output_path + "{0}_rna1_detection.tif".format(filename), dot_size= 2)
            plt.output_spot_tiffvisual(rna2_max, [rna2_spots, rna2_spots_postdecomp , rna2_cluster_spots, rna2_uncluster_spots], output_path + "{0}_rna2_detection.tif".format(filename), dot_size= 2)
            plt.output_spot_tiffvisual(suntag_max, [suntag_spots, suntag_spots_postdecomp, suntag_spots_without_artifacts], output_path + "{0}_suntag_detection.tif".format(filename), dot_size= 2)
            plt.output_spot_tiffvisual(rna1_max, [rna1_spots_postdecomp, rna2_spots_postdecomp, suntag_spots_without_artifacts], output_path + "{0}_rna1_coloc.tif".format(filename), dot_size= 3)
        if do_colocalisation_visuals and Acquisition.at[index, "fov_number"] == '01' :
            plt.colocalisation_plot(im[rna1_channel], shape, voxel_size=voxel_size, colocalisation_distance= colocalisation_distance, path_output= output_path + "{0}_cy3_coloc.tif".format(filename), spot_list1=rna1_spots_postdecomp, spot_list2= rna2_spots_postdecomp, spot_list3= suntag_spots_without_artifacts)
            plt.colocalisation_plot(im[rna2_channel], shape, voxel_size=voxel_size, colocalisation_distance= colocalisation_distance, path_output= output_path + "{0}_cy5_coloc.tif".format(filename), spot_list1=rna2_spots_postdecomp, spot_list2= rna1_spots_postdecomp, spot_list3= suntag_spots_without_artifacts)
            plt.output_spot_tiffvisual(suntag_max, [suntag_spots, suntag_spots_postdecomp, suntag_spots_without_artifacts], output_path + "{0}_suntag_detection.tif".format(filename), dot_size= 2)

        print(" Quantification...")
        Acquisition.at[index, "rna1_deconvolution_sucess"] = len(rna1_spots_postdecomp) != 0
        Acquisition.at[index, "rna2_deconvolution_sucess"] = len(rna2_spots_postdecomp) != 0
        Acquisition.at[index, "suntag_deconvolution_sucess"] = len(suntag_spots_postdecomp) != 0
        Acquisition.at[index, "colocalisation_distance"] = colocalisation_distance
        mask = segmentation.gaussian_threshold_segmentation(im[rna1_channel], sigma=1, percentile= 75, voxel_size=voxel_size)
        Colocalisation = quantification.compute_Colocalisation(Colocalisation, index, rna1_spots_postdecomp, rna2_spots_postdecomp, suntag_spots_without_artifacts, rna1_clusterd_spots= rna1_clustered_spots_dataframe, rna2_clusterd_spots= rna2_clustered_spots_dataframe, image_shape= shape, colocalisation_distance= colocalisation_distance, voxel_size=voxel_size)
        ClusteredSpots = quantification.compute_ClusteredSpots(ClusteredSpots, index, rna_clustered_spots_dataframe= rna1_clustered_spots_dataframe, sun_tag_spots= list(suntag_spots_without_artifacts), radius_nm=colocalisation_distance, image_shape= shape, voxel_size=voxel_size, rna_name= rna1, mask=mask)
        ClusteredSpots = quantification.compute_ClusteredSpots(ClusteredSpots, index, rna_clustered_spots_dataframe= rna2_clustered_spots_dataframe, sun_tag_spots= list(suntag_spots_without_artifacts), radius_nm=colocalisation_distance, image_shape= shape, voxel_size=voxel_size, rna_name = rna2, mask=mask)

print("End of pipeline : saving results.")
output_path = path_out + "/result_tables/"
os.makedirs(output_path, exist_ok= True)
Acquisition.sort_values('id').reset_index(drop=True).to_feather(output_path + "Acquisition.feather")
Colocalisation.sort_values('id').reset_index(drop=True).to_feather(output_path + "Colocalisation.feather")
ClusteredSpots.sort_values('id').reset_index(drop=True).to_feather(output_path + "ClusteredSpots.feather")
Acquisition.sort_values('id').reset_index(drop=True).to_excel(output_path + "Acquisition.xlsx")
Colocalisation.sort_values('id').reset_index(drop=True).to_excel(output_path + "Colocalisation.xlsx")
ClusteredSpots.sort_values('id').reset_index(drop=True).to_excel(output_path + "ClusteredSpots.xlsx")

print("Done.")
print('pipelinetime : ', time.process_time() - clock)