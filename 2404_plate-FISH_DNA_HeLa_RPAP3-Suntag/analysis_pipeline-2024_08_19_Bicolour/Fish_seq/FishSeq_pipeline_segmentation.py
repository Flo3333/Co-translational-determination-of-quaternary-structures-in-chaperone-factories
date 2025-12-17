import os
import dataprocessing.fishseq_preprocess_folder as prepro
import pbwrap.segmentation as segm
import numpy as np
import bigfish.plot as plot
from pbwrap.utils import open_image
from tqdm import tqdm

"""
This script aims at performing segmentation and savings results as .npy format to be used in FishSeq_pipeline_quantification.py
"""


#### USER PARAMETERS

RUN_PATH = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2024-03-05 - R2TP_Pool-1_Run1/"

SEGMENTATION_PATH = {
    'nucleus' : RUN_PATH + "/DAPI_Z-stacks/",
    'cytoplasm' : RUN_PATH + "/Cy3_Z-stacks/",
}

MODEL_DICT = {
    'nucleus' : 'nuclei',
    'cytoplasm' : 'cyto3'
}

OBJECT_SIZE_DICT = {
    'nucleus' : 140,
    'cytoplasm' : 200
}

SAVE_PATH = None

PLOT_VISUALS = True


############

# MAIN SCRIPT

############


#Reading input folder.
file_dict = prepro.assert_run_folder_integrity(run_path=RUN_PATH)
location_list = list(file_dict.keys())
print("{0} locations found.".format(len(location_list)))

#Setting output folder.
if type(SAVE_PATH) == type(None) :
    SAVE_PATH = RUN_PATH + "/segmentation_files/"
elif not os.path.isdir(SAVE_PATH) :
    raise ValueError("Incorrect path for SAVE_PATH.")
os.makedirs(SAVE_PATH, exist_ok=True)

print("Running segmentation")
for location in tqdm(location_list, total= len(location_list)) :

    #Loading file
    filename = file_dict[location][0]
    nucleus_image = open_image(SEGMENTATION_PATH['nucleus'] + "{0}/{1}".format(location,filename))
    nucleus_image = np.mean(nucleus_image, axis=0)
    cytoplasm_image = open_image(SEGMENTATION_PATH['cytoplasm'] + "{0}/{1}".format(location,filename))
    cytoplasm_image = np.mean(cytoplasm_image, axis=(0,1))

    #Segmentation
    nucleus_label = segm.Nucleus_segmentation(
        dapi=nucleus_image,
        diameter=OBJECT_SIZE_DICT['nucleus'],
        model_type= MODEL_DICT['nucleus'],
        use_gpu= True
    )
    cytoplasm_label = segm.Cytoplasm_segmentation(
        cy3=cytoplasm_image,
        dapi=nucleus_image,
        diameter=OBJECT_SIZE_DICT['cytoplasm'],
        model_type=MODEL_DICT['cytoplasm'],
        use_gpu=True
    )

    #Saving labels
    np.savez(
        file= SAVE_PATH + "{0}_segmentation".format(location),
        nucleus= nucleus_label,
        cytoplasm= cytoplasm_label,
    )

    if PLOT_VISUALS : 
        plot.plot_segmentation_boundary(
            image=cytoplasm_image,
            cell_label=cytoplasm_label,
            nuc_label=nucleus_label,
            boundary_size=3,
            contrast=True,
            path_output=SAVE_PATH + "{0}_segmentation_cyto_view.png".format(location),
            show=False
        )
        plot.plot_segmentation_boundary(
            image=nucleus_image,
            cell_label=cytoplasm_label,
            nuc_label=nucleus_label,
            boundary_size=3,
            contrast=True,
            path_output=SAVE_PATH + "{0}_segmentation_nuc_view.png".format(location),
            show=False
        )
