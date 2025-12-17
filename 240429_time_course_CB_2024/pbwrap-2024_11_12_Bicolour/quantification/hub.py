import numpy as np
import pandas as pd
import bigfish.classification as classification
import bigfish.plot as plot
import warnings
from skimage.measure import regionprops_table
from .measures import compute_mask_area, nucleus_signal_metrics
from ._errors import QuantificationError
from ..detection.centrosome import detect_centrosome

def compute_Cell(cell: dict, voxel_size, dapi_stack, acquisition_id, pipeline_name: str ,**pipeline_parameters) :

    if 'centrosome' in pipeline_name :
        new_cell = _centrosome_cell_quant(cell, voxel_size, dapi_stack, acquisition_id)
    
    else : raise ValueError("pipeline name not recognised.")

    return new_cell


def _main_cell_quant(cell, voxel_size, dapi_stack=None, compute_centrosome= False, centrosome_coords= None) :
    """
    Basic function for cell quantification using bigfish AND COMPUTING DAPI MEASUREMENTS if dapi_stack != None.

    Parameters
    ----------
        cell : dict  
            computed from `bigfish.multistack.extract_cell`
        voxel_size : tuple(z,y,x)
        dapi_stack : np.ndarray
            3D uncropped dapi stack.

    Returns
    -------
        Cell : pd.DataFrame  
    """
    
    #Extracting bigfish cell information
    if isinstance(voxel_size, (list, np.ndarray)) : voxel_size = tuple(voxel_size)
    voxel_size_yx = float(voxel_size[1])
    cell_mask: np.ndarray = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"] 
    rna_coord = cell["rna_coord"]
    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]
    label = cell["cell_id"] # is the label of this cell in cell_label

    if "transcription_site" in cell.keys() :
        ts_coord = cell["transcription_site"]
    else : ts_coord = np.empty(shape=(0,0), dtype= np.int64)

    if "foci" in cell.keys() :
        foci_coord = cell["foci"]
    else : foci_coord = np.empty(shape=(0,0), dtype= np.int64)

    if not isinstance(centrosome_coords, (np.ndarray, list)) and compute_centrosome:
        print("Warning : compute_centrosome is set to True but centrosome_coords parameter is not np.ndarray nor list. \nCompute centrosome has been set to False; centrosome_coords type : {0}".format(type(centrosome_coords)))
        compute_centrosome = False
    if type(centrosome_coords) == list : centrosome_coords = np.array(centrosome_coords, dtype=int)
    
    warnings.filterwarnings("ignore")

    features, features_names = classification.compute_features(
        cell_mask= cell_mask, 
        nuc_mask= nuc_mask, 
        ndim= 3, 
        rna_coord= rna_coord, 
        smfish= smfish, 
        foci_coord= foci_coord, 
        voxel_size_yx= voxel_size_yx, 
        centrosome_coord= centrosome_coords,
        compute_centrosome=compute_centrosome,
        compute_distance=True,
        compute_intranuclear=True,
        compute_protrusion=True,
        compute_dispersion=True,
        compute_topography=True,
        compute_foci=True,
        compute_area=True,
        return_names=True
    )
    warnings.filterwarnings("default")                                                      

    #Custom features
    cell_props_table = regionprops_table(cell_mask.astype(int), properties= ["centroid"])
    cell_coordinates = (float(cell_props_table["centroid-0"] + min_y), float(cell_props_table["centroid-1"] + min_x))
    del cell_props_table
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
 
    #Adding custom features to DataFrame
    features = list(features)
    features.extend([cell_coordinates, label, cell["bbox"], nucleus_area_px, nucleus_area_nm,])
     
    features_names += [ "cell_coordinates", "label", "bbox", "nucleus_area_px", "nucleus_area_nm"]
    
    #signal features
    if type(dapi_stack) != type(None) :
        nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi_stack, projtype= 'mip')
        nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi_stack, projtype= 'mean')
        features.extend([nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                         nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"]])
     
        features_names += ["nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                           "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]



    data = features
    header = features_names

    new_cell = pd.DataFrame(data= [data], columns= header)
    return new_cell


def _get_minimal_columns_names(nucleus_quant=False) :
    names = [
            'index_mean_distance_cell',
            'index_median_distance_cell', 'index_mean_distance_nuc',
            'index_median_distance_nuc', 'proportion_rna_in_nuc', 'nb_rna_out_nuc',
            'nb_rna_in_nuc', 'index_rna_protrusion', 'proportion_rna_protrusion',
            'protrusion_area', 'index_polarization', 'index_dispersion',
            'index_peripheral_distribution', 'index_rna_nuc_edge',
            'proportion_rna_nuc_edge', 'index_rna_nuc_radius_500_1000',
            'proportion_rna_nuc_radius_500_1000', 'index_rna_nuc_radius_1000_1500',
            'proportion_rna_nuc_radius_1000_1500', 'index_rna_nuc_radius_1500_2000',
            'proportion_rna_nuc_radius_1500_2000', 'index_rna_nuc_radius_2000_2500',
            'proportion_rna_nuc_radius_2000_2500', 'index_rna_nuc_radius_2500_3000',
            'proportion_rna_nuc_radius_2500_3000', 'index_rna_cell_radius_0_500',
            'proportion_rna_cell_radius_0_500', 'index_rna_cell_radius_500_1000',
            'proportion_rna_cell_radius_500_1000',
            'index_rna_cell_radius_1000_1500',
            'proportion_rna_cell_radius_1000_1500',
            'index_rna_cell_radius_1500_2000',
            'proportion_rna_cell_radius_1500_2000',
            'index_rna_cell_radius_2000_2500',
            'proportion_rna_cell_radius_2000_2500',
            'index_rna_cell_radius_2500_3000',
            'proportion_rna_cell_radius_2500_3000', 'proportion_rna_in_foci',
            'proportion_nuc_area', 'cell_area', 'nuc_area', 'cell_area_out_nuc',
            'cell_coordinates', 'label', 'bbox', 'nucleus_area_px',
            'nucleus_area_nm']
    
    if nucleus_quant : nucleus_names = ["nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                           "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]

    return names
        
    
def _centrosome_cell_quant(cell, voxel_size, dapi_stack, acquisition_id) :
    # Note : it appears bigfish centrosome features only work with 2D coords, while it accepts 3D coords only 2D measures are computed and worse it fails to remove the z dimension properly (as of today's version)
    # --> we have to remove the z coords
    # TODO : open pull request on github
    
    if 'centrosome_spots' in cell.keys() :
        centrosome_coords = detect_centrosome(cell=cell)
        centrosome_number = len(centrosome_coords)
        if centrosome_number == 0 : raise QuantificationError("No centrosome found")
    else :
        centrosome_coords = np.empty(shape=(0,0), dtype=np.int32)
        centrosome_number = 0

    centrosome_coords_2d = centrosome_coords[:,1:]
    min_y,min_x,max_y,max_x = cell['bbox']
    min_z = 0

    cell_res = _main_cell_quant(cell=cell, voxel_size=voxel_size, dapi_stack=dapi_stack, compute_centrosome= True, centrosome_coords=centrosome_coords_2d)
    cell_res.at[0, "centrosome_number"] = centrosome_number

    if 'centrosome_spots' in cell.keys() :
        global_centrosome_coords = []
        local_centrosome_coords = []
        for local_coords in centrosome_coords :
            global_coordinates = []
            for offset, position in zip(local_coords, [min_z, min_y, min_x]) :
                global_coordinates.append(position + offset)
            global_centrosome_coords.append(tuple(global_coordinates))
            local_centrosome_coords.append(tuple(local_coords))
        global_centrosome_coords = [tuple(global_centrosome_coords)]
        local_centrosome_coords = [tuple(local_centrosome_coords)]
    else :
        global_centrosome_coords = np.NaN
        local_centrosome_coords = np.NaN

        
    cell_res["centrosome_coords"] = global_centrosome_coords
    cell_res["centrosome_coords_local"] = local_centrosome_coords
    cell_res['acquisition_id'] = acquisition_id
    clusters_quant = _clusters_quant(cell)
    cell_res = pd.concat([cell_res, clusters_quant], axis= 1).reset_index(drop=True)

    return cell_res

def _clusters_quant(cell: dict, clusters_coords_key= 'clusters_coords', clustered_spots_key = 'clustered_spots', unclustered_spots_key = 'unclustered_spots') :

    """
    Keys : cluster_number, nucleus_cluster_number, clustered_spots_number, unclustered_spots_number, clustered_spots_fraction
    """

    ymin, xmin, ymax, xmax = cell.get('bbox')
    clustered_spots = cell.setdefault(clustered_spots_key, np.empty(shape=(0,0,0), dtype= int))
    unclustered_spots = cell.setdefault(unclustered_spots_key, np.empty(shape=(0,0,0), dtype= int))
    clusters_coords = cell.setdefault(clusters_coords_key, np.empty(shape=(0,0,0), dtype= int))

    if len(clusters_coords) == 0 :
        Z = Y = X = []
    else :
        Z,Y,X = zip(*clusters_coords)

    nucleus_mask = cell.get("nuc_mask")

    if type(clustered_spots) == type(None) : raise KeyError("clustered_spots array not found in extracted cells")
    if type(unclustered_spots) == type(None) : raise KeyError("unclustered_spots array not found in extracted cells")
    if type(clusters_coords) == type(None) : raise KeyError("clusters_coords array not found in extracted cells")

    clustered_spots_number = len(clustered_spots)
    unclustered_spots_number = len(unclustered_spots)

    res = pd.DataFrame({
        "cluster_number" : [len(clusters_coords)],
        "nucleus_cluster_number" : [sum(nucleus_mask[Y,X])],
        "clustered_spots_number" : [clustered_spots_number],
        "unclustered_spots_number" : [unclustered_spots_number],
        "clustered_spots_fraction" : [clustered_spots_number/(clustered_spots_number + unclustered_spots_number)] if clustered_spots_number + unclustered_spots_number != 0 else [np.NaN]
    })

    return res