"""
This submodule aims at realigning images using fluorescent bead markers.
"""

import numpy as np
import pandas as pd
import pbwrap.detection as detection
import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt
import warnings

class Match_Error(Exception) :
    """
    Exception class rosen when a beads matching test fails.
    """
    pass

def shift_array(arr : np.ndarray,*args) :
    indexer_new_array = []
    indexer_old_array = []
    for delta in args :
        if delta == 0 : 
            indexer_new_array.append(slice(None))
            indexer_old_array.append(slice(None))
        elif delta > 0 :
            indexer_new_array.append(slice(delta,None))
            indexer_old_array.append(slice(None,-delta))

        else :
            indexer_new_array.append(slice(None,delta))
            indexer_old_array.append(slice(-delta,None))

    if len(args) < arr.ndim :
        indexer_old_array.append(...)
        indexer_new_array.append(...)

    indexer_new_array = tuple(indexer_new_array)
    indexer_old_array = tuple(indexer_old_array)
    new_arr = np.zeros_like(arr)

    if len(args) > arr.ndim :
        raise ValueError("too many axis to shift; dim : {0}, shift : {1}".format(arr.ndim, args))
    else :
        new_arr[indexer_new_array] = arr[indexer_old_array]
        new_arr[indexer_new_array] = arr[indexer_old_array]

    return new_arr

def _drift_statistics_plots(distance_df: pd.DataFrame, path_output: str) :
    
    fig,axes = plt.subplots(1,2, sharex=True, sharey= True)
    
    distance_df = distance_df.loc[:,['channel','dz','dy','dx']].dropna(axis=1, how='all').dropna(axis=0, how='any')

    #Reference subplot (left)
    distance_df.loc[distance_df['channel'] == 'reference'].plot(
        kind='hist',
        ax=axes[0],
        bins= 10,
        alpha = 0.7,
        title= 'Reference beads'
    )

    #Drift subplot (right)
    distance_df.loc[distance_df['channel'] == 'drift'].plot(
        kind='hist',
        ax=axes[1],
        bins= 10,
        alpha = 0.7,
        title= 'Drift beads'
    )

    fig.savefig(path_output + "individual_axis_drift_distribution.png")
    plt.close(fig=fig)

def _safe_mode(ser : pd.Series) :

    res = ser.mode()

    if len(res) == 1 :
        return res
    else :
        return np.NaN

def _compute_drift_statistics(coordinates_df : pd.DataFrame, path_output= None) :
    
    aggregate_df = coordinates_df.loc[:,['channel','dz','dy','dx']]
    aggregate_df = aggregate_df.melt(
        id_vars= 'channel',
        value_vars= ['dz','dy','dx'],
        value_name= 'drift',
        var_name= 'dim'
    )
    aggregate_df = aggregate_df.groupby(['channel','dim']).agg(['mean','median','std', _safe_mode])
    if type(path_output) != type(None) : 
        _drift_statistics_plots(distance_df=coordinates_df, path_output=path_output)
        aggregate_df.rename(columns={'_safe_mode' : 'most_frequent'}).to_csv(path_output + "drift_statistics.csv")

    return aggregate_df


def _compute_axis_drift(coordinates_df : pd.DataFrame) :
    for ax in ['z','y','x'] :
        coordinates_df['d' + ax] = coordinates_df['matching_' + ax] - coordinates_df[ax]
    
    return coordinates_df

def _find_drift_value(matching_coordinates_df: pd.DataFrame, path_output = None) -> 'tuple[int]':
    
    matching_coordinates_df = _compute_axis_drift(matching_coordinates_df)
    statistics_df = _compute_drift_statistics(
        coordinates_df= matching_coordinates_df,
        path_output=path_output,
    )

    reference_drift = tuple(list(statistics_df.loc[('reference'), ('drift','_safe_mode')]))
    if any(np.isnan(reference_drift)) :
        reference_drift = tuple(list(statistics_df.loc[('reference'), ('drift','median')].astype(int)))

    drift_drift = tuple(list(statistics_df.loc[('drift'), ('drift','_safe_mode')]))
    if any(np.isnan(drift_drift)) :
        drift_drift = tuple(list(statistics_df.loc[('drift'), ('drift','median')].astype(int)))

    return reference_drift, drift_drift

def _apply_distance_threshold(
        distance_df: pd.DataFrame,
        distance_threshold : float
        ) :

    distance_df = distance_df.loc[distance_df['distance_to_closest_spot'] <= distance_threshold]

    return distance_df

def _build_coordinates_df(
        reference_beads,
        drifted_beads,
        reference_indices,
        drift_indices,
        distance_reference_from_drift,
        distance_drift_from_reference,
        dim
) :
    
    ## 1.Preparing data
    is_3D: bool = dim == 3
    distance_to_closest_spot = {
        'reference' : distance_reference_from_drift
        ,'drift' : distance_drift_from_reference
    }
    beads_number = {
        'reference' : len(reference_beads)
        ,'drift' : len(drifted_beads)
    }
    unpacked_coordinates = { #Coordinates of beads
        'reference': tuple(list(zip(*reference_beads)))
        ,'drift' : tuple(list(zip(*drifted_beads)))
        
    }
    matching_unpacked_coordinates = { #Coordinates of matching bead (meaning from the other channel)
        'reference' : [ax[unpacked_coordinates['reference']] for ax in drift_indices]
        ,'drift' : [ax[unpacked_coordinates['drift']] for ax in reference_indices]
    }

    ## 2.Bulding DataFrame
    distance_df = pd.concat([
        pd.DataFrame({
            'channel' : [channel] * beads_number[channel]
            ,'distance_to_closest_spot' : distance_to_closest_spot[channel]
            ,'z' : unpacked_coordinates[channel][0] if is_3D else [np.NaN] * beads_number[channel]
            ,'y' : unpacked_coordinates[channel][0 + int(is_3D)] #0 or 1 if True
            ,'x' : unpacked_coordinates[channel][1 + int(is_3D)] #1 or 2 if True
            ,'matching_z' : matching_unpacked_coordinates[channel][0] if is_3D else [np.NaN] * beads_number[channel]
            ,'matching_y' : matching_unpacked_coordinates[channel][0 + int(is_3D)] #0 or 1 if True
            ,'matching_x' : matching_unpacked_coordinates[channel][1 + int(is_3D)] #1 or 2 if True
            })
    for channel in ['reference', 'drift']], axis= 0)

    distance_df = distance_df.reset_index(drop=False)

    return distance_df

def _find_distance_threshold(
        distance_reference_from_drift,
        distance_drift_from_reference,
        output_path = None,
        bins = 50,
        talk=False,
) -> float :
    fig = plt.figure()
    ax = fig.gca()
    count_reference,reference_value,_ = ax.hist(distance_reference_from_drift, bins=bins, label="reference beads to drifted beads", alpha = 0.5, color='blue')
    count_drift, drift_value,_ = ax.hist(distance_drift_from_reference, bins=reference_value, label="drifted beads to reference beads", alpha = 0.5, color='red')

    ax.set_xlabel('Distance to closest bead (nm)')
    ax.set_ylabel('count')
    ax.legend()

    if type(output_path) != type(None) : fig.savefig(output_path + "distance_to_closest_bead_histogram.png")
    if talk : plt.show()
    plt.close(fig)

    reference_max_count = np.argmax(count_reference)
    drift_max_count = np.argmax(count_drift)
    reference_threshold = reference_value[reference_max_count + 1]
    drift_threshold = drift_value[drift_max_count + 1]

    distance_threshold: float = reference_threshold
    return distance_threshold


def _get_distance_from_map(
        spots : np.ndarray,
        distance_map : np.ndarray,
) -> np.ndarray :
    
    coordinates_index = tuple(list(zip(*spots)))
    distances = distance_map[coordinates_index]

    return distances

def _build_maps(
        spots : np.ndarray,
        voxel_size : tuple,
        shape : tuple
) :
    spot_signal = np.ones(shape=shape, dtype=bool)
    coordinates_index = tuple(list(zip(*spots)))
    spot_signal[coordinates_index] = 0
    distance_map, indices = distance_transform_edt(
        input=spot_signal,
        sampling=voxel_size,
        return_distances=True,
        return_indices=True,
    )

    return distance_map, indices


def _detect_beads(
    reference_bead_signal : np.ndarray,
    drifted_bead_signal : np.ndarray,
    voxel_size : 'tuple[float]',
    bead_size : 'tuple[float]',
    ndim : int,
    reference_threshold,
    reference_threshold_penalty,
    drift_threshold,
    drift_threshold_penalty,
    extended_result = False,

) :
    """
    For the reference channel we apply a higher threshold than automatic to ensure having only true spots. 
    For drifted we use a lower threshold so that all reference beads find their match, the false spots being filter out 
    when distance threshold is applied down in the process.
    """
    
    reference_beads,*reference_threshold = detection.detect_spots(
        images= reference_bead_signal,
        threshold= reference_threshold,
        threshold_penalty= reference_threshold_penalty,
        voxel_size= voxel_size,
        spot_radius= bead_size,
        return_threshold=extended_result,
        ndim=ndim
    )

    drifted_beads, *drift_threshold = detection.detect_spots(
        images= drifted_bead_signal,
        threshold= drift_threshold,
        threshold_penalty= drift_threshold_penalty,
        voxel_size= voxel_size,
        spot_radius= bead_size,
        return_threshold=extended_result,
        ndim=ndim
    )

    if len(reference_beads) < 5 : warnings.warn("Less than 5 beads were detected for reference beads channel")
    if len(drifted_beads) < 5 : warnings.warn("Less than 5 beads were detected for drifted beads channel")

    if extended_result :
        return reference_beads, drifted_beads, reference_threshold[0], drift_threshold[0]
    else :
        return reference_beads, drifted_beads 

def find_drift(
        reference_bead_signal : np.ndarray,
        drifted_bead_signal : np.ndarray,
        voxel_size : 'tuple[float]',
        bead_size : 'tuple[float]',
        reference_threshold = None,
        reference_threshold_penalty = 2,
        drift_threshold = None,
        drift_threshold_penalty = 1,
        plot_path = None,
        talk= False,
        extended_result=False,
) : 
    
    """
    Correction for spatial uniform translation using signal from fluorescent beads.
    Beads are detected using bigFish spot detection in both channels, they are then paired up using euclidian distance distribution between closest beads.
    Drift is then computed individually along all axis (2 or 3) as median drift between pairs.

    Parameters
    ----------
        reference_bead_signal : np.ndarray; ndim = {2,3},
            bead signal from reference channel.
        drifted_bead_signal : np.ndarray; ndim = reference_bead_signal.ndim,
            bead signal from channel to correct drift on.
        voxel_size : 'tuple[float]'; len = signal.ndim
            size of a pixel in nanometer (Z,Y,X)
        bead_size : 'tuple[float]',
            size of a bead in nanometer (Z,Y,X)
        reference_threshold = None,
            Optional threshold for bead detection
        reference_threshold_penalty = 10,
            Penalty to apply for automatic threshold on bead detection (recomanded > 1)
        drift_threshold = None,
            Optional threshold for bead detection
        drift_threshold_penalty = 1,
            Penalty to apply for automatic threshold on bead detection (recomanded ~ 1)
        plot_path = None,
            if passed, statics plots from distance and drift distributions.

    Returns
    -------
        corrected_signal : np.ndarray
            Signal from `drifted_bead_signal` with drift correction applied.
    """

    #Type Error
    if not isinstance(reference_bead_signal, np.ndarray) : raise TypeError("Expected np.ndarray got {0}".format(type(reference_bead_signal)))
    if not isinstance(drifted_bead_signal, np.ndarray) : raise TypeError("Expected np.ndarray got {0}".format(type(drifted_bead_signal)))
    if not isinstance(voxel_size, (tuple, list)) : raise TypeError("Expected tuple got {0}".format(type(drifted_bead_signal)))
    if not isinstance(bead_size, (tuple, list)) : raise TypeError("Expected tuple got {0}".format(type(bead_size)))

    #Value Error
    if reference_bead_signal.ndim not in [2,3] : raise ValueError("Incorrect dimension for reference bead signal, got {0}".format(reference_bead_signal.ndim))
    if drifted_bead_signal.ndim != reference_bead_signal.ndim : raise ValueError("Incorrect dimension for reference bead signal, got {0} expected {1}".format(drifted_bead_signal.ndim, reference_bead_signal.ndim))
    if len(voxel_size) != reference_bead_signal.ndim : raise ValueError("Incorrect dimension number in voxel size, got {0} expected {1}".format(len(voxel_size), reference_bead_signal.ndim))
    if len(bead_size) != reference_bead_signal.ndim : raise ValueError("Incorrect dimension number in bead size, got {0} expected {1}".format(len(bead_size), reference_bead_signal.ndim))
    
    dim = len(voxel_size)

    reference_beads, drifted_beads, *thresholds = _detect_beads(
        reference_bead_signal,
        drifted_bead_signal,
        voxel_size,
        bead_size,
        dim,
        reference_threshold=reference_threshold,
        reference_threshold_penalty=reference_threshold_penalty,
        drift_threshold=drift_threshold,
        drift_threshold_penalty=drift_threshold_penalty,
        extended_result = extended_result,
        )
    
    if talk :
        print(
            "reference beads number : ", len(reference_beads),
            "\ndrift beads nuber : ", len(drifted_beads)
        )
    
    #Finding threshold to consider beads are matching
    shape = np.max([reference_bead_signal.shape, drifted_bead_signal.shape], axis=0)
    reference_distance_map, reference_indices = _build_maps(spots= reference_beads, voxel_size=voxel_size, shape=shape)
    drifted_distance_map, drifted_indices = _build_maps(spots=drifted_beads, voxel_size=voxel_size, shape=shape)
    distance_reference_from_drift = _get_distance_from_map(spots= reference_beads, distance_map= drifted_distance_map)
    distance_drift_from_reference = _get_distance_from_map(spots= drifted_beads, distance_map= reference_distance_map)
    distance_threshold = _find_distance_threshold(
        distance_reference_from_drift, 
        distance_drift_from_reference,
        output_path= plot_path,
        talk=talk,
        )

    if talk : print ("Distance threshold : ", distance_threshold)

    coordinates_df = _build_coordinates_df(
        reference_beads,
        drifted_beads,
        reference_indices,
        drifted_indices,
        distance_reference_from_drift,
        distance_drift_from_reference,
        dim=dim
    )

    coordinates_df = _apply_distance_threshold(
        coordinates_df,
        distance_threshold
    )

    reference_drift, drift = _find_drift_value(
        matching_coordinates_df=coordinates_df,
        path_output=plot_path,
        )
    
    if extended_result :
        drift_z, drift_y, drift_x = drift
        result = {
            'drift_z' : [drift_z],
            'drift_y' : [drift_y],
            'drift_x' : [drift_x],
            'ref_bead_threshold' : [thresholds[0]],
            'drift_bead_threshold' : [thresholds[1]],
            'ref_bead_number' : [len(reference_beads)],
            'drift_bead_number' : [len(drifted_beads)],
            'found_symetric_drift' : [all(np.array(reference_drift) == -np.array(drift))]
        }
        return result
    else :
        return drift