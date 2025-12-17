import warnings
import numpy as np
import pandas as pd

from .utils import hist_maximum


def from_IntegratedSignal_spike_compute_CellularCycleGroup(Cell : pd.DataFrame, column_name= 'cellular_cycle', surface_column='nucleus area (nm^2)', dapi_signal_column='nucleus_mean_mean_signal', shift = -0.12e12, auto_shift= True, multiplicator = 2, bins= 100,) :

    """
    Add new column to Cell DF : column_name (default = 'cellular_cycle')
        Values can be either one of ['g1', 'g2' or np.nan]

    Also add or re-compute 'IntegratedSignal' column.    
    """
    if column_name in Cell.index : 
        warnings.warn("{0} column already in dataframe, returning dataframe unchanged.".format(column_name))
        return Cell

    Cell[column_name] = np.NaN
    Cell.loc[:,["IntegratedSignal"]] = Cell[dapi_signal_column] * Cell[surface_column]
    distribution = Cell["IntegratedSignal"]
    median_value = np.median(Cell['IntegratedSignal'])
    maximum_value = hist_maximum(np.histogram(distribution, bins=bins))
    box_size = abs(median_value - maximum_value)

    if auto_shift : shift = (-np.percentile(distribution, 0.5)/2)

    g2_spike = (maximum_value + shift)*multiplicator


    indexg1 = Cell[Cell["IntegratedSignal"] < maximum_value + box_size].index
    indexg2 = Cell[Cell["IntegratedSignal"] > g2_spike].index

    Cell.loc[indexg1,column_name] = 'g1'
    Cell.loc[indexg2,column_name] = 'g2'


    return Cell