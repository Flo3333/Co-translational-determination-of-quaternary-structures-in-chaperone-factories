import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CustomPandasFramework.centrosome.analysis as analysis

from CustomPandasFramework.update import from_IntegratedSignal_spike_compute_CellularCycleGroup

def _plot_IntegratedSignal_hist(ax: plt.Axes, df: pd.DataFrame) :
    pass

def _gaussian_fit() :
    pass

def classification_analysis(input_path: str, frameon=True) :
    """
    Idea : create histogram of Integrated signal per well.
    + Fit 2 gaussian (G1 + G2) and look at quality of fit ? G1 center --> Hist maximum; G2 center = 2*G1 center
    """
    
    df = analysis.prepare_data_frame(
        input_path,
        data_columns= ['nucleus_area_nm', 'nucleus_mean_mean_signal'],
        index_name= []
        )

    df = from_IntegratedSignal_spike_compute_CellularCycleGroup(df, surface_column=' nucleus_area_nm')

    figsize = (10,10)
    fig = plt.figure(figsize=figsize, frameon=frameon)
    count, bins = _plot_IntegratedSignal_hist()
    score = _gaussian_fit(df)
    plt.show()

def _g1g2_scatter_plot(ax, data) :
    pass

def cellular_cycle_centrosome_analysis(input_path, measures, frameon=True) :

    df = analysis.prepare_data_frame(
        input_path,
        data_columns= ['nucleus_area_nm', 'nucleus_mean_mean_signal'] + measures,
        index_name= []
        )

    df = from_IntegratedSignal_spike_compute_CellularCycleGroup(df, surface_column=' nucleus_area_nm')

    grouper = df.groupby(['rna', 'treatment', 'cellular_cycle'])

    for measure in measures :

        data = grouper[measure].apply(list)

        figsize = (10,10)
        fig = plt.figure(figsize=figsize, frameon=frameon)
        ax = fig.gca()
        ax = _g1g2_scatter_plot(ax, data)