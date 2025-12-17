"""
Analysis script for hela kyoto
Pipeline : Fish_seq / BC_clusterkinetics_pipeline

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pbwrap.plot as plot
import pbwrap.quantification.statistical_test as stest
import CustomPandasFramework.pvalues as pvalues
from CustomPandasFramework.Fish_seq.plot import plot_coloc, add_pvalue, save_fig_dict
from itertools import chain


FIGSIZE = (15,15)
BINS = 100

DO_STEP2 = True #< Number of cell quantified
DO_STEP3 = False # fraction rna in foci
DO_STEP4 = False # proportion spots colocalisation
EXPECTED_TREATMENT_NUMBER = 2

PATH_IN= "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/merge_folder/"
PATH_OUT= "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/analysis/"

Z_SCORE = 2.575829 # 99% confidence interval

"""
STEP1 : Open data & merge
"""
Acquisition = pd.read_feather(PATH_IN + '/Acquisition.feather')
Cell = pd.read_feather(PATH_IN + '/Cell.feather')
# Spots = pd.read_feather(PATH_IN + '/Spots.feather')

for table in [Acquisition, Cell] :
    table['acquisition_id'] = table['acquisition_id'].apply(tuple)


if 'rna2' in Cell.columns : Cell = Cell.drop(columns='rna2')

Cell = pd.merge(
    Acquisition.loc[:,['acquisition_id', 'rna1', 'rna2']],
    Cell,
    on= 'acquisition_id',
    how='right'
)


"""
STEP 2 : Number of cell quantified.
"""

if DO_STEP2 :
    data = Cell.groupby(["rna1","rna2", "treatment"])['cell_id'].count().rename('cell number')
    data.to_excel("/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/analysis/quantified_cell_number.xlsx")

"""

Step 3 :  proportion Cy5 focci (Puro vs non puro)

"""

def add_significance(ax: plt.Axes, pvalues: pd.DataFrame, significance = 0.01) :
    xtickslabels = ax.get_xticklabels()

    if len(xtickslabels) != len(pvalues) : raise ValueError("Expected same length for pvalues and x_ticks.\n pvalues : {0} x_ticks : {1}".format(len(pvalues), len(xtickslabels)))
    if isinstance(pvalues,pd.DataFrame) :
        if len(pvalues.columns) != 1 : raise ValueError("Pass Dataframe with loc on pvalues (1 col expected)")
    elif isinstance(pvalues, pd.Series) :
        pass
    else :
        raise TypeError("Expected dataframe or Series")
    
    res_labels = []

    for label, pvalue in zip(xtickslabels, pvalues) :
        label = label.get_text()
        if pvalue <= significance :
            res_labels.append(
                label + "\n*"
            )
        else :
            res_labels.append(label)
    
    ax.set_xticklabels(res_labels)

    return ax

if DO_STEP3 :
    os.makedirs(PATH_OUT + '/proportion_RNA_focci/', exist_ok=True)

    columns = ['acquisition_id', 'treatment', 'rna1', 'rna2'] + ['rna1_number', 'rna1_clustered_number','rna2_number', 'rna2_clustered_number']
    data = Cell.loc[:, columns]


    data['rna1_proportion_rna_in_foci'] = data['rna1_clustered_number'] / data['rna1_number']
    data['rna2_proportion_rna_in_foci'] = data['rna2_clustered_number'] / data['rna2_number']
    data = data.dropna()
    assert len(data['treatment'].unique()) == EXPECTED_TREATMENT_NUMBER #+/- puro

    pvalues_rna1 = pvalues.compute_Ttest(
        df = data,
        group_keys= ["rna1","treatment"],
        measure_keys= ["rna1_proportion_rna_in_foci", "rna1_number"],

    )

    pvalues_rna2 = pvalues.compute_Ttest(
        df = data,
        group_keys= ["rna2","treatment"],
        measure_keys= ["rna2_proportion_rna_in_foci", "rna2_number"],

    )
    

    #Proportion Cy3 in foci
    sub_data = data.groupby(['rna1','treatment'])['rna1_proportion_rna_in_foci'].apply(list)
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data=sub_data,
        ax=ax,
        title= 'Proportion Cy3 in foci',
        xlabel= 'RNA',
        ylabel= 'fraction of RNA in foci'
    )
    ax.set_ylim(bottom= 0, top=1)
    ax = add_significance(ax, pvalues_rna1['rna1_proportion_rna_in_foci_pvalue'])

    fig.savefig(PATH_OUT + '/proportion_RNA_focci/proportion_Cy3_in_foci.png')
    plt.close()


    #Proportion Cy5 in foci
    sub_data = data.groupby(['rna1','treatment'])['rna1_proportion_rna_in_foci'].apply(list)
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data=sub_data,
        ax=ax,
        title= 'Proportion Cy3 in foci',
        xlabel= 'RNA',
        ylabel= 'fraction of RNA in foci'
    )
    ax.set_ylim(bottom= 0, top=1)
    ax = add_significance(ax, pvalues_rna2['rna1_proportion_rna_in_foci_pvalue'])
    fig.savefig(PATH_OUT + '/proportion_RNA_focci/proportion_Cy3_in_foci.png')
    plt.close()

    #Proportion Cy5 in foci
    sub_data = data.groupby(['rna2','treatment'])['rna2_proportion_rna_in_foci'].apply(list)
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data=sub_data,
        ax=ax,
        title= 'Proportion Cy5 in foci',
        xlabel= 'RNA',
        ylabel= 'fraction of RNA in foci'
    )
    ax.set_ylim(bottom= 0, top=1)
    ax = add_significance(ax, pvalues_rna2['rna2_proportion_rna_in_foci_pvalue'])
    fig.savefig(PATH_OUT + '/proportion_RNA_focci/proportion_Cy5_in_foci.png')
    plt.close()

    #Distributions

    os.makedirs(PATH_OUT + '/distributions/', exist_ok=True)
    
    #RNA 1 COUNT PER CELL
    data = Cell.groupby(['rna1','treatment'])['rna1_number'].apply(list).rename('single molecule count')
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data = data,
        ax=ax,
        title= "Number of Cy3 single molecule per cell",
        xlabel= 'Experiment date',
        ylabel= 'single molecule count'
    )

    ax = add_significance(ax, pvalues_rna1['rna1_number_pvalue'])
    fig.savefig(PATH_OUT + '/distributions/Cy3_count_per_cell.png')
    plt.close()

    #RNA 2 COUNT PER CELL
    data = Cell.groupby(['rna2','treatment'])['rna2_number'].apply(list).rename('single molecule count')
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data = data,
        ax=ax,
        title= "Number of Cy3 single molecule per cell",
        xlabel= 'RNA',
        ylabel= 'single molecule count'
    )

    ax = add_significance(ax, pvalues_rna2['rna2_number_pvalue'])
    fig.savefig(PATH_OUT + '/distributions/Cy5_count_per_cell.png')
    plt.close()

    #RNA 2 COUNT IN FOCI PER CELL
    data = Cell.groupby(['rna2','treatment'])['rna2_clustered_number'].apply(list).rename('single molecule count in foci')
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data = data,
        ax=ax,
        title= "Number of Cy5 single molecule in foci per cell",
        xlabel= 'RNA',
        ylabel= 'single molecule count'
    )

    fig.savefig(PATH_OUT + '/distributions/Cy5_in_foci_count_per_cell.png')
    plt.close()



"""

Step 4 : fraction of colocalising spot (Cy3 with Cy5 ~ rna1 vs rna2)

"""

if DO_STEP4 :
    mask = (Cell['rna2'].str.upper() != 'EMPTY')
    Cell_coloc = Cell[mask]
    Cell_coloc.columns = Cell_coloc.columns.str.replace('rna1','Cy3')
    Cell_coloc.columns = Cell_coloc.columns.str.replace('rna2','Cy5')

    pvalues_col = [
        "Cy3_Cy5_colocalisation_count",
        "Cy3_Cy5_free_colocalisation_count",
        "Cy3_Cy5_clustered_colocalisation_count",
        "Cy5_Cy3_colocalisation_count",
        "Cy5_Cy3_clustered_colocalisation_count",
        "Cy5_Cy3_free_colocalisation_count",
        ]
    
    
    #Plot for Cy3 colloc with all Cy5 populations
    cy3_coloc_dict = plot_coloc(
        dataframe=Cell_coloc,
        obj1_key='Cy3',
        obj2_key='Cy5',
        additional_grouping_keys= ['treatment'],
        do_clustered_pop=True,
    )

    #pvalues
    pvalues_df = pvalues.compute_Ttest(
        df = Cell_coloc,
        group_keys= ["Cy3","Cy5","treatment"],
        measure_keys= pvalues_col,
    )

    cy3_coloc_dict = add_pvalue(
        fig_dict=cy3_coloc_dict,
        pvalues=pvalues_df,
        do_clustered_pop=True
    )

    save_fig_dict(
        fig_dict=cy3_coloc_dict,
        obj_name="Cy3",
        output_folder=PATH_OUT
    )



    #Plot for Cy5 colloc with all Cy3 populations
    cy5_coloc_dict = plot_coloc(
        dataframe=Cell_coloc,
        obj1_key='Cy5',
        obj2_key='Cy3',
        additional_grouping_keys= ['treatment'],
        do_clustered_pop=True,
    )

    #pvalues
    pvalues_df = pvalues.compute_Ttest(
        df = Cell_coloc,
        group_keys= ["Cy5","Cy3","treatment"],
        measure_keys= pvalues_col,
    )

    cy5_coloc_dict = add_pvalue(
        fig_dict=cy5_coloc_dict,
        pvalues=pvalues_df,
        do_clustered_pop=True
    )

    save_fig_dict(
        fig_dict=cy5_coloc_dict,
        obj_name="Cy5",
        output_folder=PATH_OUT
    )