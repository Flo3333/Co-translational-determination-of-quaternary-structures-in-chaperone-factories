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
from CustomPandasFramework.Fish_seq.plot import plot_coloc
from itertools import chain


FIGSIZE = (15,15)
BINS = 100

DO_STEP2 = False #< Number of cell quantified
DO_STEP3 = False # fraction rna in foci
DO_STEP4 = True # proportion spots colocalisation
EXPECTED_TREATMENT_NUMBER = 2

PATH_IN= "/home/flo/Documents/scp/merge_folder/"
PATH_OUT= "/home/flo/Documents/scp/analysis/"

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
    quit()



"""

Step 4 : fraction of colocalising spot (Cy3 with Cy5 ~ rna1 vs rna2)

"""

def plot_pvalue_bar(ax: plt.Axes, x_sample1, y_sample1, x_sample2, y_sample2) :

    y_ticks = ax.get_yticks()
    if len(y_ticks) >= 2 : y_step = y_ticks[1] - y_ticks[0]
    else : y_step = 0

    y_sample1 +=y_step*0.2
    y_sample2 +=y_step*0.2

    ceiling = max(y_sample1 + y_step , y_sample2 + y_step,)

    ax.plot([x_sample1 - 0.1, x_sample1 + 0.1, x_sample1,x_sample1,x_sample2, x_sample2, x_sample2 - 0.1, x_sample2 + 0.1], [y_sample1, y_sample1, y_sample1, ceiling, ceiling, y_sample2, y_sample2, y_sample2], 'r')
    ax.text(
        x = x_sample1 + (x_sample2 - x_sample1)/2,
        y= ceiling,
        s= "*",
        c='red',
        fontdict={'size' : 20}
        )

    return ax


def coloc_bar_plot(data: pd.DataFrame, ax: plt.Axes, couple:str) :

    data = data.dropna()
    
    data = data.melt(
        id_vars= 'treatment',
        value_vars= ['fraction_coloc', 'fraction_free_coloc','fraction_clustered_coloc']
    )

    treatments_names = data['treatment'].unique()
    assert len(treatments_names) <= EXPECTED_TREATMENT_NUMBER #+/- puro

    data['variable'] = data['variable'].replace('fraction_coloc','all')
    data['variable'] = data['variable'].str.replace('fraction_','')
    data['variable'] = data['variable'].str.replace('_coloc','')
    data_list = data.groupby(['variable','treatment'])['value'].apply(list)

    if  len(treatments_names) == 2 :
        pvalue = stest.t_test(
            samples1= data_list.loc[('all', 'puro')],
            samples2= data_list.loc[('all', 'unt')],
            equal_variance=False
        )[0]
    else : 
        pvalue = np.NaN
    data_mean = data.groupby(['variable','treatment'])['value'].mean()
    data_std = data.groupby(['variable','treatment'])['value'].std()
    data_sample_size_sqrt = data_list.apply(len).apply(np.sqrt)
    data_std = (data_std * Z_SCORE) / data_sample_size_sqrt

    xticks_pos = np.arange(1, 1+len(treatments_names)*4, 4)


    #all pop bars
    ax.bar(
        x= xticks_pos,
        height= data_mean.loc['all'],
        yerr= data_std['all'],
        capsize= 4,
        align= 'center',
        width= 2,
        color= 'magenta',
        label= 'Cy5 all',
        alpha= 0.5
    )
    #free pop bars
    ax.bar(
        x= xticks_pos -0.5,
        height= data_mean.loc['free'],
        yerr= data_std['free'],
        capsize= 4,
        align= 'center',
        width= 0.75,
        color= 'red',
        label= 'Cy5 free',
        alpha= 0.75
    )

    #clustered pop bars
    ax.bar(
        x= xticks_pos +0.5,
        height= data_mean.loc['clustered'],
        yerr= data_std['clustered'],
        capsize= 4,
        align= 'center',
        width= 0.75,
        color= 'blue',
        label= 'Cy5 clustered',
        alpha = 0.75
    )


    if pvalue <= 0.01 :
        plot_pvalue_bar(
            ax,
            xticks_pos[0],
            (data_mean.loc['all'].iat[0] + data_std.loc['all'].iat[0]),
            xticks_pos[1],
            (data_mean.loc['all'].iat[1] + data_std.loc['all'].iat[1]),

        )
    

    ax.set_xticks(xticks_pos, treatments_names)
    ax.legend()
    ax.set_ylim(bottom= 0)

    return ax


if DO_STEP4 :

    #Plot for Cy5 colloc with all Cy3 population
    plot_coloc(
        dataframe=Cell,
        obj1_key='rna1',
        obj2_key='rna2',
        output_folder=PATH_OUT,
        additional_grouping_keys= ['treatment'],
        do_clustered_pop=True,
    )

    quit()






    
    
    data['fraction_coloc'] = data['rna2_rna1_colocalisation_count'] / data['rna2_number']
    data['fraction_free_coloc'] = data['rna2_free_rna1_colocalisation_count'] / data['rna2_number']
    data['fraction_clustered_coloc'] = data['rna2_clustered_rna1_colocalisation_count'] / data['rna2_number']

    coloc_couples = data.value_counts(subset=['rna1','rna2']).index.unique()
    fig = plt.figure(figsize=(40,10))
    axes = fig.subplots(1,4)
    # axes = chain(*axes) #flatten list
    
    i=0
    miny = 0 # will set the global scale for not b-cat RNA
    maxy = 0
    for couple, ax in zip(coloc_couples, axes) :
        rna1,rna2 = couple
        sub_data_idx = data.query("rna1 == '{0}' and rna2 == '{1}'".format(rna1,rna2)).index
        sub_data = data.loc[sub_data_idx]

        ax = coloc_bar_plot(
            data=sub_data,
            ax=ax,
            couple=couple
        )

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title("Cy5 : {0} \nCy5 : {1}".format(couple[1], couple[0]))

        if couple != ('B-CAT-Cy5', 'B-CAT-Cy5') :
            
            min_x,max_x,min_y, max_y = ax.axis()
            miny = min(miny, min_y)
            maxy = max(maxy, max_y)

            if i > 0 : #meaning is not the leftmost plot.
                ax.axes.get_yaxis().set_visible(False)
                ax.spines['left'].set_visible(False)
        i += 1
    
    for couple, ax in zip(coloc_couples, axes) :
        if couple == ('B-CAT-Cy5', 'B-CAT-Cy5') : continue
        min_x,max_x,min_y, max_y = ax.axis()
        ax.axis([min_x,max_x,miny, maxy])

    fig.suptitle("Colocalisation Cy5 with Cy3 (all)")
    fig.savefig(PATH_OUT + '/colocalisation/fraction_coloc_Cy5_Cy5.png')
    plt.close()

    #Plot for Cy5 colloc with Cy3 clustered

    fig = plt.figure(figsize=(40,20))
    axes = fig.subplots(2,4)
    axes = list(chain(*axes)) #flatten list

    for target, line_index in zip(['clustered', 'free'], [0,1]) :

        columns = ['treatment','rna1', 'rna2','rna2_number', 'rna2_free_number', 'rna2_clustered_number', 'rna2_rna1_{0}_colocalisation_count'.format(target), 'rna2_free_rna1_{0}_colocalisation_count'.format(target),'rna2_clustered_rna1_{0}_colocalisation_count'.format(target)]
        data = Cell.loc[:,columns]
        data['fraction_coloc'] = data['rna2_rna1_{0}_colocalisation_count'.format(target)] / data['rna2_number']
        data['fraction_free_coloc'] = data['rna2_free_rna1_{0}_colocalisation_count'.format(target)] / data['rna2_number']
        data['fraction_clustered_coloc'] = data['rna2_clustered_rna1_{0}_colocalisation_count'.format(target)] / data['rna2_number']

        coloc_couples = data.value_counts(subset=['rna1','rna2']).index.unique()

        miny = 0 # will set the global scale for not b-cat RNA
        maxy = 0
        i=0
        for couple, ax in zip(coloc_couples, axes[line_index*4 : line_index*4+4]) :
            rna1,rna2 = couple
            sub_data_idx = data.query("rna1 == '{0}' and rna2 == '{1}'".format(rna1,rna2)).index
            sub_data = data.loc[sub_data_idx]

            ax = coloc_bar_plot(
                data=sub_data,
                ax=ax,
                couple=couple
            )

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            if line_index == 0 : ax.set_title("Cy5 : {0} \nCy3 : {1}".format(couple[1], couple[0]), fontdict={'size' : 20})

            if couple != ('B-CAT-Cy5', 'B-CAT-Cy5') : 

                min_x,max_x,min_y, max_y = ax.axis()
                miny = min(miny, min_y)
                maxy = max(maxy, max_y)

            if i > 0 : #meaning is not the leftmost plot.
                if couple != ('B-CAT-Cy5', 'B-CAT-Cy5') : 
                    ax.axes.get_yaxis().set_visible(False)
                    ax.spines['left'].set_visible(False)
            else :
                ax.set_ylabel("with {0}".format(target), fontdict= {"size" : 20})

            i += 1

        for couple, ax in zip(coloc_couples, axes[line_index*4 : line_index*4+4]) :
            if couple == ('B-CAT-Cy5', 'B-CAT-Cy5') : continue
            min_x,max_x,min_y, max_y = ax.axis()
            final_axis = ax.axis([min_x,max_x,miny, maxy])

        fig.suptitle("Colocalisation Fractions Cy5 with\nCy5 cluster (top)\nCy5 free spots(bottom)", size= 20 )
        fig.savefig(PATH_OUT + '/colocalisation/fraction_coloc_Cy5_Cy5_clusters_or_free.png')
        plt.close()