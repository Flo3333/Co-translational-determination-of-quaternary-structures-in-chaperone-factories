"""
Analysis script for Soha quantifications
Pipeline : Fish_seq / BC_clusterkinetics_pipeline

You have tu run first the merge/cleaning script : soha_quantification_merge.py
This analysis restricts cells (counted as phenotype_positive) to cell showing numerous clusters from APC channel (rna1).

Notes : 

→ proportion ARN qui colocalisent APC(Vert) ARN(Orange)

→ proportion ARN focci (Puro vs non puro)

→ Distribution intensité des spots APC/ARN

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pbwrap.plot as plot
import pbwrap.quantification.statistical_test as stest
from itertools import chain

FIGSIZE = (10,10)
BINS = 100

DO_STEP2 = True # Spots intensity distribution
DO_STEP3 = True # fraction rna in foci
DO_STEP4 = True # proportion spots colocalisation
EXPECTED_TREATMENT_NUMBER = 2

PATH_IN= "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/merge_folder/"
PATH_OUT= "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/analysis/"

Z_SCORE = 2.575829 # 99% confidence interval

"""
STEP1 : Open data
"""
Acquisition = pd.read_feather(PATH_IN + '/Acquisition.feather')
Cell = pd.read_feather(PATH_IN + '/Cell.feather')
Spots = pd.read_feather(PATH_IN + '/Spots.feather')

for table in [Acquisition, Cell, Spots] :
    table['acquisition_id'] = table['acquisition_id'].apply(tuple)

"""
STEP 2 : Spots intensity distribution / Mean rna in foci number distribution
"""

def intensity_histogram_plot(data:pd.DataFrame, logscale = False) :
    fig = plt.figure(figsize=(FIGSIZE[0], FIGSIZE[1]*2))
    axes = fig.subplots(2,1)
    treatments = data['treatment'].unique()

    assert len(treatments) == EXPECTED_TREATMENT_NUMBER # +/- puromycin
    colors = ['green', 'blue']

    for ax, treatment, color in zip(axes, treatments, colors) :
        sub_data = data.loc[data['treatment'] == treatment]
        ax.hist(sub_data['intensity'], bins=BINS, label= treatment, alpha= 0.7, color=color)
        ax.set_xlabel('spot intensity')
        ax.set_ylabel('count')
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        ax.set_title('Spot number {0}'.format(len(sub_data)))
        ax.legend()
    
    #Same scale
    xmin1,xmax1,ymin1,ymax1 = axes[0].axis()
    xmin2,xmax2,ymin2,ymax2 = axes[1].axis()
    xmin = min(xmin1,xmin2)
    xmax = max(xmax1,xmax2)
    if logscale :
        ymin1 = ymin2 = 1
        for ax in axes : ax.set_yscale('log') 

    axes[0].axis([xmin,xmax,ymin1,ymax1])
    axes[1].axis([xmin,xmax,ymin2,ymax2])

    return fig



if DO_STEP2 :
    os.makedirs(PATH_OUT + '/spot_intensity_histogram/', exist_ok=True)

    #Merge for treatment
    Spots = pd.merge(Spots, Acquisition.loc[:,['acquisition_id', 'treatment']], how='left', on= 'acquisition_id')
    Spots['spot_type'] = Spots['spot_type'].replace({
        'rna1' : 'APC',
        'rna2' : 'RNA'
    })
    spot_types = Spots['spot_type'].unique()

    for spot_type in spot_types :
        data = Spots.loc[Spots['spot_type'] == spot_type]

        fig = intensity_histogram_plot(data)
        fig.suptitle(spot_type)
        fig.savefig(PATH_OUT + '/spot_intensity_histogram/{0}_spot_intensity_histogram.png'.format(spot_type))
        plt.close()
        fig = intensity_histogram_plot(data, logscale=True)
        fig.suptitle(spot_type)
        fig.savefig(PATH_OUT + '/spot_intensity_histogram/{0}_logscale_spot_intensity_histogram.png'.format(spot_type))
        plt.close()

    

del Spots

"""
Cell merge
"""

if 'rna2' in Cell.columns : Cell = Cell.drop(columns='rna2')

Cell = pd.merge(
    Acquisition.loc[:,['acquisition_id', 'rna1', 'rna2','date']],
    Cell,
    on= 'acquisition_id',
    how='right'
)


"""

Step 3 :  proportion ARN focci (Puro vs non puro)

"""

def mean_rna_per_foci(data:pd.DataFrame, logscale = False) :
    fig = plt.figure(figsize=(FIGSIZE[0], FIGSIZE[1]*2))
    axes = fig.subplots(2,1)
    treatments = data['treatment'].unique()

    data['mean_rna_per_foci'] = data['rna2_clustered_number'] / data['rna2_cluster_number']
    data['mean_rna_per_foci'] = data['mean_rna_per_foci'].replace(np.inf, np.NaN)

    assert len(treatments) <= EXPECTED_TREATMENT_NUMBER # +/- puromycin
    colors = ['green', 'blue']

    for ax, treatment, color in zip(axes, treatments, colors) :
        sub_data = data.loc[data['treatment'] == treatment]
        ax.hist(sub_data['mean_rna_per_foci'], bins=10, label= treatment, alpha= 0.7, color=color, density=True)
        ax.set_xlabel('mean rna per foci')
        ax.set_ylabel('distribution density')
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        ax.set_title('Cell number {0}'.format(len(sub_data)))
        ax.legend()
    
    #Same scale
    xmin1,xmax1,ymin1,ymax1 = axes[0].axis()
    xmin2,xmax2,ymin2,ymax2 = axes[1].axis()
    xmin = min(xmin1,xmin2)
    xmax = max(xmax1,xmax2)
    if logscale :
        ymin1 = ymin2 = 1
        for ax in axes : ax.set_yscale('log') 

    axes[0].axis([xmin,xmax,ymin1,ymax1])
    axes[1].axis([xmin,xmax,ymin2,ymax2])

    return fig

if DO_STEP3 :

    os.makedirs(PATH_OUT + '/mean_rna_per_foci/', exist_ok=True)
    bcat_ARN_idx = Cell.query("rna2 == 'B-CAT-ARN'").index
    fig = mean_rna_per_foci(Cell.loc[bcat_ARN_idx], logscale=False)
    fig.savefig(PATH_OUT + '/mean_rna_per_foci/logscale_mean_b-cat_RNA_per_foci_histogram.png')
    plt.close()

    os.makedirs(PATH_OUT + '/proportion_ARN_focci/', exist_ok=True)

    columns = ['acquisition_id', 'treatment', 'rna1', 'rna2'] + ['rna1_number', 'rna1_clustered_number','rna2_number', 'rna2_clustered_number']
    data = Cell.loc[:, columns]


    data['rna1_proportion_rna_in_foci'] = data['rna1_clustered_number'] / data['rna1_number']
    data['rna2_proportion_rna_in_foci'] = data['rna2_clustered_number'] / data['rna2_number']

    data = pd.concat([
        data.loc[:,['treatment','rna1', 'rna1_proportion_rna_in_foci']].rename(columns={'rna1' : 'rna', 'rna1_proportion_rna_in_foci' : 'proportion_rna_in_foci'}),
        data.loc[:,['treatment','rna2', 'rna2_proportion_rna_in_foci']].rename(columns={'rna2' : 'rna', 'rna2_proportion_rna_in_foci' : 'proportion_rna_in_foci'}),
    ],
    axis=0)

    data = data.dropna()
    assert len(data['treatment'].unique()) == EXPECTED_TREATMENT_NUMBER #+/- puro

    data = data.groupby(['rna','treatment'])['proportion_rna_in_foci'].apply(list)

    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data=data,
        ax=ax,
        title= 'Proportion RNA in foci',
        xlabel= 'RNA',
        ylabel= 'fraction of RNA in foci'
    )
    ax.set_ylim(bottom= 0, top=1)

    fig.savefig(PATH_OUT + '/proportion_ARN_focci/proportion_ARN_in_foci.png')
    plt.close()

    #Distributions

    os.makedirs(PATH_OUT + '/distributions/', exist_ok=True)
    
    #RNA 1 COUNT PER CELL
    data = Cell.groupby(['date','treatment'])['rna1_number'].apply(list).rename('single molecule count')
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data = data,
        ax=ax,
        title= "Number of APC single molecule per cell",
        xlabel= 'Experiment date',
        ylabel= 'single molecule count'
    )

    fig.savefig(PATH_OUT + '/distributions/APC_count_per_cell.png')
    plt.close()

    #RNA 2 COUNT PER CELL
    data = Cell.groupby(['rna2','treatment'])['rna2_number'].apply(list).rename('single molecule count')
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data = data,
        ax=ax,
        title= "Number of RNA single molecule per cell",
        xlabel= 'RNA',
        ylabel= 'single molecule count'
    )

    fig.savefig(PATH_OUT + '/distributions/RNA_count_per_cell.png')
    plt.close()

    #RNA 2 COUNT IN FOCI PER CELL
    data = Cell.groupby(['rna2','treatment'])['rna2_clustered_number'].apply(list).rename('single molecule count in foci')
    fig = plt.figure(figsize= FIGSIZE)
    ax = fig.gca()
    ax = plot.distribution_super_plot(
        data = data,
        ax=ax,
        title= "Number of RNA single molecule in foci per cell",
        xlabel= 'RNA',
        ylabel= 'single molecule count'
    )

    fig.savefig(PATH_OUT + '/distributions/RNA_in_foci_count_per_cell.png')
    plt.close()




"""

Step 4 : fraction of colocalising spot (APC with ARN ~ rna1 vs rna2)

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
            samples1= data_list.loc[('all', 'puromycin')],
            samples2= data_list.loc[('all', 'untreated')],
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
        label= 'ARN all',
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
        label= 'ARN free',
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
        label= 'ARN clustered',
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
    os.makedirs(PATH_OUT + '/colocalisation/', exist_ok=True)

    #Plot for ARN colloc with all APC population
    columns = ['treatment','rna1', 'rna2','rna2_number', 'rna2_free_number', 'rna2_clustered_number', 'rna2_rna1_colocalisation_count', 'rna2_free_rna1_colocalisation_count','rna2_clustered_rna1_colocalisation_count']
    data = Cell.loc[:,columns]
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
        ax.set_title("ARN : {0} \nAPC : {1}".format(couple[1], couple[0]))

        if couple != ('B-CAT-APC', 'B-CAT-ARN') :
            
            min_x,max_x,min_y, max_y = ax.axis()
            miny = min(miny, min_y)
            maxy = max(maxy, max_y)

            if i > 0 : #meaning is not the leftmost plot.
                ax.axes.get_yaxis().set_visible(False)
                ax.spines['left'].set_visible(False)
        i += 1
    
    for couple, ax in zip(coloc_couples, axes) :
        if couple == ('B-CAT-APC', 'B-CAT-ARN') : continue
        min_x,max_x,min_y, max_y = ax.axis()
        ax.axis([min_x,max_x,miny, maxy])

    fig.suptitle("Colocalisation ARN with APC (all)")
    fig.savefig(PATH_OUT + '/colocalisation/fraction_coloc_APC_ARN.png')
    plt.close()

    #Plot for ARN colloc with APC clustered

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
            if line_index == 0 : ax.set_title("ARN : {0} \nAPC : {1}".format(couple[1], couple[0]), fontdict={'size' : 20})

            if couple != ('B-CAT-APC', 'B-CAT-ARN') : 

                min_x,max_x,min_y, max_y = ax.axis()
                miny = min(miny, min_y)
                maxy = max(maxy, max_y)

            if i > 0 : #meaning is not the leftmost plot.
                if couple != ('B-CAT-APC', 'B-CAT-ARN') : 
                    ax.axes.get_yaxis().set_visible(False)
                    ax.spines['left'].set_visible(False)
            else :
                ax.set_ylabel("with {0}".format(target), fontdict= {"size" : 20})

            i += 1

        for couple, ax in zip(coloc_couples, axes[line_index*4 : line_index*4+4]) :
            if couple == ('B-CAT-APC', 'B-CAT-ARN') : continue
            min_x,max_x,min_y, max_y = ax.axis()
            final_axis = ax.axis([min_x,max_x,miny, maxy])

        fig.suptitle("Colocalisation Fractions ARN with\nAPC cluster (top)\nAPC free spots(bottom)", size= 20 )
        fig.savefig(PATH_OUT + '/colocalisation/fraction_coloc_ARN_APC_clusters_or_free.png')
        plt.close()