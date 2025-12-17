import os
import numpy as np
import pandas as pd
from typing import Literal
from itertools import combinations, product
from CustomPandasFramework.pvalues import compute_Ttest
import matplotlib.pyplot as plt
import warnings

def _get_title_sufix() :

    res = {
        '' : "(ALL)",
        '_free' : "(FREE)",
        '_clustered' : "(CLUSTERED)",
    }

    return res

def _get_axes_list(axes, do_clustered_pop, row) -> 'list[plt.Axes]' :
    if do_clustered_pop :
        return axes[row]
    else :
        return axes

def _get_sufixes(do_clustered_pop) :
    sufixes = [""]
    if do_clustered_pop : 
        sufixes.append("_clustered")
        sufixes.append("_free")
    
    return sufixes


def _coloc_get_target_columns(
        dataframe : pd.DataFrame, 
        obj1_key, 
        obj2_key,
        do_clustered_pop,
        ) :

    sufixes = _get_sufixes(do_clustered_pop)

    key_list = []

    for pop1, pop2 in product(sufixes,sufixes) :

        key1 = "{2}{0}_{3}{1}_colocalisation_count".format(pop1,pop2, obj1_key, obj2_key)
        key3 = "{2}{0}_{3}{1}_colocalisation_count".format(pop1,pop2, obj2_key, obj1_key)

        if key1 not in dataframe.columns : raise KeyError(key1)
        if key3 not in dataframe.columns : raise KeyError(key3)

        key_list += [key1,key3,]

    return key_list

def _compute_coloc_fractions(
        dataframe : pd.DataFrame,
        target_columns : 'list[str]',
        obj1_key,
        obj2_key,
        additional_grouping_keys,

) :
    columns = [
        obj1_key, 
        obj2_key, 
        "{0}_number".format(obj1_key),
        "{0}_number".format(obj2_key)
        ]
    
    if type(additional_grouping_keys) == list : columns += additional_grouping_keys

    fraction_frame = dataframe.loc[:, columns]

    new_columns = []

    for column in target_columns :
        population_key = column.split('_')[0]
        assert population_key == obj1_key or population_key == obj2_key
        new_col = "colocalisation_fraction_{0}".format(column.replace('_colocalisation_count', ''))
        fraction_frame[new_col] = dataframe[column] / dataframe["{0}_number".format(population_key)]
        if any(fraction_frame[new_col] > 1) : warnings.warn("fraction > 1 found in {0}".format(new_col))
        new_columns.append(new_col)

    return fraction_frame, new_columns

def _plot_pvalue_bar(ax: plt.Axes, x_sample1, y_sample1, x_sample2, y_sample2) :

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

def coloc_bar_plot(data: pd.DataFrame, ax: plt.Axes, couple:str, expected_treatment_number) :

    Z_SCORE = 2.575829 # 99% confidence interval

    
    data = data.melt(
        id_vars= 'treatment',
        value_vars= ['fraction_coloc', 'fraction_free_coloc','fraction_clustered_coloc']
    )

    treatments_names = data['treatment'].unique()
    assert len(treatments_names) <= expected_treatment_number #+/- puro

    data['variable'] = data['variable'].replace('fraction_coloc','all')
    data['variable'] = data['variable'].str.replace('fraction_','')
    data['variable'] = data['variable'].str.replace('_coloc','')
    data_list = data.groupby(['variable','treatment'])['value'].apply(list)

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

    ax.set_xticks(xticks_pos, treatments_names)
    ax.legend()
    ax.set_ylim(bottom= 0)

    return ax


def _plot_bars(
        ax : plt.Axes,
        sub_data : pd.DataFrame,
        targets,
        obj1,
        obj2,
        obj1_key,
        do_clustered_pop
) :
    
    ax_data = sub_data.xs(obj2, axis=0, level=1, drop_level=False)
    labels = ax_data.index.get_level_values(-1)
    xticks_pos = np.arange(len(labels))
    x_range = [xticks_pos[0] - 1, xticks_pos[-1] + 1]

    ax.bar(
        x= xticks_pos,
        height= ax_data.xs(targets['target_column_all'], axis=1, level=0)['mean'],
        yerr= ax_data.xs(targets['target_column_all'], axis=1, level=0)['error'],
        capsize=4,
        width= 0.8,
        color= 'magenta',
        label= '{0} all'.format(obj1_key),
        alpha= 0.5
        )
    
    if do_clustered_pop :
        #free pop bars
        ax.bar(
            x= xticks_pos -0.2,
            height= ax_data.xs(targets['target_column_free'], axis=1, level=0)['mean'],
            yerr= ax_data.xs(targets['target_column_free'], axis=1, level=0)['error'],
            capsize=4,
            align= 'center',
            width= 0.3,
            color= 'red',
            label= '{0} free'.format(obj1_key),
            alpha= 0.75
        )

        #clustered pop bars
        ax.bar(
            x= xticks_pos + 0.2,
            height= ax_data.xs(targets['target_column_clustered'], axis=1, level=0)['mean'],
            yerr= ax_data.xs(targets['target_column_clustered'], axis=1, level=0)['error'],
            capsize=4,
            align= 'center',
            width= 0.3,
            color= 'blue',
            label= 'clusted all'.format(obj1_key),
            alpha = 0.75
        ) 
        return ax

def plot_coloc(
        dataframe : pd.DataFrame,
        obj1_key : Literal["rna1", "rna2", "suntag"],
        obj2_key : Literal["rna1", "rna2", "suntag"],
        output_folder : str,
        additional_grouping_keys : 'list[str]' = None,
        do_clustered_pop = False,
        obj1_col = None,
        obj2_col = None,
        figsize = (10,10),
        Z_SCORE = 2.575829, # 99% confidence interval


) :
    
    #Integrity of columns
    if obj1_key == obj2_key : raise ValueError("obj1_key and obj2_key must be different.")
    if obj1_key not in ["rna1", "rna2", "suntag"] : raise ValueError('obj1_key must be of the following : ["rna1", "rna2", "suntag"]')
    if obj2_key not in ["rna1", "rna2", "suntag"] : raise ValueError('obj2_key must be of the following : ["rna1", "rna2", "suntag"]')
    target_columns = _coloc_get_target_columns(
        dataframe=dataframe,
        obj1_key=obj1_key,
        obj2_key=obj2_key,
        do_clustered_pop=do_clustered_pop,
    )

    #Optional if col name doesn't correspond to obj key
    if type(obj1_col) == str :
        if obj1_col not in dataframe.columns : raise ValueError("{0} column not found.".format(obj1_col))
        dataframe = dataframe.rename(columns={obj1_col : obj1_key})

    if type(obj2_col) == str :
        if obj2_col not in dataframe.columns : raise ValueError("{0} column not found.".format(obj2_col))
        dataframe = dataframe.rename(columns={obj2_col : obj2_key})

    #Preparing suffix and title legend
    sufixes = _get_sufixes(do_clustered_pop)
    title_suffix = _get_title_sufix()
    
    #Computing colocalisation fractions
    fraction_frame, target_columns = _compute_coloc_fractions(
        dataframe=dataframe,
        target_columns=target_columns,
        obj1_key=obj1_key,
        obj2_key=obj2_key,
        additional_grouping_keys=additional_grouping_keys
    )

    #Grouping data
    fraction_frame = fraction_frame.groupby([obj1_key, obj2_key] + additional_grouping_keys).agg(
        dict.fromkeys(target_columns, ['mean','std','count'])
    )

    #Computing confidence interval as 'error' key
    error_bar = fraction_frame.xs('std',axis=1,level=1) * Z_SCORE / fraction_frame.xs('count',axis=1,level=1)
    error_bar.columns = pd.MultiIndex.from_product([error_bar.columns, ['count']])
    fraction_frame.loc[:,(slice(None),"count")] = error_bar
    fraction_frame = fraction_frame.rename({'count' : 'error'}, axis=1)

    #Plotting
    res_dict = {}
    obj1_list = list(fraction_frame.index.get_level_values(0).unique())
    couple_list = list(zip(fraction_frame.index.get_level_values(0),fraction_frame.index.get_level_values(1)))

    subplot_rows_number = 1 + 2*do_clustered_pop
    for obj1 in obj1_list :

        #create folder
        plot_path = output_folder + "/colocalization/{0}/".format(obj1_key)
        os.makedirs(plot_path, exist_ok=True)

        #Select data for group of plots.
        #FIGURE LEVEL
        sub_data = fraction_frame.loc[[obj1],:]
        obj2_list = sub_data.index.get_level_values(1).unique()
        subplot_col_number = len(obj2_list)
        if type(additional_grouping_keys) != type(None) :
            conditions_list = list(sub_data.index.get_level_values(-1).unique())
            bar_number = len(conditions_list)
        else : 
            conditions_list = []
            bar_number = 0
            
        #Start plot
        fig = plt.figure(figsize= (subplot_col_number*figsize[0], subplot_rows_number*figsize[1]))
        axes = fig.subplots(nrows=subplot_rows_number, ncols= subplot_col_number, squeeze=False, sharey=True)

        #ROW LEVEL
        for row, obj2_sufix in enumerate(sufixes) :
            ax_list = _get_axes_list(axes, do_clustered_pop, row)

            targets = {
                'target_column_all' : "colocalisation_fraction_{0}_{1}{2}".format(obj1_key, obj2_key,obj2_sufix),
                'target_column_free' : "colocalisation_fraction_{0}_free_{1}{2}".format(obj1_key, obj2_key,obj2_sufix),
                'target_column_clustered' : "colocalisation_fraction_{0}_clustered_{1}{2}".format(obj1_key, obj2_key,obj2_sufix),
                'target_column_all_sym' : "colocalisation_fraction_{0}{1}_clustered_{2}".format(obj2_key, obj2_sufix, obj1_key),
                'target_column_free_sym' : "colocalisation_fraction_{0}{1}_clustered_{2}".format(obj2_key, obj2_sufix, obj1_key),
                'target_column_clustered_sym' : "colocalisation_fraction_{0}{1}_clustered_{2}".format(obj2_key, obj2_sufix, obj1_key),
            }
                
            for obj2, ax in zip(obj2_list, ax_list) :

                has_symetric = (obj2,obj1) in couple_list

                ax_data = sub_data.xs(obj2, axis=0, level=1, drop_level=False)
                labels = ax_data.index.get_level_values(-1)
                xticks_pos = np.arange(len(labels))
                x_range = [xticks_pos[0] - 1, xticks_pos[-1] + 1]

                _plot_bars(
                    ax=ax,
                    sub_data=sub_data,
                    targets=targets
                ) 

                ax.set_xlabel(ax_data.index.names[-1])
                ax.set_ylabel("colocalization fraction")
                ax.set_title("{0} --> {1} {2}".format(obj1,obj2, title_suffix[obj2_sufix]))
                ax.set_xticks(xticks_pos,labels)
                ax.set_xlim(left=x_range[0], right=x_range[1])
                ax.legend()
                    
        fig.subplots_adjust(hspace=0.4)
        fig.savefig(plot_path + "/{0}_colocalization.png".format(obj1))
        res_dict[obj1] = fig

    quit()
    return res_dict