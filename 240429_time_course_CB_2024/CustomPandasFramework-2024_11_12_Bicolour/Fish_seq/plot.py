import os
import numpy as np
import pandas as pd
from typing import Literal
from itertools import combinations, product
from CustomPandasFramework.pvalues import compute_Ttest
from matplotlib.offsetbox import AnchoredText
from pbwrap.plot import add_pvalue_star
import matplotlib.pyplot as plt
import warnings

# Single molecules in foci
def compute_fraction_in_foci(
        data : pd.DataFrame,
        keys : list,
        ) :
    
    for key in keys : 
        data['{0}_proportion_in_foci'.format(key)] = data['{0}_clustered_number'.format(key)] / data['{0}_number'.format(key)]
    
    return data


# plot_coloc

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

def _get_symetric_couples(dataframe : pd.DataFrame) :
    dataframe = dataframe.swaplevel(0,1,axis=0).sort_index()
    couples_list = pd.unique(list(zip(
        dataframe.index.get_level_values(0),
        dataframe.index.get_level_values(1),
    )))
    return list(couples_list)

def _harmonize_sym_data(ax_data : pd.DataFrame, ax_data_sym: pd.DataFrame) -> 'tuple[pd.DataFrame,pd.DataFrame]':
    """
    Handles case where ax_data and ax_data_sym have different index because of 'additional grouping key' in groupby. Such case could happen if different drugs have been used for symmetric couple.
    """
    
    ax_data = ax_data.sort_index()
    ax_data_sym = ax_data_sym.sort_index().swaplevel(0,1)

    if ax_data.index.equals(ax_data_sym.index) :
        return ax_data, ax_data_sym
    else :

        full_index : pd.Index = ax_data.index.append(ax_data_sym.index)
        full_index = full_index.unique().sort_values()
        ax_data = ax_data.reindex(index=full_index, method=None)
        ax_data_sym = ax_data_sym.reindex(index=full_index, method=None)

        return ax_data, ax_data_sym.swaplevel(0,1)
        
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

    # wfraction_frame = fraction_frame.groupby([obj1_key, obj2_key] + additional_grouping_keys).agg(
    #         dict.fromkeys(target_columns, ['mean','std','count'])
    #     )

    new_columns = []

    for column in target_columns :
        population_key = column.split('_')[0]
        assert population_key == obj1_key or population_key == obj2_key
        new_col = "colocalisation_fraction_{0}".format(column.replace('_colocalisation_count', ''))
        fraction_frame[new_col] = dataframe[column] / dataframe["{0}_number".format(population_key)]
        if any(fraction_frame[new_col] > 1) : warnings.warn("fraction > 1 found in {0}".format(new_col))
        new_columns.append(new_col)


    return fraction_frame, new_columns

def _plot_bars(
        ax : plt.Axes,
        ax_data : pd.DataFrame,
        targets,
        obj1_key,
        obj2_key,
        do_clustered_pop,
        symetric_plot,
        has_sym,
        ax_data_sym: pd.DataFrame=None,
) :
    
    assert (has_sym and type(ax_data_sym != type(None)) or (not has_sym and type(ax_data_sym == type(None)))), "If data has symetric rna pair ax_data_sym should be passed, else it shouldn't."
    
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
            label= '{0} clustered'.format(obj1_key),
            alpha = 0.75
        )

    if has_sym :
        ax.bar(
        x= xticks_pos,
        height= -ax_data_sym.xs(targets['target_column_all_sym'], axis=1, level=0)['mean'],
        yerr= ax_data.xs(targets['target_column_all_sym'], axis=1, level=0)['error'],
        capsize=4,
        width= 0.8,
        color= 'magenta',
        label= None,
        alpha= 0.5
        )

        if do_clustered_pop :
            #free pop bars
            ax.bar(
                x= xticks_pos -0.2,
                height= -ax_data_sym.xs(targets['target_column_free_sym'], axis=1, level=0)['mean'],
                yerr= ax_data.xs(targets['target_column_free_sym'], axis=1, level=0)['error'],
                capsize=4,
                align= 'center',
                width= 0.3,
                color= 'red',
                label= None,
                alpha= 0.75
            )

            #clustered pop bars
            ax.bar(
                x= xticks_pos + 0.2,
                height= -ax_data_sym.xs(targets['target_column_clustered_sym'], axis=1, level=0)['mean'],
                yerr= ax_data.xs(targets['target_column_clustered_sym'], axis=1, level=0)['error'],
                capsize=4,
                align= 'center',
                width= 0.3,
                color= 'blue',
                label= None,
                alpha = 0.75
            )


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if has_sym :
        left, right, bottom, top  = ax.axis()
        bottom = top = max(
            abs(bottom),
            abs(top), 
            ax_data.xs(targets['target_column_all'], axis=1, level=0)['mean'].max(),
            ax_data_sym.xs(targets['target_column_free_sym'], axis=1, level=0)['mean'].max(),
            ax_data_sym.xs(targets['target_column_all_sym'], axis=1, level=0)['mean'].max(),
            ax_data_sym.xs(targets['target_column_clustered_sym'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_clustered'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_free'], axis=1, level=0)['mean'].max(),
            )

        ax.set_ylim(bottom= -bottom*1.1, top= top*1.1)
        ax.plot([left,right], [0,0], '--k',label=None, alpha =0.7)
        ax.axhspan(-1, 0, facecolor='orange', alpha=0.3)

    elif symetric_plot:
        left, right, bottom, top  = ax.axis()
        bottom = top = max(
            abs(bottom),
            abs(top), 
            ax_data.xs(targets['target_column_all'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_clustered'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_free'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_all'], axis=1, level=0)['mean'].max(),
            )
        ax.set_ylim(bottom= -bottom*1.1, top= top*1.1)
        ax.axhspan(-1, 0, facecolor='gray', alpha=0.3)
        ax.plot([left,right], [0,0], 'k',label=None, alpha =0.7)
        

    else :
        left, right, bottom, top  = ax.axis()
        top = max(
            abs(top), 
            ax_data.xs(targets['target_column_all'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_clustered'], axis=1, level=0)['mean'].max(),
            ax_data.xs(targets['target_column_free'], axis=1, level=0)['mean'].max(),
            )
        ax.set_ylim(bottom=0, top=top)
        ax.axhspan(-1, 0, facecolor='gray', alpha=0.3)
        ax.plot([left,right], [0,0], 'k',label=None, alpha =0.7)


    return ax

def plot_coloc(
        dataframe : pd.DataFrame,
        obj1_key : str,
        obj2_key : str,
        additional_grouping_keys : 'list[str]' = None,
        do_clustered_pop = False,
        obj1_col = None,
        obj2_col = None,
        figsize = (10,10),
        Z_SCORE = 2.575829, # 99% confidence interval,
        verbose = False,


) -> 'dict[plt.Figure]':
    """
    Return dict with key = obj1_obj2 containing matplotlib.figure object ready to be saved or for further customisation.
    """
    
    #Integrity of columns
    if obj1_key == obj2_key : raise ValueError("obj1_key and obj2_key must be different.")
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
    res_dict = {} # init return object
    obj1_list = list(fraction_frame.index.get_level_values(0).unique())
    sym_couple_list = _get_symetric_couples(fraction_frame)

    subplot_rows_number = 1 + 2*do_clustered_pop
    for obj1 in obj1_list :

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
        symetric_plot = False #Init symetric plot
        fig = plt.figure(figsize= (subplot_col_number*figsize[0], subplot_rows_number*figsize[1]))
        axes = fig.subplots(nrows=subplot_rows_number, ncols= subplot_col_number, squeeze=False, sharey=True)

        #ROW LEVEL
        for row, obj2_sufix in enumerate(sufixes) :
            ax_list = _get_axes_list(axes, do_clustered_pop, row)

            targets = {
                'target_column_all' : "colocalisation_fraction_{0}_{1}{2}".format(obj1_key, obj2_key,obj2_sufix),
                'target_column_free' : "colocalisation_fraction_{0}_free_{1}{2}".format(obj1_key, obj2_key,obj2_sufix),
                'target_column_clustered' : "colocalisation_fraction_{0}_clustered_{1}{2}".format(obj1_key, obj2_key,obj2_sufix),
                'target_column_all_sym' : "colocalisation_fraction_{0}{1}_{2}".format(obj2_key, obj2_sufix, obj1_key),
                'target_column_free_sym' : "colocalisation_fraction_{0}{1}_{2}_free".format(obj2_key, obj2_sufix, obj1_key),
                'target_column_clustered_sym' : "colocalisation_fraction_{0}{1}_{2}_clustered".format(obj2_key, obj2_sufix, obj1_key),
            }
            
            for obj2, ax in zip(obj2_list, ax_list) :
                
                title = "{0} --> {1} {2}".format(obj1,obj2, title_suffix[obj2_sufix])
                has_symetric = bool((obj1,obj2) in sym_couple_list)
                if has_symetric : 
                    title += "\n {0} --> {1}".format(obj1_key, obj2_key)
                    symetric_plot = True

                ax_data = sub_data.xs(obj2, axis=0, level=1, drop_level=False)
                if has_symetric : 
                    ax_data_sym = fraction_frame.xs(key=obj2,axis=0,level=0, drop_level=False).xs(key=obj1,axis=0,level=1, drop_level=False)
                    ax_data, ax_data_sym = _harmonize_sym_data(ax_data, ax_data_sym)

                else :
                    ax_data_sym = None



                labels = ax_data.index.get_level_values(-1)
                xticks_pos = np.arange(len(labels))
                x_range = [xticks_pos[0] - 1, xticks_pos[-1] + 1]

                ax.set_xlabel(ax_data.index.names[-1])
                ax.set_ylabel("colocalization fraction")
                ax.set_title(title)
                ax.set_xticks(xticks_pos,labels)
                ax.set_xlim(left=x_range[0], right=x_range[1])

                _plot_bars(
                    ax=ax,
                    ax_data=ax_data,
                    targets=targets,
                    obj1_key=obj1_key,
                    obj2_key=obj2_key,
                    do_clustered_pop=do_clustered_pop,
                    symetric_plot = symetric_plot,
                    has_sym=has_symetric,
                    ax_data_sym = ax_data_sym,

                ) 
                if has_symetric : ax.add_artist(
                    AnchoredText("{0} --> {1} {2}\n{3} --> {4}".format(obj1,obj2, title_suffix[obj2_sufix], obj2_key, obj1_key), loc=3)
                    )
                ax.legend(loc=1)

                #Verbose : show which index and which columns are targeted for each plot to check coherence with graph legend.
                if verbose : ax.add_artist(
                    AnchoredText("{0}\n{1}".format(ax_data.index, targets), loc=2)
                )
        fig.subplots_adjust(hspace=0.4)
        res_dict[obj1] = fig
        plt.close()

    return res_dict

def add_pvalue(
        fig_dict : dict,
        pvalues : pd.DataFrame,
        do_clustered_pop : bool,
        significance = 0.01

) :
    
    obj1_key,obj2_key = pvalues.index.names

    for obj1 in pvalues.index.get_level_values(0).unique() :
        fig:plt.Figure = fig_dict[obj1]
        fig = plt.figure(fig)

        obj2_idx = pvalues.loc[obj1].index.get_level_values(0).unique()
        ax_list = fig.axes
        col_num_shift = len(obj2_idx)

        for col_idx, obj2 in enumerate(obj2_idx) :
            ax_all = ax_list[col_idx]
            axes: 'list[plt.Axes]' = [ax_all]
            
            pvalue_all = pvalues.loc[(obj1,obj2)].at["{0}_{1}_colocalisation_count_pvalue".format(obj1_key, obj2_key)]
            pvalues_list = [pvalue_all]
            
            
            if do_clustered_pop :
                ax_clustered = ax_list[col_idx + col_num_shift]
                ax_free = ax_list[col_idx + 2*col_num_shift]
                axes +=  [ax_clustered, ax_free]

                pvalue_clustered = pvalues.loc[(obj1,obj2)].at["{0}_{1}_clustered_colocalisation_count_pvalue".format(obj1_key, obj2_key)]
                pvalue_free = pvalues.loc[(obj1,obj2)].at["{0}_{1}_free_colocalisation_count_pvalue".format(obj1_key, obj2_key)]
                pvalues_list += [pvalue_clustered, pvalue_free]

            for ax, pvalue in zip(axes, pvalues_list) :
                if pvalue <= significance :
                    ax = add_pvalue_star(ax)
        fig_dict[obj1] = fig
        plt.close()

    return fig_dict

# Output
def save_fig_dict(
        fig_dict : dict,
        obj_name : str,
        output_folder
        
        ) :
    #create folder
    plot_path = output_folder + "/colocalization/{0}/".format(obj_name)
    os.makedirs(plot_path, exist_ok=True)


    for obj in fig_dict.keys() :
        fig = fig_dict[obj]
        fig = plt.figure(fig)
        fig.savefig(plot_path + "/{0}_colocalization.png".format(obj))
        plt.close()
