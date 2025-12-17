import os
import pandas as pd
import numpy as np
import pbwrap.plot.box as box
import pbwrap.plot.violin as violin
import pbwrap.plot.bar as bar
import pbwrap.plot.stats as stats
import matplotlib.pyplot as plt
import pbwrap.quantification.statistical_test as s_test

from pbwrap.plot.utils import plot_horizontal_bar
from pbwrap.plot.utils import get_colors_list, make_color_frame

def open_data_and_merge(table_path, Acquisition_columns=None, Cell_columns= None, rna_list = None) :

    if type(Acquisition_columns) != type(None) :
        if not 'id' in Acquisition_columns : Acquisition_columns += ['id']

    if type(Cell_columns) != type(None) :
        for key in ['id', 'AcquisitionId'] :
            if not key in Cell_columns : Cell_columns += [key]

    Acquisition = pd.read_feather(table_path + 'Acquisition.feather', columns=Acquisition_columns)
    Acquisition = Acquisition.fillna('None')
    if type(rna_list) == list :
        idx = Acquisition.query("rna in {0}".format(rna_list)).index
        Acquisition = Acquisition.loc[idx,:]
    elif type(rna_list) == str :
        idx = Acquisition.query("rna == '{0}'".format(rna_list)).index
        Acquisition = Cell.loc[idx,:]

    for key in Acquisition.columns : 
        if key in Cell_columns and key != 'id' : Cell_columns.remove(key)
    Cell = pd.read_feather(table_path + 'Cell.feather', columns= Cell_columns)
    idx = (Acquisition.groupby('id')['fov'].count()[Acquisition.groupby('id')['fov'].count() > 1]).index
    
    #Rename specific treatment for order on plots
    Acquisition['treatment'] = Acquisition['treatment'].replace("dTAG8h", "dTAG08h")


    res = pd.merge(Acquisition, Cell, how= 'inner', left_on='id', right_on='AcquisitionId', validate= '1:m').drop('id_x', axis= 1).rename(columns={'id_y' : 'id'})
    res.set_index('id', verify_integrity=True)

    return res

def prepare_data_frame(table_path, data_columns, index_keys, rna_list = None) :
    """
    Open Acquisition and Cell dataframes located at table_path and prepare a merged dataframe for analysis.
    """
    
    base_ids = ['id', 'AcquisitionId']
    if 'id' in index_keys : base_ids.remove('id')
    if 'AcquisitionId' in index_keys : base_ids.remove('AcquisitionId')
    columns = base_ids + index_keys + data_columns
    Cell = open_data_and_merge(table_path, Cell_columns= columns, rna_list=rna_list)

    return Cell
    
def Cell_extract(table_path, output_path,  extract_columns, index_keys = ['rna', 'treatment', 'fov', 'id'], rna_list= None) :
    """
    
    """

    Cell = prepare_data_frame(table_path, extract_columns, index_keys, rna_list=rna_list)
    Cell = pd.concat([Cell.loc[:, ['rna', 'treatment', 'fov', 'id']], Cell.drop(['rna', 'treatment', 'fov', 'id'], axis= 1)], axis= 1) #Re-order columns
    Cell.sort_values(index_keys).to_excel(output_path + 'Cell_data_extract.xlsx')

def individual_box_plots(table_path, output_path, box_plots_columns, rna_list = None, index_keys = ['rna', 'treatment', 'fov', 'id'], folder_name = 'individual_box_plots/') :
    """
    
    """
    
    if len(index_keys) < 3 : raise ValueError("index_keys must contain at least 3 elements : 1 plot per [0], 1 box per [1] and data_grouping from [2:]")
    Cell = prepare_data_frame(table_path, box_plots_columns, index_keys, rna_list=rna_list)
    Cell = Cell.set_index(index_keys, verify_integrity=True)
    output_path += folder_name

    for column in box_plots_columns :
        print(' ',column)
        data = Cell.groupby(index_keys[:2], axis= 0)[column].apply(list)
        os.makedirs(output_path + column,exist_ok= True)
        for rna in data.index.get_level_values(0).unique() :
            data_rna = data.loc[rna,:]
            cell_number = [len(res)for res in data_rna]
            labels = [treatment + '\n {0} cells'.format(cell_num) for treatment, cell_num in zip(list(data_rna.index), cell_number)]
            box.box_plot(data_rna,show=False, ylabel= column, title= rna, path_output= output_path + column + '/{1}_{0}'.format(rna, column), labels = labels, showfliers= True)

def individual_violin_plots(table_path, output_path, measures, rna_list = None, folder_name = 'individual_violin_plots/', frameon=True) :
    """
    
    """
    Cell = prepare_data_frame(table_path, measures, [], rna_list=rna_list)
    color_df = make_color_frame(Cell['treatment'].unique())
    if not folder_name.endswith('/') : folder_name += '/'
    output_path += folder_name
    


    for column in measures :
        print(' ',column)
        data = Cell.dropna(subset= column).groupby(['rna', 'treatment'])[column].apply(list)
        os.makedirs(output_path + column,exist_ok= True)
        for rna in data.index.get_level_values(0).unique() :
            data_rna = data.loc[rna,:]
            cell_number = [len(res)for res in data_rna]
            labels = [treatment + '\n {0} cells'.format(cell_num) for treatment, cell_num in zip(list(data_rna.index), cell_number)]
            colors = pd.merge(data_rna, color_df, left_on= 'treatment', right_on='labels')['colors']

            FIGSIZE = 15
            fig = plt.figure(figsize= (FIGSIZE, FIGSIZE),frameon=frameon)
            ax = plt.gca()
            violin.violin_plot(ax, data_rna, labels= labels, colors=colors, ylabel= column, title= "{0}\n{1}".format(rna, column))

            plt.savefig(output_path + "{0}/{1}_{0}".format(column, rna))
            plt.close()

def individual_violin_rna_number_per_cell(table_path, output_path, rna_list=None, frameon= True) :

    #plot settings
    FIGSIZE = 15

    Cell = prepare_data_frame(table_path, ['nb_rna_out_nuc','nb_rna_in_nuc'], index_keys=[], rna_list=rna_list)
    Cell['rna_number'] = Cell['nb_rna_out_nuc'] + Cell['nb_rna_in_nuc']

    rna_grouping = Cell.groupby(['rna', 'treatment'], axis= 0)
    output_path += 'rna_number_per_cell/individual_distributions/'
    os.makedirs(output_path, exist_ok= True)

    data = rna_grouping['rna_number'].apply(list)

    #color_frame
    color_frame = make_color_frame(data.index.get_level_values(1).unique())

    rna_list = data.index.get_level_values(0).unique()
    for rna in rna_list :
        sub_data = data.loc[rna]
        colors = pd.merge(sub_data, color_frame, left_on= 'treatment', right_on= 'labels')['colors']

        fig = plt.figure(figsize=(FIGSIZE, FIGSIZE), frameon= frameon)
        ax = fig.gca()
        violin.violin_plot(ax, distributions=sub_data, labels= sub_data.index, colors=colors, y_axis= (-0, None), showextrema=True, showmedians=True, ylabel='rna per cell', title= ("{0}\nrna per cell distribution".format(rna)))

        path = output_path + "{0}_rna_per_cell_distribution".format(rna)
        plt.savefig(path)
        plt.close()

def rna_number_per_cell_plot(table_path, output_path, rna_list=None, frameon=True) :

    Cell = prepare_data_frame(table_path, ['nb_rna_out_nuc','nb_rna_in_nuc'], index_keys=[], rna_list=rna_list)
    Cell['rna_number'] = Cell['nb_rna_out_nuc'] + Cell['nb_rna_in_nuc']

    #Bar plot code --> To turn into overall_bar_plot ??
    rna_grouping = Cell.groupby(['rna', 'treatment'], axis= 0)
    output_path += 'rna_number_per_cell/'
    os.makedirs(output_path, exist_ok= True)

    data = rna_grouping['rna_number'].mean()
    errors = rna_grouping['rna_number'].std()

    rna_list = data.index.get_level_values(0).unique()
    treatment_list = data.index.get_level_values(1).unique()
    color_df = make_color_frame(treatment_list)
    colors = pd.merge(data.reset_index(), color_df, left_on= 'treatment', right_on='labels').sort_values(['rna', 'treatment'])['colors']

    data_list = data.rename("rna_per_cell").reset_index().groupby('rna')["rna_per_cell"].apply(list)
    errors_list = errors.rename("error").reset_index().groupby('rna')["error"].apply(list)

    #plot settings
    fig = plt.figure(figsize= (len(data), 10), frameon=frameon)
    ax = plt.gca()

    ax = bar.bar_plot(ax,
                      data= data_list,
                      errors= errors_list,
                      labels= rna_list,
                      colors=colors,
                      xlabel= 'rna',
                      ylabel= 'rna per cell',
                      title= "Number of rna per cell",
                      multi_bar_plot= True
                      )
    

    handles = bar._make_bar_handles(color_df['colors'])
    plt.legend(handles, color_df.index)

    plt.savefig(output_path + "Number_of_rna_per_cell")
    plt.close()

def overall_plot(table_path, output_path, measures, rna_list=None, folder_name = 'overall_plots', show_cell_numbers=True, frameon=True, **kargs) :
    """
    Overal violin plots
    """

    Cell = prepare_data_frame(table_path,measures, index_keys=[], rna_list=rna_list)

    for measure in measures :
        print(measure) 
        rna_grouping = Cell.replace(np.inf, np.NaN).dropna(subset=measure).groupby(['rna', 'treatment'], axis= 0)
        if not folder_name.endswith('/') : folder_name += '/'
        path = output_path + '{0}'.format(folder_name)
        os.makedirs(path, exist_ok= True)

        cell_number = rna_grouping[measure].apply(list).apply(len)
        _stats = rna_grouping[measure].apply(list).reset_index().groupby('rna')[measure].apply(list)
        data = rna_grouping[measure].apply(list).rename("data").reset_index().groupby('rna')["data"].apply(list) # making list of list...


        #plot settings
        treatment_list = cell_number.index.get_level_values(1).unique()
        rna_list = cell_number.index.get_level_values(0).unique()
        if show_cell_numbers :
            cell_number_list = [list(cell_number.loc[rna]) for rna in rna_list]
            labels = ["{0}\n{1}".format(rna, cell_numbers) for rna, cell_numbers in zip(rna_list, cell_number_list)]
        else : labels = [str(rna) for rna in rna_list]

        figsize = (
            len(cell_number_list) * 2 if len(cell_number_list) > 5 else 10,
            10
        )
        fig = plt.figure(figsize=figsize, frameon=frameon)
        ax = plt.gca()

        #Preparing colors
        color_df = make_color_frame(labels= treatment_list)
        colors = pd.merge(cell_number.reset_index(), color_df, left_on= 'treatment', right_on='labels').sort_values(['rna', 'treatment'])['colors']

        violin.violin_plot(ax,
                           data,
                           labels= labels,
                           colors= colors,
                           xlabel= 'rna [cell numbers]',
                           ylabel= measure,
                           title= '{0} general distribution'.format(measure),
                           multi_violin_plot= True,
                           linewith= 0.8,
                           showmeans=True,
                           **kargs
                           )
        

        #Horizontal lines
        if 'index' in measure :
            plot_horizontal_bar(y_value= 1)

        plt.legend(treatment_list)
        plt.xticks(rotation= 5)
        plt.savefig(path + measure)
        plt.close()

def cell_distribution(table_path, output_path, measure_list, rna_list=None, frameon=True) :
    Cell = prepare_data_frame(table_path, measure_list, index_keys= ['id'], rna_list=rna_list)
    Cell = Cell.set_index(['rna', 'treatment', 'id'], verify_integrity= True)
    
    #Preparing colors
    treatment_index = Cell.index.get_level_values(1).unique()
    color_df = pd.DataFrame({
        'treatment' : treatment_index,
        'color' : get_colors_list(len(treatment_index))
    }).set_index('treatment')

    for measure in measure_list :
        print("  {0}".format(measure))
        path = output_path + 'cell_distribution_plots/{0}/'.format(measure)
        os.makedirs(path, exist_ok= True)
        
        for rna in Cell.index.get_level_values(0).unique() :
            data_grouper = Cell.loc[rna,:].groupby('treatment', axis=0)[measure].apply(list)
            treatment_index = data_grouper.index

            #plot settings
            SUBPLOT_SIZE = 10
            subplots_number = len(treatment_index)
            fig = plt.figure(figsize= (SUBPLOT_SIZE * subplots_number, SUBPLOT_SIZE + 10), frameon=frameon)
            ax: 'list[plt.Axes]' = fig.subplots(1, subplots_number)
            if type(ax) != np.ndarray : ax = np.array([ax])
            yaxis = []
            
            #plotting histogram
            for subplot, treatment in zip(ax, treatment_index) :
                data = data_grouper.loc[treatment]
                cell_number = len(data)
                color = color_df.at[treatment, 'color']
                
                subplot.hist(data, color=color)
                subplot.set_title("{0}\nCell number : {1}".format(treatment, cell_number))
                subplot.set_ylabel("count")
                yaxis.append(subplot.axis()[3])
            
            #Uniformising scale
            ymax = max(yaxis)
            for subplot, treatment in zip(ax, treatment_index) :
                axis = list(subplot.axis())
                axis[2] = 0
                axis[3] = ymax
                subplot.axis(axis)

            plt.suptitle("{0}\n{1} distribution".format(rna, measure), fontsize= 'xx-large')
            plt.savefig(path + "{0}_cell_distribution".format(rna))
            plt.close()

def p_value(table_path, output_path, measure_list, significance= 0.05, folder_name='p_value', rna_list= None, show_cell_numbers= True, frameon=True, h0_treatment= 'DMSO') :
    Cell = prepare_data_frame(table_path,measure_list, index_keys=[], rna_list=rna_list)

    for measure in measure_list :
        print(measure) 
        rna_grouping = Cell.replace(np.inf, np.NaN).dropna(subset=measure).groupby(['rna', 'treatment'], axis= 0)
        if not folder_name.endswith('/') : folder_name += '/'
        path = output_path + '{0}'.format(folder_name)
        os.makedirs(path, exist_ok= True)

        cell_number = rna_grouping[measure].apply(list).apply(len)
        data = rna_grouping[measure].apply(list).rename("data")
        enlisted_data = data.reset_index().groupby('rna')["data"].apply(list) # making list of list...


        #plot settings
        treatment_list = cell_number.index.get_level_values(1).unique()
        rna_list = cell_number.index.get_level_values(0).unique()
        if show_cell_numbers :
            cell_number_list = [list(cell_number.loc[rna]) for rna in rna_list]
            labels = ["{0}\n{1}".format(rna, cell_numbers) for rna, cell_numbers in zip(rna_list, cell_number_list)]
        else : labels = [str(rna) for rna in rna_list]

        #Preparing colors
        color_df = make_color_frame(labels= treatment_list)
        colors = pd.merge(cell_number.reset_index(), color_df, left_on= 'treatment', right_on='labels').sort_values(['rna', 'treatment'])['colors']

        figsize = (len(cell_number_list) * 2, 20) if len(cell_number_list) > 5 else (10,20)

        fig = plt.figure(figsize= figsize, frameon=frameon)
        
        #measure distribution plot
        ax = plt.subplot(2,1,1)

        ax = violin.violin_plot(ax,
                           enlisted_data,
                           labels= labels,
                           colors= colors,
                           xlabel= 'rna [cell numbers]',
                           ylabel= measure,
                           title= '{0} general distribution'.format(measure),
                           multi_violin_plot= True,
                           linewith= 0.8,
                           showmeans=True,
                           )
        
        #Horizontal lines
        if 'index' in measure :
            plot_horizontal_bar(y_value= 1)
        legend = plt.legend(list(treatment_list))
        icons, texts = legend.get_patches(), legend.get_texts()
        texts = [text.get_text() for text in texts]
        plt.legend(icons + [plt.scatter(0,0, s= 20, c='white', linewidths= 1, edgecolors= 'black')], texts + ['distribution mean'])
        
        #p-value plot ANOVA
        ax = plt.subplot(2,2,3)
        title = "ANOVA test ({0} significance)".format(significance)
        ax = stats.ANOVA_test_plot(ax, data=data, significance=significance, xlabel=title, title=None, group= True)
        plt.xticks(rotation=45, ha='right')
        
        #p-value plot Tukey-hsd
        ax = plt.subplot(2,2,4)
        title = "Tukey_hsd test ({0} significance)\n H0 : {1}".format(significance, h0_treatment)
        ax = stats.Tukey_hsd_plot(ax, data=data, h0_treatment= h0_treatment, significance=significance, colors=colors, xlabel=title, title=None)

        plt.xticks(rotation=45, ha='right')
        plt.savefig(path + '{0}_p_value_superplot'.format(measure))
        plt.close()

def Tukey_hsd_tiles(table_path, output_path, measure_list, rna_list, significance= 0.01, folder_name= 'Tukey_tiles', **kwargs) :
    if not folder_name.endswith('/') : folder_name += '/'
    path = output_path + '{0}'.format(folder_name)
    os.makedirs(path, exist_ok= True)
    
    Cell = prepare_data_frame(table_path, measure_list, index_keys=[], rna_list=rna_list)
    
    for measure in measure_list :
        print(" " + measure)
        data = Cell.loc[:,["id","treatment","rna"] + [measure]]

        for rna in rna_list :
            print("     " + rna)
            if len(data[data['rna'] == rna]['treatment'].unique()) == 1: continue
            rna_data = data[data['rna'] == rna].dropna(subset= measure).groupby("treatment")[measure].apply(list)

            #plot
            figsize = (
                len(rna_data) * 2 if len(rna_data) > 10 else 20,
                20
            )
            plt.figure(figsize=figsize,
                       frameon= kwargs.get('frameon'))
            
            #distribution
            treatment_list = list(rna_data.index)
            cell_number_list = list(rna_data.apply(len))
            labels = ["{0} \n[{1}]".format(treatment, cell_number) for treatment, cell_number in zip(treatment_list, cell_number_list)]

            ax = plt.subplot(2,1,1)
            ax = violin.violin_plot(ax,
                               rna_data,
                               labels= labels,
                               xlabel= 'rna [cell numbers]',
                               ylabel= measure,
                               title= '{0} [{1}]'.format(measure, rna),
                               multi_violin_plot= False,
                               linewith= 0.8,
                               showmeans=True,
                               )

            #Alexandergovern
            ax = plt.subplot2grid((2,6), (1,5), colspan=1)
            title = "AlexanderGovern test ({0} significance)".format(significance)
            ax = stats.alexandergovern_test_plot(ax, data=rna_data, significance=significance, xlabel=title, title=None, xtick_label= [rna])
            plt.xticks(rotation=45, ha='right')            

            #Gameshowell tile
            ax = plt.subplot2grid((10,6), (6,0), colspan=4, rowspan=4)
            stats.pairwise_stest_tile(
                ax,
                rna_data,
                measure=None,
                groupby_key=None,
                title= None,
                xlabel= "Games-Howell pairwise p-values",
                significance=significance,
                test= 'gameshowell'
                )

            plt.savefig(path + '{0}_{1}_Tukey_Tile'.format(measure, rna))
            plt.close()