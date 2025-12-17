import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pbwrap.plot as plot
import pbwrap.plot.bar as bar
import pbwrap.plot.violin as violin
import pbwrap.plot.scatter as scatter
import pbwrap.plot.stats as stats_plot
import pbwrap.quantification.statistical_test as stats_test
from pbwrap.plot.utils import make_color_frame

def prepare_cololocalisation_data(tables_path: str, gene_list= None) : 

    if not tables_path.endswith('/') : tables_path += '/'
    Acquisition_columns = ['id','rna1','rna2','treatment','fov_number']
    Colocalisation_columns = None

    Acquisition = pd.read_feather(tables_path + 'Acquisition.feather', columns= Acquisition_columns)
    Colocalisation = pd.read_feather(tables_path + 'Colocalisation.feather', columns= Colocalisation_columns)

    if type(gene_list) != type(None) :
        keep_idx = Acquisition.query('rna1 in {0} or rna2 in {0}'.format(gene_list)).index
        Acquisition = Acquisition.loc[keep_idx,:]

    coloc_df = pd.merge(Acquisition, Colocalisation, left_on='id', right_on='AcquisitionId').drop(['id_x', 'id_y'], axis= 1)
    coloc_df = coloc_df.sort_values(['rna1','rna2','treatment','fov_number']).reset_index(drop=True)
    
    assert len(Acquisition) == len(coloc_df)

    return coloc_df    

def compute_colocalisation_fractions(coloc_df: pd.DataFrame, assertion_test=True) :

    rna1_tot = coloc_df['rna1_number']
    rna2_tot = coloc_df['rna2_number']
    fraction_keys = ['rna1_colocalising_with_rna2_fraction', 'rna1_colocalising_with_suntag_fraction', 'rna2_colocalising_with_rna1_fraction', 'rna2_colocalising_with_rna1_fraction', 'rna2_colocalising_with_suntag_fraction', 'rna1_closer_1000nm_rna2_fraction', 'rna1_closer_2000nm_rna2_fraction']
    coloc_df['rna1_colocalising_with_rna2_fraction'] = coloc_df['rna1_rna2_colocalisation_count'] / rna1_tot
    coloc_df['rna1_colocalising_with_suntag_fraction'] = coloc_df['rna1_suntag_colocalisation_count'] / rna1_tot
    coloc_df['rna2_colocalising_with_rna1_fraction'] = coloc_df['rna2_rna1_colocalisation_count'] / rna2_tot
    coloc_df['rna2_colocalising_with_suntag_fraction'] = coloc_df['rna2_suntag_colocalisation_count'] / rna2_tot
    coloc_df['rna1_closer_1000nm_rna2_fraction'] = coloc_df['rna1_number_1000nm_around_rna2'] / rna1_tot
    coloc_df['rna1_closer_2000nm_rna2_fraction'] = coloc_df['rna1_number_2000nm_around_rna2'] / rna1_tot
    
    if assertion_test :
        for key in fraction_keys : assert all(coloc_df[key] <= 1), "{0} values were found > 1".format(key)

    return coloc_df

def colocalisation_extract(tables_path: str, save_path:str, gene_list= None) :

    extract = prepare_cololocalisation_data(tables_path, gene_list)
    extract = compute_colocalisation_fractions(extract)

    if type(save_path) != type(None) : 
        if not save_path.endswith('/') : save_path += '/'
        os.makedirs(save_path + '/excel_extracts/', exist_ok= True)
        extract.to_excel(save_path + '/excel_extracts/' + 'colocalisation_extract.xlsx')
    
def prepare_clusteredspots_data(tables_path, gene_list=None) :
    
    if not tables_path.endswith('/') : tables_path += '/'

    Acquisition_columns = ['id','treatment','fov_number']
    ClusteredSpots_columns = None

    Acquisition = pd.read_feather(tables_path + 'Acquisition.feather', columns= Acquisition_columns)
    Cluster = pd.read_feather(tables_path + 'ClusteredSpots.feather', columns= ClusteredSpots_columns)
    
    #gene_list filter
    if type(gene_list) != type(None) :
        keep_idx = Cluster.query('rna in {0}}'.format(gene_list)).index
        Cluster = Cluster.loc[keep_idx,:]
    merge_df = pd.merge(Cluster, Acquisition, how='left', left_on='AcquisitionId', right_on='id').drop('id_y', axis=1).rename(columns={'id_x' : 'id'}) #This will filter acquisition with rna not in rna_list

    return merge_df

def clusteredspots_data_mean_quantif(clusteredspots_data, groupby_index = ['rna','treatment','fov_number']) :

    #Sorting free and clustered data
    free_idx = clusteredspots_data.query('cluster_id.isna()').index
    clustered_idx = clusteredspots_data.query('not cluster_id.isna()').index
    clusteredspots_data.loc[free_idx, 'cluster_id'] = -1

    free_grouper = clusteredspots_data.loc[free_idx,:].groupby(groupby_index)
    clustered_grouper = clusteredspots_data.loc[clustered_idx,:].groupby(groupby_index)

    #Quantifications
    free_df = pd.concat(
        [
        free_grouper['id'].count().rename('free_spot_number'),
        free_grouper['colocalising_suntag_number'].mean().rename('mean_suntag_per_free_spot'),
        free_grouper['colocalising_suntag_number'].std().rename('std_suntag_per_free_spot'),
        free_grouper['colocalising_suntag_number'].median().rename('median_suntag_per_free_spot'),
        free_grouper['colocalising_suntag_number'].sum().rename('total_suntag_colocalising_free_spot')
    ],
    axis=1
    )

    clustered_df = pd.concat(
        [
        clustered_grouper['id'].count().rename('clustered_spot_number'),
        clustered_grouper['colocalising_suntag_number'].mean().rename('mean_suntag_per_clustered_spot'),
        clustered_grouper['colocalising_suntag_number'].std().rename('std_suntag_per_clustered_spot'),
        clustered_grouper['colocalising_suntag_number'].median().rename('median_suntag_per_clustered_spot'),
        clustered_grouper['colocalising_suntag_number'].sum().rename('total_suntag_colocalising_clustered_spot')
    ],
    axis=1
    )

    df = pd.merge(free_df, clustered_df, how= 'inner', validate= '1:1', on= groupby_index)
    df = df.reset_index(drop=False).sort_values(groupby_index)
    
    if groupby_index == ['rna','treatment','fov_number'] : assert len(df) == len(free_df) == len(clustered_df)

    return df

def cluster_quantification_extract(tables_path, save_path, gene_list= None) :
    
    extract = prepare_clusteredspots_data(tables_path, gene_list)
    extract = clusteredspots_data_mean_quantif(extract)
    if not save_path.endswith('/') : save_path += '/'
    os.makedirs(save_path + '/excel_extracts/', exist_ok= True)
    extract.to_excel(save_path + '/excel_extracts/' + 'clusters_extract.xlsx')

def colocalisation_plots(tables_path, save_path, measure_list, gene_list=None, **kargs) :

    if not save_path.endswith('/') : save_path += '/'
    os.makedirs(save_path + 'colocalisation_plots/', exist_ok=True)

    #Data_load
    colocalisation_df = prepare_cololocalisation_data(tables_path, gene_list)

    grouper = colocalisation_df.groupby(['rna1','rna2', 'treatment','fov_number']).sum(numeric_only=True)
    data = compute_colocalisation_fractions(grouper)

    #colors
    treatment_list = colocalisation_df['treatment'].unique()
    color_df = make_color_frame(treatment_list)
    print(colocalisation_df.columns)

    for measure in measure_list :
        if measure not in colocalisation_df.columns :
            print('{0} measure was not found in colocalisation data. Computing next measurement.'.format(measure))
            continue
        else : print(measure)
        measure_data = data.reset_index().groupby(['rna1','rna2','treatment'])[measure].apply(list)
        if 'suntag' in measure  and 'rna1' in measure:
            labels = [rna1 for rna1 in measure_data.index.get_level_values(0).unique()]
            xlabel = 'rna1'
        elif 'suntag' in measure  and 'rna2' in measure : 
            labels = [rna2 for rna2 in measure_data.index.get_level_values(1).unique()]
            xlabel = 'rna2'
        else : 
            labels = ['{0}\n/\n{1}'.format(rna1,rna2) for rna1,rna2 in zip(measure_data.index.get_level_values(0).unique(),measure_data.index.get_level_values(1).unique())]
            xlabel = '(rna1 / rna2)'

        #Plot
        colors = pd.merge(measure_data, color_df, left_on= 'treatment', right_on= 'labels')['colors']
        enlisted_measure_data = measure_data.reset_index().groupby(['rna1','rna2'])[measure].apply(list)

        figsize = (
            len(enlisted_measure_data)*2 if len(enlisted_measure_data) > 5 else 10,
            20
        )

        plt.figure(
            figsize= kargs.setdefault("figsize", figsize),
            frameon= kargs.get("frameon")
            )
        
        #data plot
        ax = plt.subplot(2,1,1)
        ax = scatter.vertical_scatter(
            ax,
            distributions=enlisted_measure_data,
            labels= labels,
            colors= colors,
            ylabel= measure,
            title= None,
            showmeans= True,
            show_std= True,
            multi_distributions=True  
        )

        handles = bar._make_bar_handles(color_df['colors']) + [plt.scatter(0,0, s= 20, marker= 'D', c = 'white', linewidths= 1, edgecolors= 'black')]
        plt.legend(handles, list(color_df.index) + ['Sample mean'])

        #p_value plot
        ax = plt.subplot(2,1,2)
        ax = stats_plot.p_value_plot(
            ax,
            data= measure_data.apply(np.mean),
            statistical_test= stats_test.chi_squared_test,
            significance= 0.01,
            colors= colors,
            xlabel= 'Chi-squared test',
            ignore_level_test = True
            )

        plt.savefig(save_path + 'colocalisation_plots/' + '{0}_total_bar_plot'.format(measure))
        plt.close()

def compute_clustered_spots_fractions(clusteredspots_df : pd.DataFrame) :
    clusteredspots_df['clustered_spot_fraction'] = clusteredspots_df['clustered_spot_number'] / (clusteredspots_df['free_spot_number'] + clusteredspots_df['clustered_spot_number'])
    return clusteredspots_df

def compute_p_values() :
    pass

def clusters_plots(tables_path, save_path, gene_list=None, **kargs) :

    if not save_path.endswith('/') : save_path += '/'
    os.makedirs(save_path + 'clusteredspots_plots/', exist_ok=True)

    df = prepare_clusteredspots_data(tables_path, gene_list=gene_list)
    df = clusteredspots_data_mean_quantif(df, groupby_index= ['rna','treatment'])
    df = compute_clustered_spots_fractions(df)

    measure_list = ['free_spot_number', 'clustered_spot_number', 'mean_suntag_per_free_spot', 'mean_suntag_per_clustered_spot','clustered_spot_fraction']
    std_list = [None, None, 'std_suntag_per_free_spot', 'std_suntag_per_clustered_spot',None]
    assert len(measure_list) == len(std_list)

    #colors
    treatment_list = df['treatment'].unique()
    color_df = make_color_frame(treatment_list)
    
    for measure, std in zip(measure_list, std_list) :
        if measure not in df.columns :
            print('{0} measure was not found in colocalisation data. Computing next measurement.')
            continue
        else : print(measure)
        data = df.reset_index().groupby(['rna'])[measure].apply(list)
        error = df.reset_index().groupby(['rna'])[std].apply(list) if std != None else None

        #Plot
        colors = pd.merge(df.reset_index(), color_df, left_on= 'treatment', right_on= 'labels').sort_values(['rna','treatment'])['colors']

        plt.figure(
            figsize= kargs.setdefault("figsize", (len(data)*2, 10)),
            frameon= kargs.get("frameon")
            )
        ax = plt.gca()

        bar.bar_plot(
            ax,
            data= data,
            errors=error,
            colors=colors,
            labels= data.index,
            ylabel= str(measure),
            y_axis= (0,None),
            multi_bar_plot= True
        )

        handles = bar._make_bar_handles(color_df['colors'])
        plt.legend(handles, color_df.index)
        plt.savefig(save_path + 'clusteredspots_plots/' + '{0}_bar_plot'.format(measure))
        plt.close()

def closest_suntag_distance(tables_path, save_path, gene_list=None, **kargs) :

    if not save_path.endswith('/') : save_path += '/'
    os.makedirs(save_path + 'suntagdistance/', exist_ok=True)

    df = prepare_clusteredspots_data(tables_path, gene_list)
    df['normalised_distance_closest_suntag'] = df['distance_closest_suntag'] / df['median_suntag_distance']
    df = df.groupby(['rna', 'treatment'])['normalised_distance_closest_suntag'].apply(list)
    data = df.reset_index().groupby('rna')['normalised_distance_closest_suntag'].apply(list)
    
    #
    treatment_list = df.index.get_level_values(1).unique()

    #colors
    colors_df = make_color_frame(treatment_list)
    colors = pd.merge(df.reset_index(), colors_df, left_on= "treatment", right_on='labels').sort_values(['rna','treatment'])['colors']

    fig = plt.figure(
        figsize= kargs.setdefault("figsize", (len(df), 10)),
        frameon= kargs.get("frameon")
    )
    ax = plt.gca()

    violins = violin.violin_plot(
        ax,
        distributions=data,
        labels= data.index,
        colors= colors,
        y_axis= (0,3),
        ylabel= 'Normalised distance to closest suntag',
        title= 'Distance RNAs to suntag',
        linewith= 1,
        showextrema= False,
        multi_violin_plot=True
        ) 
    
    #Horizontal bar
    plot.plot_horizontal_bar(y_value= 1)

    plt.legend(treatment_list)
    plt.savefig(save_path + 'suntagdistance/' + 'normalised_distance_to_closest_suntag')
    plt.close()