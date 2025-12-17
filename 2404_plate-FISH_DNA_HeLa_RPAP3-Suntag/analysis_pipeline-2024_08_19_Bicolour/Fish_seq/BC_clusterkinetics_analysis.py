import os
import pandas as pd
import matplotlib.pyplot as plt
import pbwrap.plot.superplot as superplot

group = ['date', 'treatment'] # Heloise data
# group = ['cell_line', 'treatment'] # Manon data


def open_clusteredspots(input_path, conditions) :
    Spots = pd.read_feather(input_path + 'ClusteredSpots.feather')
    drop_idx = Spots.query('cell_label == 0').index
    Spots = Spots.drop(drop_idx, axis= 0)
    Spots = _add_experiment_conditions(input_path, Spots, conditions)

    return Spots

def _add_experiment_conditions(input_path, Spots: pd.DataFrame, conditions) :
    Acquisition = pd.read_feather(input_path + "Acquisition.feather")

    Spots = pd.merge(Acquisition.loc[:,conditions], 
             Spots, 
             left_on='id', 
             right_on='AcquisitionId', 
             how= 'inner', 
             validate= '1:m').drop('id_x', axis= 1).rename(columns={'id_y' : 'id'})
    
    return Spots

def _plot(data, xlabel, ylabel, output_path=None, fig_name='fig', axis=None) :
    
    if type(output_path) != type(None) :
        if not output_path.endswith('/') : output_path += '/'

    fig = plt.figure(figsize= (10,10))
    ax = fig.gca()
    if len(data) > 0 :
        ax = superplot.distribution_super_plot(
            data, 
            ax,
            xlabel= xlabel,
            ylabel= ylabel
            )
        
        #Auto axis or set axis or part set axis
        xmin, xmax, ymin, ymax = plt.axis()
        if type(axis) != type(None) :
            xmin_passed,xmax_passed,ymin_passed,ymax_passed = axis
            if xmin_passed == None : xmin_passed = xmin 
            if xmax_passed == None : xmax_passed = xmax 
            if ymin_passed == None : ymin_passed = ymin 
            if ymax_passed == None : ymax_passed = ymax
        else :  xmin_passed, xmax_passed, ymin_passed, ymax_passed = xmin, xmax, ymin, ymax
        plt.axis([xmin_passed, xmax_passed, ymin_passed, ymax_passed])

    # plt.xticks(rotation=60, )
    plt.tight_layout()

    if type(output_path) == type(None) :
        plt.show()
    else :
        plt.savefig(output_path + fig_name)
        plt.close(fig=fig)

def cluster_per_cell(Spots : pd.DataFrame, output_path= None) :
    data = Spots.groupby(group + ['fov_number', 'cell_label'])['cluster_id'].nunique(dropna=True).rename('cluster_number')
    data = data.reset_index().groupby(group)['cluster_number'].apply(list)

    _plot(data, 'Cell Line\n[cell number]', 'Cluster per cell number', output_path, fig_name= 'cluster_per_cell.png', axis= [None,None,0,None])

def nuclear_cluster(Spots: pd.DataFrame, output_path=None) :
    data = Spots.groupby(group + ['fov_number', 'cell_label', 'cluster_id'])['is_nuclear'].min().rename('is_nuclear')
    grouper = data.reset_index().groupby(group + ['fov_number', 'cell_label'])['is_nuclear']
    count_nuclear = grouper.sum().rename('nuclear_cluster_count')
    count_total = grouper.count().rename('cluster_count')
    data = (count_nuclear / count_total).rename("nuclear_cluster_fraction")
    data = data.reset_index().groupby(group)['nuclear_cluster_fraction'].apply(list)

    _plot(data, 'Cell Line\n[cell number]', 'Nuclear cluster fraction', output_path, fig_name= 'nuclear_cluster_per_cell.png', axis=[None,None,0,None])
    
def proportion_rna_in_cluster(Spots : pd.DataFrame, output_path=None):
    nb_in_cluster = Spots.groupby(group + ['fov_number', 'cell_label'])['cluster_id'].count()
    nb_total = Spots.groupby(group + ['fov_number', 'cell_label'])['id'].count()
    data = (nb_in_cluster / nb_total).rename("proportion_rna_cluster").reset_index().groupby(group)['proportion_rna_cluster'].apply(list)

    _plot(data, 'Cell Line\n[cell number]', 'rna propotion in cluster', output_path=output_path, fig_name= 'proportion_rna_in_clusters')

def rna_per_cell(Spots : pd.DataFrame, output_path = None) :
    data = Spots.groupby(group + ['fov_number', 'cell_label'])['id'].count()
    data = data.rename("rna_per_cell").reset_index().groupby(group)['rna_per_cell'].apply(list)

    _plot(data, 'Cell Line\n[cell number]', 'rna per cell', output_path=output_path, fig_name= 'number_of_rna_per_cell')

def spot_per_cluster(Spots, output_path=None):
    data = Spots.groupby(group + ['fov_number', 'cell_label', 'cluster_id'])['id'].count()
    data = data.rename("rna_per_cluster").reset_index().groupby(group)['rna_per_cluster'].apply(list)
    print(data)

    _plot(data, 'Cell Line\n[cluster number]', 'rna per cluster', output_path=output_path, fig_name= 'number_of_rna_per_cluster', axis=[None,None,None,None])



### MAIN SCRIPT ###

# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/Cinétique des foyers/240206_hct_h9_polr2a/pipeline_output/20240305_14-20-31/result_tables/"
# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/Cinétique des foyers/240206_hct_h9_polr2a/pipeline_output/20240306_15-23-11/result_tables/"
# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/20240306_16-14-06/result_tables/"
# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/20240307_12-17-29/result_tables/"
# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/20240318_11-13-02/result_tables/" # Manon dataset with nuclear cluster count
# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/20240321_10-48-30/result_tables/" # Manon dataset with nuclear cluster count
# input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/20240325_14-52-03/result_tables/" # Heloise dataset
input_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/pipeline_output/20240402_11-39-05/result_tables/" # Heloise dataset; cy3 very good


output_path = input_path.replace('result_tables', 'result_plots')
os.makedirs(output_path, exist_ok=True)

conditions = ['id', 'date', 'treatment', 'fov_number'] # Heloise
# conditions = ['id', 'cell_line', 'treatment', 'fov_number'] # Manon

Spots = open_clusteredspots(input_path=input_path, conditions=conditions)
print(Spots)
print(Spots.value_counts(subset=['date', 'rna', 'treatment']).sort_index())


for rna in Spots['rna'].unique() :
    path_out = output_path + "/{0}/".format(rna)
    os.makedirs(path_out, exist_ok=True)

    Spots_loc = Spots.loc[Spots['rna'] == rna]

    cluster_per_cell(Spots_loc, path_out)
    try :
        nuclear_cluster(Spots_loc, path_out)
    except Exception as e :
        print("nuclear cluster failed : \n{0}".format(e))
    proportion_rna_in_cluster(Spots_loc, output_path=path_out)
    rna_per_cell(Spots_loc, path_out)
    spot_per_cluster(Spots_loc,path_out)