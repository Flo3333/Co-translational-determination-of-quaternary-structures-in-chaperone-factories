"""
Data pre-processing script for Soha quantifications
Pipeline : Fish_seq / BC_clusterkinetics_pipeline


This analysis restricts cells (counted as phenotype_positive) to cell showing numerous clusters from APC channel (rna1).
This script updates Cell table by deleting cells that were manually deleted in the positive phenotype visual.

"""

import pandas as pd
import os,re


merged_data_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/Soha quantification/Merged_results/"

REVERSE = False


#############
if not REVERSE :
    if not os.path.isfile(merged_data_path + '/Cell_beforeupdate.feather') :
        Cell = pd.read_feather(merged_data_path + '/Cell.feather')
        Cell.to_feather(merged_data_path + '/Cell_beforeupdate.feather')
    else :
        Cell = pd.read_feather(merged_data_path + '/Cell_beforeupdate.feather')

    datasets = os.listdir(merged_data_path + 'positive_phenotype_visuals/')

    remaning_files = []
    for dataset in datasets :
        if os.path.isdir(merged_data_path + 'positive_phenotype_visuals/' + dataset) : 
            remaning_files += (os.listdir(merged_data_path + 'positive_phenotype_visuals/' + dataset))

    df = pd.DataFrame({'cell_id' : remaning_files})

    def _get_cell_id(filename) :
        regex = r".*_(\d+).png"
        try : 
            cell_id = re.findall(regex, filename)[0]
        
        except Exception as e:
            print(filename)

        return float(cell_id)


    id_list = list(df['cell_id'].apply(_get_cell_id))

    drop_idx = Cell.query('treatment == "untreated" and cell_id not in {0}'.format(id_list)).index
    Cell = Cell.drop(drop_idx, axis=0)
    if 'rna2' in Cell : Cell = Cell.drop(columns= 'rna2')
    print("{0} cells were deleted.".format(len(drop_idx)))
    if 'level_0' in Cell.columns : Cell = Cell.drop(columns= 'level_0', axis=1)
    Cell.reset_index(drop=False).to_feather(merged_data_path + '/Cell.feather')

    Acquisition = pd.read_feather(merged_data_path + '/Acquisition.feather', columns= ['acquisition_id', 'rna1', 'rna2'])
    Acquisition['acquisition_id'] = Acquisition['acquisition_id'].apply(tuple)
    Cell['acquisition_id'] = Cell['acquisition_id'].apply(tuple)
    Cell = pd.merge(Acquisition, Cell, on= 'acquisition_id')
    Cell.value_counts(subset=['rna1', 'rna2', 'treatment']).rename('count').to_excel(merged_data_path + '/cell_summary_after_cleaning.xlsx')

else :
    Cell = pd.read_feather(merged_data_path + '/Cell_beforeupdate.feather')
    Cell.to_feather(merged_data_path + '/Cell.feather')