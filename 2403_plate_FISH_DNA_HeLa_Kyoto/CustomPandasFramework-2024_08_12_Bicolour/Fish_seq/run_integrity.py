"""
Module to check data integrity for Bicolor/colocalisation/seq-fish pipelines
"""

import pandas as pd
import numpy as np
from CustomPandasFramework.computer_interface import get_datetime
from ..integrity import get_referencement_relation
from itertools import combinations
    

def _spots_summarize(Spots:pd.DataFrame) :
    spots_sum = Spots.groupby(["spot_type", "cell_id"])['spots_id'].count().rename('count')
    
    return spots_sum

def _cell_spots_summarize(Cell:pd.DataFrame) :
    summ = Cell.loc[:, #Reducing col number
                   ['cell_id'] + ["{0}_nb_rna_in_nuc".format(obj) for obj in ['rna1','rna2','suntag']] + ["{0}_nb_rna_out_nuc".format(obj) for obj in ['rna1','rna2','suntag']]
                   ].fillna(0)
    
    for obj in ['rna1','rna2','suntag'] : # Getting full count of rnas
        summ['{0}'.format(obj)] = summ["{0}_nb_rna_in_nuc".format(obj)] + summ["{0}_nb_rna_out_nuc".format(obj)]

    rna_total_columns = ["{0}".format(obj) for obj in ['rna1','rna2','suntag']]
    summ = summ.loc[:,['cell_id'] + rna_total_columns] 

    summ = summ.melt(value_vars= rna_total_columns, id_vars='cell_id', value_name='count',var_name='spot_type') #Passing columns as index
    summ = summ.set_index(['spot_type', 'cell_id'], verify_integrity=True)

    return summ.sort_index()

def _colloc_spots_summarize(Cell:pd.DataFrame) :

    name_dict = {
        'rna1_number' : 'rna1',
        'rna2_number' : 'rna2',
        'suntag_number' : 'suntag',
    }

    summ = Cell.loc[:,['cell_id','rna1_number', 'rna2_number', 'suntag_number']].rename(columns=name_dict).fillna(0)
    summ = summ.melt(id_vars='cell_id', value_vars=['rna1', 'rna2', 'suntag'], value_name= 'count', var_name= 'spot_type')
    summ = summ.set_index(['spot_type', 'cell_id'], verify_integrity=True)

    return summ.sort_index()

def _get_colloc_count_columns(spot_types, populations) :
    columns =[]

    for spot_type1  in spot_types:
        for spot_type2 in spot_types :
            if spot_type2 == spot_type1 : continue
            ref_columns = [
                '{0}_{1}_colocalisation_count',
                '{0}_clustered_{1}_colocalisation_count',
                '{0}_free_{1}_colocalisation_count',
                ]

            for population in populations :
                    new_columns = [column.replace('{0}' ,spot_type1 ).replace('{1}', spot_type2 + population) for column in ref_columns]
                    columns.extend(new_columns)
        columns.append('{0}_number'.format(spot_type1))
        columns.append('{0}_clustered_number'.format(spot_type1))
        columns.append('{0}_free_number'.format(spot_type1))

    return columns

def _get_df_columns(spot_types, populations) :

    columns = []

    for spot_type in spot_types :
        for population in populations :
            columns.append(spot_type + population)
    
    return columns

def _colloc_count_consistency(Cell:pd.DataFrame) :
    spot_types = ['rna1', 'rna2', 'suntag']
    populations = ['', '_clustered', '_free']

    df_columns =  _get_df_columns(spot_types, populations)
    integrity_df = pd.DataFrame(columns=spot_types, index= df_columns, dtype=bool)
    columns = ['acquisition_id', 'cell_id'] + _get_colloc_count_columns(spot_types, populations)
    cell_view = Cell.loc[:,columns].reset_index(drop=True).infer_objects(copy=False)

    fail_index = []

    for spot_type1 in spot_types :
        for population in populations :
            for spot_type2 in spot_types :
                if spot_type1 == spot_type2 : continue
                pop_max = cell_view['{0}_number'.format(spot_type1)].fillna(0)
                total = cell_view['{0}_{1}_colocalisation_count'.format(spot_type1, spot_type2 + population)].fillna(0)
                clustered = cell_view['{0}_clustered_{1}_colocalisation_count'.format(spot_type1, spot_type2 + population)].fillna(0)
                free = cell_view['{0}_free_{1}_colocalisation_count'.format(spot_type1, spot_type2 + population)].fillna(0)

                counts_consistency = (total == clustered + free) & (total <= pop_max)
                
                integrity_df.at[spot_type2 + population, spot_type1] = all(counts_consistency)

                fail_index += list(cell_view[~counts_consistency].index)

    fail_view = cell_view.loc[fail_index]

    return integrity_df, fail_view

def _filename_duplicate_check(Acquisition: pd.DataFrame, filename_key = "filename") :
    count = Acquisition.value_counts(subset= filename_key)
    if any(count > 1) : print("Duplicate found in Acquisition : several acquisitions with same filename.")

    return count[count>1]


def full_integrity_check(Acquisition: pd.DataFrame, Cell: pd.DataFrame, Spots: pd.DataFrame) :
    
    Acquisition_is_empty = len(Acquisition) == 0
    Cell_is_empty = len(Cell) == 0
    Spots_is_empty = len(Spots) == 0
    
    Acquisition_Cell_relation = get_referencement_relation(Acquisition, 'acquisition_id', Cell, 'acquisition_id')
    Acquisition_Spots_relation = get_referencement_relation(Acquisition, 'acquisition_id', Spots, 'acquisition_id')
    Spots_Cell_relation = get_referencement_relation(Spots, 'cell_id', Cell, 'cell_id')

    sum_spots = _spots_summarize(Spots)
    sum_cell = _cell_spots_summarize(Cell)
    sum_coloc = _colloc_spots_summarize(Cell)

    spots_count_consistency_1 = sum_cell.eq(sum_coloc).all().iat[0]
    spots_count_consistency_2 = (len(sum_spots) == len(pd.merge(left= sum_spots.reset_index(), right=sum_cell.reset_index(), on=['spot_type', 'cell_id', 'count'], how='inner')))
    spots_count_consistency = spots_count_consistency_1 and spots_count_consistency_2

    log = {
        "Acquisition_is_empty" : Acquisition_is_empty,
        "Cell_is_empty" : Cell_is_empty,
        "Spots_is_empty" : Spots_is_empty,
        "Acquisition_Cell_relation" : Acquisition_Cell_relation,
        "Acquisition_Spots_relation" : Acquisition_Spots_relation,
        "Spots_Cell_relation" : Spots_Cell_relation,
        "spots_count_consistency" : spots_count_consistency # check total number of spot is consistent
    }

    return log


def write_integrity_log(path: str, log : dict) :
    if not path.endswith('/') : path += '/'
    time = get_datetime()

    with open(path + "data_integrity_log.txt", 'w+') as txt : 
        txt.write("This log was created on {0}\n\n".format(time))
        
        for test, res in log.items() :
            txt.write("{0} : {1}\n".format(test,res))
    
    return True