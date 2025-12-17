"""
Data pre-processing script for BC_clusterkinetics.py outputs
Pipeline : Fish_seq / BC_clusterkinetics_pipeline

Merge run results together along with some data checks and name cleansing.

"""

import os
import functools
import pandas as pd
import matplotlib.pyplot as plt
from CustomPandasFramework.computer_interface import get_datetime
from CustomPandasFramework.Fish_seq.run_integrity import _colloc_count_consistency, _filename_duplicate_check

PATH_LIST = [
    "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/pipeline_output/2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto_20240626_11-41-05/result_tables/",
    "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/pipeline_output/2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto_20240628_10-58-42/result_tables/",
    "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/pipeline_output/2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto_20240628_13-38-25/result_tables/",
    "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/pipeline_output/2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto_20240701_10-52-09/result_tables/",
]

DROP_RNA = [
    [('SCAF4',''), ('JUP','R2A'), ('JUP','EP400'), ('JUP','EFTUD2'), ('R2B','R2A'), ('R2A','R2B'), ('JUP','')],
    [('JUP','EFTUD2'), ('JUP','R2A')],
    [],
    [],
]


PATH_OUT = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/R2TP/2403_plate_FISH_DNA_HeLa_Kyoto/merge_folder/"
DROP_NEGATIVE_CELLS = False
MIN_NUMBER_FOCI= 3
MIN_NUMBER_FOCI_ARGAP= 10
EXPECTED_TREATMENT_NUMBER = 2 # puro, unt
EXPECTED_RNA1 = ['R1A', 'R2A', 'R1A', 'JUP', 'BRD8', 'EP400', 'PRPF8', 'EFTUD2', 'R2B', 'R2A', 'SCAF4', 'GOLGA4', 'FLAP-only', 'SRCAP', 'CTL', 'A-CATENIN']
EXPECTED_RNA2 = ['R2A', 'R1A', 'BRD8', 'EP400', 'EFTUD2', 'PRPF8', 'R2B','R1A','EMPTY']
EXPECTED_SUNTAG = False

"""

1. Data Opening&Merging

"""

print("Merging data...")

Acquisition_columns = None
Cell_columns = None
Spots_columns = None

# ------

tables = {
'Acquisition' : pd.DataFrame(),
'Cell' : pd.DataFrame(),
'Spots' : pd.DataFrame(), #Bigger by far
}

def _id_to_tuple(id, run_id) : 
    return (run_id, id)

run_id = 0
cell_id_max = 0

if len(DROP_RNA) != len(PATH_LIST) and len(DROP_RNA) != 0 :
    raise ValueError("len PATH LIST and DROP RNA must match. check input parameters")
elif len(DROP_RNA) == 0 :
    DROP_RNA = [[]]*len(PATH_LIST)

for acquisition_path in PATH_LIST:
    if not acquisition_path.endswith('/') : acquisition_path += '/'
    id_maker = functools.partial(_id_to_tuple, run_id=run_id)

    run = acquisition_path.replace('/result_tables/','').split('/')[-1]

    for table in tables.keys() :
        new_table = pd.read_feather(acquisition_path + table + '.feather')
        new_table['acquisition_id'] = new_table['acquisition_id'].apply(id_maker)
        
        if table == 'Acquisition' :
            new_table['run'] = run
        else:
            new_table['cell_id'] += cell_id_max

        tables[table] = pd.concat([
            tables[table],
            new_table
        ])
    cell_id_max = tables['Cell']['cell_id'].max()
    run_id +=1

Acquisition = tables['Acquisition']
Cell = tables['Cell']
Spots = tables['Spots']
del tables
Acquisition = Acquisition.reset_index(drop=True)

Acquisition['couple'] = list(zip(Acquisition['rna1'], Acquisition['rna2']))

#Dropping acquisition from drop lists
for acquisition_path, drop_list in zip(PATH_LIST, DROP_RNA):
    run = acquisition_path.replace('/result_tables/','').split('/')[-1]
    drop_idx = Acquisition[(Acquisition['run'] == run) & (Acquisition['couple'].isin(drop_list))].index

    count = len(Acquisition)
    Acquisition = Acquisition.drop(drop_idx)
    print("{0} acquisitions dropped from passed drop list : {1}.".format(count - len(Acquisition), drop_list))

#Dropping failled/skipped acquisition
count = len(Acquisition)
if EXPECTED_RNA1 : Acquisition = Acquisition.loc[Acquisition['rna1_deconvolution_sucess']]
if EXPECTED_RNA2 : Acquisition = Acquisition.loc[Acquisition['rna2_deconvolution_sucess']]
if EXPECTED_SUNTAG : Acquisition = Acquisition.loc[Acquisition['suntag_deconvolution_sucess']]
print("{0} empty acquisitions dropped.".format(count - len(Acquisition)))

Acquisition = Acquisition.reset_index(drop=True)
Acquisition = Acquisition.drop('couple', axis=1)


if not 'treatment' in Acquisition.columns :
    Acquisition['treatment'] = "untreated"

Cell = pd.merge(Cell, Acquisition.loc[:,['acquisition_id']], on= 'acquisition_id', how='inner')
Spots = pd.merge(Spots, Acquisition.loc[:,['acquisition_id']], on= 'acquisition_id', how='inner')

"""

1.5 Data Integrity

"""
os.makedirs(PATH_OUT, exist_ok=True)

integrity_df, fail_view = _colloc_count_consistency(Cell)
fail_view.to_excel(PATH_OUT + '/detail_colloc_count_consistency.xlsx')
integrity_df.to_excel(PATH_OUT + '/sumup_colloc_count_consistency.xlsx')
filename_duplicate = _filename_duplicate_check(Acquisition=Acquisition)
filename_duplicate.to_excel(PATH_OUT + '/filename_duplicate.xlsx')

if len(fail_view) > 0 : 
    print("Warning : Inconsistency in spot coloc")
    print(integrity_df)


"""

2. Data cleaning

"""
print("cleaning data...")

#treatment
if 'treatment' in Acquisition :
    Acquisition['treatment'] = Acquisition['treatment'].str.lower()
    Acquisition['treatment'] = Acquisition['treatment'].replace({
        "without puromycin" : 'untreated',
        "with puromycin" : "puromycin"
    })
    print("Remaining treatments : ")
    print(Acquisition.value_counts(subset='treatment'))

else : Acquisition['treatment'] = 'untreated'

if len(Acquisition['treatment'].unique()) != EXPECTED_TREATMENT_NUMBER :
    print("Treatment\n",Acquisition['treatment'].unique())
    raise ValueError("Incorrect number of treatment remaining after cleaning")

#RNA1
Acquisition['rna1'] = Acquisition['rna1'].str.upper()
Acquisition['rna1'] = Acquisition['rna1'].replace({
    'FLAP' : 'FLAP-only'
})

print("Remaining rna1 : ")
print(Acquisition.value_counts(subset='rna1'))

if not Acquisition.query("rna1 not in {0}".format(EXPECTED_RNA1)).empty: 
    raise ValueError("Incorrect rna1 remaining after cleaning")

#RNA2
Acquisition['rna2'] = Acquisition['rna2'].str.upper()
Acquisition['rna2'] = Acquisition['rna2'].replace({
    '' : 'EMPTY',
    'ONLY' : 'EMPTY',
})

print("Remaining rna2 : ")
print(Acquisition.value_counts(subset='rna2'))

if not Acquisition.query("rna2 not in {0}".format(EXPECTED_RNA2)).empty : 
    raise ValueError("Incorrect rna2 remaining after cleaning")


#Adding treatment to Cell df
before_cellnumber = len(Cell)
Cell = pd.merge(
    left=Cell,
    right=Acquisition.loc[:,['acquisition_id', 'treatment', 'rna2']],
    how='inner',
    on= 'acquisition_id'
)

assert len(Cell) == before_cellnumber

"""

3. Number of rna1 foci histograms

"""
os.makedirs(PATH_OUT + '/rna1_foci_number_histogram', exist_ok=True)

fig = plt.figure(figsize=(10, 20))
axes = fig.subplots(2,1)
treatments = Acquisition['treatment'].unique()
colors = ['green', 'blue']

for ax, treatment, color in zip(axes, treatments, colors) :
    sub_data = Cell.loc[Cell['treatment'] == treatment]
    ax.hist(sub_data['rna1_cluster_number'], bins=30, label= treatment, alpha= 0.7, color= color)
    ax.set_xlabel('foci number')
    ax.set_ylabel('count')
    ax.set_ylim(bottom=1)
    ax.set_xlim(left=0)
    # ax.set_yscale('log')
    ax.legend()
    
#Same scale
xmin1,xmax1,ymin1,ymax1 = axes[0].axis()
xmin2,xmax2,ymin2,ymax2 = axes[1].axis()
xmin = min(xmin1,xmin2)
xmax = max(xmax1,xmax2)
ymin = min(ymin1,ymin2)
ymax = max(ymax1,ymax2)
axes[0].axis([xmin,xmax,ymin,ymax])
axes[1].axis([xmin,xmax,ymin,ymax])
plt.tight_layout()

fig.savefig(PATH_OUT + '/rna1_foci_number_histogram/rna1_cluster_number_histogram.png')
plt.close()

fig = plt.figure(figsize=(10, 20))
axes = fig.subplots(2,1)
treatments = Acquisition['treatment'].unique()
colors = ['green', 'blue']


for ax, treatment, color in zip(axes, treatments, colors) :
    sub_data = Cell.loc[Cell['treatment'] == treatment]
    ax.hist(sub_data['rna1_clustered_number'], bins=30, label= treatment, alpha= 0.7, color= color)
    ax.set_xlabel('APC in foci number')
    ax.set_ylabel('count')
    ax.set_ylim(bottom=1)
    ax.set_xlim(left=0)
    ax.set_yscale('log')
    ax.legend()
    
#Same scale
xmin1,xmax1,ymin1,ymax1 = axes[0].axis()
xmin2,xmax2,ymin2,ymax2 = axes[1].axis()
xmin = min(xmin1,xmin2)
xmax = max(xmax1,xmax2)
ymin = min(ymin1,ymin2)
ymax = max(ymax1,ymax2)
axes[0].axis([xmin,xmax,ymin,ymax])
axes[1].axis([xmin,xmax,ymin,ymax])
plt.tight_layout()

os.makedirs(PATH_OUT + '/APC_spot_in_foci_number_histogram', exist_ok=True)
fig.savefig(PATH_OUT + '/APC_spot_in_foci_number_histogram/APC_spot_in_foci_number_histogram.png')
plt.close()


"""

4. Droping phenotype negative cells.

"""

if DROP_NEGATIVE_CELLS : 

    print(Cell.value_counts(subset='rna2'))

    print(Cell.columns)
    drop_idx = Cell[
        (Cell['rna1_cluster_number'] < MIN_NUMBER_FOCI) & (Cell['treatment'] == 'untreated') & (Cell['rna2'] != "ARHGAP21")
        # ((Cell['rna1_clustered_number'] <= 200)) & (Cell['treatment'] == 'untreated')
        ].index

    drop_idx_ARGAP = Cell[
        (Cell['rna1_cluster_number'] < MIN_NUMBER_FOCI_ARGAP) & (Cell['treatment'] == 'untreated') & (Cell['rna2'] == "ARHGAP21")
        # ((Cell['rna1_clustered_number'] <= 200)) & (Cell['treatment'] == 'untreated')
        ].index

    Cell = Cell.drop(drop_idx, axis=0)
    Cell = Cell.drop(drop_idx_ARGAP, axis=0)
    print(Cell.value_counts(subset='rna2'))

    print("{0} untreated cells were drop with min cluster limit : {1}".format(len(drop_idx), MIN_NUMBER_FOCI))

"""

5. Saving results

"""

print("Saving results...")

datetime = get_datetime()

with open(PATH_OUT + "/merge_log.txt", 'w+') as log :
    log.write("Merge log\n")
    log.write("This log was created on {0}\n".format(datetime))

    log.write("\nMerged runs:\n")
    for run in PATH_LIST : log.write('"{0}"\n'.format(run))

    log.write("\nParameters:\n")
    log.write("Min number of rna1 cluster for phenotype positive foci : "+ str(MIN_NUMBER_FOCI))
    log.write("Min number of rna1 cluster for phenotype positive foci : "+ str(MIN_NUMBER_FOCI_ARGAP))
    log.write("\nExpected rna1 : "+ str(EXPECTED_RNA1))
    log.write("\nExpected rna2 : "+ str(EXPECTED_RNA2))
    log.write("\nNumber of treatment expected : "+ str(EXPECTED_TREATMENT_NUMBER))

Acquisition.reset_index(drop=True).to_feather(PATH_OUT + '/Acquisition.feather')
Cell.reset_index(drop=True).to_feather(PATH_OUT + '/Cell.feather')
Spots.reset_index(drop=True).to_feather(PATH_OUT + '/Spots.feather')
Acquisition.reset_index(drop=True).to_excel(PATH_OUT + '/Acquisition.xlsx')

print("done")