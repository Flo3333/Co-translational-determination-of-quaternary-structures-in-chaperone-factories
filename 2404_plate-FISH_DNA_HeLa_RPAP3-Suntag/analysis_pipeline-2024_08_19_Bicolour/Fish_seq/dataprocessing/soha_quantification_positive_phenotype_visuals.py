"""
Data pre-processing script for Soha quantifications
Pipeline : Fish_seq / BC_clusterkinetics_pipeline


This analysis restricts cells (counted as phenotype_positive) to cell showing numerous clusters from APC channel (rna1).
This script aims at giving a visual of them.

"""

import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import bigfish.stack as stack
import czifile as czi
import bigfish.plot as plot

CHANNEL = 1

merge_tables_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/Soha quantification/Merged_results/"

dataset_path_list = [
    "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/Soha quantification/b-cat_APC_231211"

]

Cell = pd.read_feather(merge_tables_path + "Cell.feather", columns=['cell_id','rna1_bbox','acquisition_id', 'rna1_cluster_number'])
Acquistion = pd.read_feather(merge_tables_path + "Acquisition.feather", columns=['dataset','treatment','filename','acquisition_id', 'date', 'rna1', 'rna2'])

Cell['acquisition_id'] = Cell['acquisition_id'].apply(tuple)
Acquistion['acquisition_id'] = Acquistion['acquisition_id'].apply(tuple)

df = pd.merge(Cell, Acquistion, on='acquisition_id')
assert len(df) == len(Cell), "duplication or deletion of cells."
df = df.loc[df['treatment'] == "untreated"]


path = merge_tables_path + "/positive_phenotype_visuals/"
os.makedirs(path, exist_ok=True)
extract = df.groupby(['rna1','rna2','treatment'])['cell_id'].count().rename('count')
print(extract)
extract.to_excel(path + "positive_phenotype_cell_count.xlsx")

positive_cell_number = len(df)
count = 0
for cell_index in tqdm(df.index) : 
    cell = df.loc[cell_index]
    filename = cell.at['filename']

    found = False
    image_path=""
    for dataset_path in dataset_path_list :
        if os.path.isfile(dataset_path + '/' + filename) :
            image_path = dataset_path + '/' + filename
            found = True
            break
    
    if not found : 
        continue

    image = czi.imread(image_path)
    image = np.squeeze(image[CHANNEL]) #remove dim of len 1
    image = stack.mean_projection(image)
    
    miny,minx,maxy,maxx = cell.at['rna1_bbox']
    foci_number = cell.at['rna1_cluster_number']
    cell_id = cell.at['cell_id']

    image = image[miny:maxy,minx:maxx]
    path = merge_tables_path + "/positive_phenotype_visuals/" + dataset_path.split('/')[-1] + '/'
    os.makedirs(path, exist_ok=True)
    plot.plot_images([image], titles=["foci number : {0}".format(foci_number)], contrast=True, show=False, path_output= path + "{0}_{1}".format(filename, cell_id))
    count +=1

print("Done. {0}/{1} cells found.".format(count, positive_cell_number))


