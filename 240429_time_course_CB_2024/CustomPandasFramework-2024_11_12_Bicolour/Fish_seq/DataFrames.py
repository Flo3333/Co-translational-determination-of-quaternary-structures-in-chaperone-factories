import pandas as pd


def create_EmptyColocalisation() :
    header = ['id', 'AcquisitionId', 'rna1_number',  'rna2_number', 'rna1_rna2_colocalisation_count', 'rna2_rna1_colocalisation_count']
    return pd.DataFrame(columns= header)

def create_EmptyClusteredSpots() :
    header = ['id', 'AcquisitionId', 'z', 'y', 'x', 'cluster_id', 'colocalising_suntag_number']
    return pd.DataFrame(columns= header)