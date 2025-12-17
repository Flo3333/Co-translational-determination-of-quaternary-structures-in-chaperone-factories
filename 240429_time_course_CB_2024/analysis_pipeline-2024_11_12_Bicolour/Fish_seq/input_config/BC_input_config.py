"""
Config file for clusterkinetics pipeline and Fish_seq_pipeline allowing it to handle variety of file namings and channel availabilities.
For each data set there are 4 lines to add in this config file. Data set names start with acquisition date followed by data set name.


1. REGEX : string
    --> r"{ENTER THE REGEX}" : regex used to read file name there must be as much capture groups () than columns passed in 2.

2. COLUMNS : list[str]
    --> names of data columns created in the Acquisition table from filename informations in regex capture groups (1).
    FOR A SET OF MONOCHANNEL IMAGE : 'channel' key is expected

3. CHANNELS : dict{str:int}
    3 keys are needed : 'dapi', 'rna1', 'rna2' for clusterkinetics
    2 keys are expected for Fish_seq : 'nuc', 'cyto + 1key per 'SINGLE_MOLECULE' to detect (5.)

    --> FOR MULTICHANNEL IMAGE : Each key must have the reference to the corresponding channel in the image.
        if no rna2 pass None or make a regex that will give rna2 empty.
    
    --> FOR A SET OF MONOCHANNEL IMAGE : Each key refers to the 'channel' column found by the regex.
        if no rna2 pass None or make a regex that will give rna2 empty.

4. GROUP_KEYS : list[str]
    keys (from columns) on which to group fov detection threshold computation.
    
        example : 
            if you have different rna and treatment and you want one threshold per couple (rna, treament) put in config :
            ['rna', 'treatment']
            if you want one threshold per rna but the same on all treatment :
            ['rna']

5. SINGLE_MOLECULE
    list of names for single molecules type to detect.
    Examples :
        --> rna1, rna2, suntag
        --> cy3, cy5, Alexa648,
        --> c1, c2, c3, c4, c5...
        ...
    
    Each element in this list is expected to be found as a key in 'CHANNELS' (3.).
    
"""

#1.
REGEX = {
    # dataset : r"(\d{6})_(\w+-\w+(?:-\w+)?)_(\w+)_(\w+)-(\d{2})\.tiff"
    "240206_hct_h9_polr2a" : r"(\d{6})_(\w+-\w+(?:-\w+)?)_(\w+)_(\w+)-(\d{2})\.tiff"
    ,"2023-06-29" : r"(\d{6})_(unknown)_(\w+)_(\w+)_([-\w]+)-(\d{2})\.czi"
    ,"2023-05-24" : r"(\d{6})_(unknown)_(\w+)_(\w+)_([-\w]+)-(\d{2})\.czi"
    ,"2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto" : r"(r\w{5})-(Kyoto)_(\w+?)_([\w-]+?)(?:_(\w+))?f(\d{2})-(.+)-sk1fk1fl1.tiff"
    ,"b-cat quantif 230723 n1 b-cat bac APC" : r"(\d{6})\s(n\d+)\s([\w-]+).*fitc\s([\w-]+).*cy3\s(.+)-(\d{2}).czi"
    ,"b-cat quantif 230725 n1 b-cat bac APC" :r"(\d{6})\s(n\d+)\s([\w-]+).*fitc\s([\w-]+).*cy3\s(.+)-(\d{2}).czi"
    ,"b-cat quantif 230724 n1 b-cat bac APC" : r"(\d{6})\s(n\d+)\s([\w-]+).*FITC\s([\w-]+).*cy3\s(.+)-(\d{2}).czi"
    ,"b-cat quantif 230905 n1 b-cat bac APC" : r"(\d{6})\s(n\d+)\s([\w-]+).*fitc\s([\w-]+).*cy3\s(.+)-(\d{2}).czi"
    ,"b-cat_APC_231211" : r"(\d{6})\s([\w-]+).*FITC\s([\w-]+).*cy3-(\d{2}).czi"
    ,"240429_time_course_CB" : r"(\d{6})_([^_]+)_([^_]+)_(.+)-(\d{2})\..*"
    ,"2404_plate-FISH_DNA_HeLa_RPAP3-Suntag" : r"r\w+-([^_]+)_([^_]+)_([^_]+)_?([^_]*)f(\d{2})-([^-]+)-.*.tiff"
    ,"240627_centrosome_golgi_RNA_FISH" : r"r\w+-([^_]+)_([^_]+)_([^_]+)_([^_]+)f(\d{2})-([^-]+)-.*.tiff"
}

#2.
COLUMNS = {
    #dataset :['date', 'cell_line', 'rna1', 'treatment', 'fov']
    "240206_hct_h9_polr2a" :['date', 'cell_line', 'rna1', 'treatment', 'fov_number']
    ,"2023-06-29" : ['date', 'cell_line', 'rna1', 'rna2', 'treatment', 'fov_number']
    ,"2023-05-24" : ['date', 'cell_line', 'rna1', 'rna2', 'treatment', 'fov_number']
    ,"2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto" : ['well', 'cell_line', 'treatment', 'rna1', 'rna2', 'fov_number', 'channel']
    ,"b-cat quantif 230723 n1 b-cat bac APC" : ['date', 'experiment', 'rna1', 'rna2', 'treatment', 'fov']
    ,"b-cat quantif 230725 n1 b-cat bac APC" : ['date', 'experiment', 'rna1', 'rna2', 'treatment', 'fov']
    ,"b-cat quantif 230724 n1 b-cat bac APC" : ['date', 'experiment', 'rna1', 'rna2', 'treatment', 'fov']
    ,"b-cat quantif 230905 n1 b-cat bac APC" : ['date', 'experiment', 'rna1', 'rna2', 'treatment', 'fov']
    ,"b-cat_APC_231211" : ['date', 'rna1', 'rna2', 'fov']
    ,"240429_time_course_CB" : ['date','rna1','cell_line','treatment','fov'] #RNA 1 unique in this set
    ,"2404_plate-FISH_DNA_HeLa_RPAP3-Suntag" : ['suntag','treatment','rna1','rna2','fov','channel']
    ,"240627_centrosome_golgi_RNA_FISH" : ['cell_line','rna1','rna2','treatment','fov','channel']
}

#3.
CHANNELS = {
    #set of monochannel --> dataset : {'dapi' : 'DAPI', 'rna1' : 'cy3', 'rna2' : 'cy5', suntag : 'suntag'}
    #multichannel --> dataset : {'dapi' : 0, 'rna1' : 1, 'rna2' : 2, 'suntag' : 3}
    "240206_hct_h9_polr2a" : {"dapi" : 0,"rna1" : 1,"rna2" : None}
    ,"2023-06-29" : {"dapi" : 0,"rna1" : 2,"rna2" : 3}
    ,"2023-05-24" : {"dapi" : 0,"rna1" : 2,"rna2" : 3}
    ,"2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto" : {"dapi" : 'DAPI', "rna1" : 'Cy3', "rna2" : "Alexa 647"}
    ,"b-cat quantif 230723 n1 b-cat bac APC" : {'dapi' : 0, 'rna1' : 1, 'rna2' : 2 }
    ,"b-cat quantif 230725 n1 b-cat bac APC" : {'dapi' : 0, 'rna1' : 1, 'rna2' : 2 }
    ,"b-cat quantif 230724 n1 b-cat bac APC" : {'dapi' : 0, 'rna1' : 1, 'rna2' : 2 }
    ,"b-cat quantif 230905 n1 b-cat bac APC" : {'dapi' : 0, 'rna1' : 1, 'rna2' : 2 }
    ,"b-cat_APC_231211" : {'dapi' : 0, 'rna1' : 1, 'rna2' : 2 }
    ,"240429_time_course_CB" : {'dapi' : 0, 'rna1' : 1}
    ,"2404_plate-FISH_DNA_HeLa_RPAP3-Suntag" : {'dapi' : 'DAPI', 'rna1' : 'Cy3', 'suntag' : 'EGFP', 'rna2' : 'Alexa 647'}
    ,"240627_centrosome_golgi_RNA_FISH" : {'dapi' : 'DAPI', 'rna1' : 'Cy3', 'rna2' : 'EGFP'}
}

#4.
GROUP_KEYS = {
    #dataset : ['date','rna1','treatment']
    "240206_hct_h9_polr2a" : ['date', 'rna1']
    ,"2023-06-29" : ['date', 'rna1','rna2']
    ,"2023-05-24" : ['date', 'rna1','rna2']
    ,"2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto" : ['rna1', 'rna2']
    ,"b-cat quantif 230723 n1 b-cat bac APC" : ['date', 'rna1', 'rna2']
    ,"b-cat quantif 230725 n1 b-cat bac APC" : ['date', 'rna1', 'rna2']
    ,"b-cat quantif 230724 n1 b-cat bac APC" : ['date', 'rna1', 'rna2']
    ,"b-cat quantif 230905 n1 b-cat bac APC" : ['date', 'rna1', 'rna2']
    ,"b-cat_APC_231211" : ['date', 'rna1', 'rna2']
    ,"240429_time_course_CB" : ['cell_line']
    ,"2404_plate-FISH_DNA_HeLa_RPAP3-Suntag" : ['rna1','rna2', 'suntag']
    ,"240627_centrosome_golgi_RNA_FISH" : ['cell_line','rna1','rna2']
}

#5.
SINGLE_MOLECULES = {
    #dataset : ['cy3', 'cy5', 'suntag']

}

#6. OPTIONAL
THRESHOLDS = {
    #dataset : {group_keys_tuple : thresholds_tuple (example : 70,None,50,None) would work with 4 type of single molecules to detect ordered like in single_molecules},
   "b-cat quantif 230723 n1 b-cat bac APC" : {
                                            ('230723', 'b-cat', 'ires') : (450,1000,None),
                                            ('230723', 'b-cat', 'ARHAP21h') : (450,350,None), 
                                            ('230723', 'b-cat', 'JUP') : (450,325,None), 
                                            ('230723', 'b-cat', 'POLR2b') : (450,350,None)
                                            },
    "b-cat quantif 230725 n1 b-cat bac APC" : {
                                            ('230725', 'b-cat', 'ires') : (450,1000,None),
                                            ('230725', 'b-cat', 'ARHAP21h') : (450,350,None), 
                                            ('230725', 'b-cat', 'JUP') : (450,325,None), 
                                            ('230725', 'b-cat', 'POLR2B') : (450,350,None),
                                            ('230725', 'b-cat', 'ARGAP21') : (450,350,None)
                                            },
    "b-cat quantif 230724 n1 b-cat bac APC" : {
                                            ('230724', 'b-cat', 'ARHAP21h') : (450,350,None), 
                                            ('230724', 'b-cat', 'POLR2B') : (450,350,None),
                                            },

    "b-cat quantif 230905 n1 b-cat bac APC" : {
                                            ('230905', 'b-cat', 'ires') : (900,700,None),
                                            },
    "b-cat_APC_231211" : {
                        ('231118', 'b-cat', 'CASK') : (500, 1103, None),
                        ('231118', 'b-cat', 'DLG1') : (500, 1200, None),
                        ('231118', 'b-cat', 'MYO5B') : (500, 750, None),
                        ('231211', 'b-cat', 'AMER1') : (500, 992, None),
                        ('231211', 'b-cat', 'MYO5B') : (500, 894, None),
                        ('231211', 'b-cat', 'DLG1') : (500, 1215, None),
                        ('231211', 'b-cat', 'CASK') : (500, 796, None),
    },
    "2024-04-04_2403_plate_FISH_DNA_HeLa_Kyoto" : {
        ('R1A','R2A') : (39, 75, None),
        ('R2A','R1A') : (60, 80, None),
        ('R1A','R2A') : (39, 125, None),
        ('JUP','R2A') : (120, 83, None),
        ('BRD8','EP400') : (120, 130, None),
        ('EP400' , 'BRD8') : (80, 120, None),
        ('JUP','EP400') : (75,81,None),
        ('PRPF8', 'EFTUD2') : (200, 110,None),
        ('EFTUD2' , 'PRPF8') : (70,200,None),
        ('JUP','EFTUD2') : (120,120,None),
        ('R2B','R2A') : (90,120,None),
        ('R2A','R2B') : (90,120,None),
        ('JUP','') : (120,46,None),
        ('SCAF4','')  : (70,32,None),
        ('GOLGA4','') : (70,34,None),
        ('FLAP','only') : (78.75,71, None),
        ('SRCAP','') : (70,35,None),
        ('Ctl','empty') : (46.5,34,None),
        ('A-CATENIN','') : (40,40,None),
        

    },
    "240429_time_course_CB" : {
        "kyoto" : (330,None,None),
        "hek" : (430,None,None),
        "hct" : (400,None,None),
    },

    "2404_plate-FISH_DNA_HeLa_RPAP3-Suntag" : { #Suntag threshold was set as global parameter to 50
        ('JUP','R2A','RPAP3-ST') : (66,103.2, None)
        ,('JUP','EP400','RPAP3-ST') : (66,94, None)
        ,('JUP','EFTUD2','RPAP3-ST') : (66.8,95, None)
        ,('JUP','','RPAP3-ST') : (66,67.2, None)
        ,('Ctl','empty','RPAP3-ST') : (25,68.4, None)
        ,('R1A','R2A','RPAP3-ST') : (39.6,81.6, None)
        ,('BRD8','EP400','RPAP3-ST') : (70,94, None)
        ,('PRPF8','EFTUD2','RPAP3-ST') : (62.4,95, None)
        ,('GOLGA4','','RPAP3-ST') : (80,67.2, None)
        ,('R2A','R1A','RPAP3-ST') : (81.6,78, None)
        ,('R2A','R2B','RPAP3-ST') : (81.6,65, None)
        ,('EP400','BRD8','RPAP3-ST') : (80,110, None)
        ,('EFTUD2','PRPF8','RPAP3-ST') : (100,290, None)
        ,('FLAP X only','','RPAP3-ST') : (86.4,186, None)
        ,('SCAF4','','RPAP3-ST') : (30,70, None)
        ,('SRCAP','','RPAP3-ST') : (70,70, None)
    },

}
