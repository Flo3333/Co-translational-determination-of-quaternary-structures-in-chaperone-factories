""" Making DataFrames structures aimed for expermintal datal sorting in our team."""

import pandas as pd


#Empty Dataframe_making
def newframe_Oligopools():
    """Returns an empty pandas DataFrame with expected data shape for storing Oligopool data acquired experimentally.
    
    Returns
    -------
    
        Oligopool : pd.Dataframe
            Empty Oligopool dataframe.

    Datashape
    ---------
        "id"
        "name"
    """

    Oligopool = pd.DataFrame({       
        "id" : [],
        "name" : [],
    })

    return Oligopool

def newframe_Plates() :
    """Returns an empty pandas DataFrame with expected data shape for storing Plate data acquired experimentally.
    
    Returns
    -------
    
        Plates : pd.Dataframe
            Empty TargetGenes dataframe.

    Datashape
    ---------
        "id"
        "name"
        "Oligopoolid"

    """

    Plates = pd.DataFrame({
        "id" : [],
        "name" : [],
        "Oligopoolid" : []
    })

    return Plates

def newframe_Sets():
    """Returns an empty pandas DataFrame with expected data shape for storing Sets data acquired experimentally."""

    res = pd.DataFrame({
        "id" : [],
        "name" : []

    })

def newframe_Wells() : # probablement useless
    """Returns an empty pandas DataFrame with expected data shape for storing Wells data acquired experimentally.
    
    Returns
    -------
        Wells : pd.Dataframe
            Empty Wells dataframe.

    DataShape
    ---------
        "id"
        "OligopoolId"
        "PlateId"
        "TargetGeneId"
        "coordinates"
        "other relevant data to add"
    """

    Wells = pd.DataFrame({
        "id" : [],
        "OligopoolId" : [],
        "PlateId" : [],
        "TargetGeneId" : [],
        "coordinates" : [],
        "other relevant data to add" : []
    })

    return Wells

def newframe_TargetGenes() :
    """Returns an empty pandas DataFrame with expected data shape for storing TargetGenes data acquired experimentally.
    
    Returns
    -------
    
        TargetGenes : pd.Dataframe
            Empty TargetGenes dataframe.
    """

    TargetGenes = pd.Dataframe({
        "id" : [],
        "OligopoolId" : [],
        "PlateId" : [],
        "name" : [],
        "other cara to add" : []
    })

    return TargetGenes

def newframe_Sequences():
    """Returns an empty pandas DataFrame with expected data shape for storing Sequences data acquired experimentally.
    
    Returns
    -------
    
        Sequences : pd.Dataframe
            Empty Squences dataframe.
    """

    Sequences = pd.DataFrame({
        "id" : [],
        "OligopoolId" : [],
        "GenesId" : [],
        "name" : [],
        "other relevant data" : []

    })
    
    return Sequences



