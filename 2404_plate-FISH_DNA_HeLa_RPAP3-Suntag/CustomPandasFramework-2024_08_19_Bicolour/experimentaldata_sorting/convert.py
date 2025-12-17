"""Convert cleaned data from experimental files cleaning to pandas DataFrames that can be stored"""

from .DataFrames import newframe_Oligopools, newframe_Plates



def convert_Oligopools(TargetGenes) :
    """Makes a DataFrame that can be added to a Oligopools DataFrame from a TargetGenes
    
    Parameters
    ----------
        list_Plan : pd.DataFrame
            Results from cleaning.TargetGenes
            
    Return
    ------
        frame_Oligopools : pd.DataFrame
    
    """

    Oligopools_temp = TargetGenes.value_counts(subset= "n° oligopool")
    Oligopools_temp = Oligopools_temp.reset_index(drop= False).drop(0,1)
    Oligopools_temp = Oligopools_temp.rename(columns= {"n° oligopool" : "name"})

    frame_Oligopools = newframe_Oligopools()

    name_list = list(Oligopools_temp["name"])
    id_list = list(range(0, len(name_list)))

    frame_Oligopools["id"] = id_list
    frame_Oligopools["name"] = name_list

    return frame_Oligopools 




def convert_Plates(list_Plan) :
    """Makes a DataFrame that can be added to a Plates DataFrame from a list_PlatePlan
    
    Parameters
    ----------
        list_Plan : pd.DataFrame
            Results from cleaning.list_PlatePlan
            
    Return
    ------
        Plates : pd.DataFrame
    
    """
    new_Plates = newframe_Plates()

    name_list = list(list_Plan["name"])
    id_list = list(range(0,len(name_list)))

    new_Plates["id"] = id_list
    new_Plates["name"] = name_list

    return new_Plates



