import pandas as pd
from .datashape import check_expectedcolumns, check_parameter

def get_referencement_relation(localDataFrame,local_key, ForeignFrame, foreign_key ):
    """Checks that keycolumn 1 from DataFrame1 forms a (order1, order2) referencement relationship with key column2 from DataFrame2.
    Raise Exception if it fails.
    
    Parameters
    ----------
        DataFrame1 : pd.DataFrame
        keycolumn1 : str
        order1 : str
            Either "1", "N" or "n"
        DataFrame2 : pd.DataFrame
        keycolumn2 : str
        order2 ; str
            Either "1", "N" or "n"
    """
    
    #Integrity
    check_parameter(localDataFrame= (pd.DataFrame), local_key= (str), ForeignFrame= (pd.DataFrame), foreign_key= (str))
    check_expectedcolumns(localDataFrame, [local_key])
    check_expectedcolumns(ForeignFrame, [foreign_key])
    if localDataFrame.empty : return 'empty local df'
    if ForeignFrame.empty : return 'empty foreign frame'

    #removing NULL rows
    tempForeignFrame = ForeignFrame.copy()
    tempForeignFrame = tempForeignFrame.dropna(subset= [foreign_key]).loc[:,[foreign_key]]

    tempLocalFrame = localDataFrame.copy()
    tempLocalFrame = localDataFrame.dropna(subset = [local_key]).loc[:, [local_key]]

    #group by keys
    group_fk = tempForeignFrame.value_counts(subset= foreign_key).reset_index().drop(0, axis= 1)
    group_lk = tempLocalFrame.value_counts(subset= local_key).reset_index().drop(0, axis= 1)

    #join tables
    left_join = pd.merge(left= group_lk, right= tempForeignFrame, how= 'left', left_on= local_key, right_on= foreign_key)
    right_join = pd.merge(left= tempLocalFrame, right= group_fk, how= 'right', left_on= local_key, right_on= foreign_key)


    #match testing
    if len(left_join) == len(group_lk) : foreign_order = "1"
    else : foreign_order = "N"

    if len(right_join) == len(group_fk) : local_order = "1"
    else : local_order = "N"

    if any(pd.isna(left_join[local_key])) : raise Warning("All foreign keys don't have a match in the local DataFrame.")
    if any(pd.isna(right_join[foreign_key])) : raise Warning("All local keys don't have a match in the foreign DataFrame.")

    referencement_relation = (local_order, foreign_order)
    return referencement_relation