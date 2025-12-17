import pandas as pd

def open_tables(result_path, table_names: list, rna_list= None, rna_key = 'rna', main_table_name = 'Acquisition') :
    """
    rna_key is the key to rna names in the main_table_name, also all tables have a ref to main_table_name following the name : key = main_table_name + 'Id'.
    """
    
    if not result_path.endswith('/') : result_path += '/'
    
    if type(rna_list) != type(None) :
        if main_table_name in table_names : table_names.remove(main_table_name)
        table = pd.read_feather(result_path + "{0}.feather".format(main_table_name))
        query = table.query('{0} in {1}'.format(rna_key, rna_list))
        keep_id = query['id']
        keep_idx = query.index
        tables = {main_table_name : table.loc[keep_idx,:].reset_index(drop=True)}
         
        other_tables = {
            name : pd.read_feather(result_path + "{0}.feather".format(name)).query('{0} in {1}'.format(main_table_name + 'Id', list(keep_id)), inplace=False).reset_index(drop=True)
            for name in table_names
        }
        tables.update(other_tables)

    else :
        tables = {name : pd.read_feather(result_path + "{0}.feather".format(name)) for name in table_names}

    return tables

def merge_tables(tables: dict, new_tables: dict) :

    for table_name, new_table in new_tables.items() :
        table = tables.get(table_name)
        tables[table_name] = pd.concat([table, new_table], axis= 0)
    
    return tables

def number_common_elements(list1, list2) :
    return sum([elmt1 in list2 for elmt1 in list1])


def look_for_duplicates(table:pd.DataFrame, new_table: pd.DataFrame, key_list) :

    duplicates = {}

    for key in key_list :

        if not key in table.columns or not key in new_table.columns : 
            print("{0} not found in tables when looking for duplicates.".format(key))
            continue

        key_elements = table[key].unique()
        key_new_elements = new_table[key].unique()
        duplicates_number = number_common_elements(key_elements, key_new_elements)

        if duplicates_number > 0 : duplicates.update({key : duplicates_number})
    
    return duplicates

def save_merged_tables(tables: 'dict[pd.DataFrame]', path_out: str) :

    if not path_out.endswith('/') : path_out += '/'
    
    for table_name, table in tables.items() :
        table.reset_index(drop= True).to_feather(path_out + table_name +'.feather')