import os

def recursive_listdir(fullpath: str)-> 'list[str]':
    if not os.path.isdir(fullpath) : raise NotADirectoryError('{0} does not refer to a directory.'.format(fullpath))
    if not fullpath.endswith('/') : fullpath += '/'
    fullname_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(fullpath)) for f in fn]
    for index in range(0,len(fullname_list)):
        count = 1
        cara_num = len(fullname_list[index]) - 1
        filename = fullname_list[index][cara_num]
        while not filename.startswith('/') :
            filename = fullname_list[index][cara_num-count] + filename
            count += 1
        fullname_list[index] = filename[1:]
    
    return fullname_list



def get_elmtindex(elmt,List) :
    """Returns index (position) of elmt in list
    
    Parameters
    ----------
        elmt : any
        list : list
    Returns
    -------
        res : int
    """
    for idx in range(0,len(List)) :
        if List[idx] == elmt : yield idx


def rm_value_from_list(value,lis : list) :
    index = get_elmtindex(value,lis)
    count = 0
    for i in index :
        if i-count == len(lis)-1 :
            lis = lis[:i-count]
            count += 1
        else:
            lis = lis[: i - count] + lis[i + 1 - count :]
            count += 1 
    return lis

def list_has_dir(path_list:list) -> 'list[bool]':
    truth_table = []
    for path in path_list :
        truth_table += [os.path.isdir(path)]
    return truth_table

