import inspect
import datetime as dt
from .utils import check_parameter

def add_datetime_folder(path:str) :
    check_parameter(path = str)
    datetime = get_datetime()
    if not path.endswith('/') : 
        path_out = path + '/'
    else : path_out = path
    return path_out + datetime + '/'


def get_datetime():
    """
    Return current date time with format : 'YYYYMMDD_HH-mm-SS'
    """
    return dt.datetime.now().strftime("%Y%m%d_%H-%M-%S")



def _get_varname(var):
    #To be used  within a function.
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var][0]



def print_parameters(path_out, *parameters, printDateTime= True):
    """
    Print parameters into a .txt file
    
    Parameters
    ----------
        path_out : str
            full_path to saving location.
        *parameters : args(int, float, str)
            parameters to print into text file.
    """

    if path_out[len(path_out)-1] == '/' : path_out += 'parameters.txt'
    elif path_out[len(path_out)-4 : len(path_out)] != '.txt' : path_out += '.txt'

    parameter_file = open(path_out, "w")

    #Header
    parameter_file.write("PARAMETERS\n")
    parameter_file.write("\n############\n")

    if printDateTime:
        datetime = get_datetime()
        parameter_file.write(datetime)
        parameter_file.write("############\n\n")
    
    lines= []
    for parameter in parameters :
        name = _get_varname(parameter)
        lines += ["{0} : {1}\n".format(name, parameter)]
    parameter_file.writelines(lines)

    parameter_file.close()


def print_dict(dic: dict, path_out):
    if path_out[len(path_out) -1] == '/' : path_out += 'dic.txt'
    elif path_out[len(path_out)-4 : len(path_out)] != '.txt' : path_out += '.txt'

    lines = []
    for key, value in dic.items():
        lines += "{0} : {1}\n".format(key, value)


    dict_file = open(path_out, "w")
    dict_file.writelines(lines)
    dict_file.close()

def dict_to_lines(di : dict) :
    for key, value in di.items() : 
        yield "{0} : {1}\n".format(key, value)