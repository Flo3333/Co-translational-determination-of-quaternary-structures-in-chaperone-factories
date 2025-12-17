class MissingColumnError(Exception) : 
    pass

class MergeError(Exception):
    pass

class InputError(Exception):
    """
    Error was encountered while reading files in input.
    """
    pass