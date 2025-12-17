# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>


"""

This subpackages handles framework related to the bicolour work : 
 1. colocalisation project
 2. cluster kinetics project

"""

from .plot import _coloc_get_target_columns, _compute_coloc_fractions
from .colocalisation import compute_colocalisation, get_coloc_objects_names,get_coloc_object_dict