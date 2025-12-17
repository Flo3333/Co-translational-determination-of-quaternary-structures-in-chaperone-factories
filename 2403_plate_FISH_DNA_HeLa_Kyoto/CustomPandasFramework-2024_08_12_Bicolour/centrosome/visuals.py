import bigfish.plot as plot
import numpy as np

def centrosome_cell_plot(cell: dict, new_Cell, path_output:str, boundary_size = 1, show=False, contrast= True) :
    
    cell_coord = cell["cell_coord"]
    nuc_coord = cell["nuc_coord"]
    cell_mask = cell["cell_mask"]
    nucleus_mask = cell["nuc_mask"]
    centrosome_image = cell["centrosome_mip"]
    label = cell["cell_id"]
    df = new_Cell.copy().reset_index(drop=True)
    spot_coords = cell['rna_coord']
    centrosome_coords = np.array(df.at[0,"centrosome_coords_local"], dtype= int)
    
    if path_output.endswith('/') : 
        path = path_output + 'centrosome_cell_{0}'.format(label)
    
    plot.plot_cell(ndim= 3, 
                   cell_coord= cell_coord,
                   nuc_coord= nuc_coord,
                   image=centrosome_image,
                   rna_coord= spot_coords,
                   other_coord= centrosome_coords,
                   cell_mask=cell_mask,
                   nuc_mask = nucleus_mask,
                   show=show,
                   contrast=contrast,
                   boundary_size=boundary_size,
                   path_output=path
                   )