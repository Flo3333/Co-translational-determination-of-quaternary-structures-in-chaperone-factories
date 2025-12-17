import czifile as czi
import bigfish.stack as stack
import numpy as np
import os

"""
Specific script for /media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/Cinétique des foyers/240206_hct_h9_polr2a dat
Aim to clean useless channels and rename file accordingly for past pipeline.

"""



path = '/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/Cinétique des foyers/240206_hct_h9_polr2a/'
out_path = '/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/Cinétique des foyers/240206_hct_h9_polr2a_cleaned/'

file_list = os.listdir(path)

for filename in file_list :
    im = czi.imread(path + filename)
    shape = im.shape
    channel_number = shape[0]
    new_im = np.zeros(shape = (2, shape[1], shape[2], shape[3]))
    new_im[0] = im[0,:,:,:,0]
    new_im[1] = im[-1,:,:,:,0]

    stack.save_image(new_im.astype(np.uint16), out_path + filename.replace('.czi',".tiff"))