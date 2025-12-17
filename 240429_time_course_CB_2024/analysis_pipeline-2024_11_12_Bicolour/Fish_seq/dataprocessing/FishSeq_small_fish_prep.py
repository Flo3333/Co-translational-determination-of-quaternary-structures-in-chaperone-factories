RUN_PATH = "/media/floricslimani/SSD4To/SSD_floricslimani/Fish_seq/Davide/2024-08-12 - SeqFISH - HeLa - Puro - R2TP1-2_Run7"


import os
import numpy as np
from pbwrap.utils import open_image
from tqdm import tqdm
import bigfish.stack as stack

im_path = RUN_PATH + "/FISH_Z-stacks/Location-01/img000_000_000000_0000000000.ome.tif"
os.makedirs(RUN_PATH + "/small_fish_sample/",exist_ok=True)

image_stack = open_image(im_path)

cycle = 0
for image in tqdm(image_stack) :
    image = np.swapaxes(image, axis1=0, axis2=-1)
    for index, color in enumerate(image) :
        stack.save_image(image=color, path= RUN_PATH + "/small_fish_sample/cycle_{0}_color_{1}.tiff".format(cycle,index))
    cycle +=1