"""
Try to open all images in a folder and output file with all the filenames that it couldn't open.
"""

import bigfish.stack as stack
import os
from tqdm import tqdm

PATH = '/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/input/'
if not PATH.endswith('/') : PATH += '/'
filename_list = os.listdir(PATH)

error_filenames = []
print(filename_list[1817])
for file in tqdm(filename_list[1817:], desc= "Progression : ", total= len(filename_list)) :
    try :
        stack.read_image(PATH + file)
    except Exception :
        print(file)
        filename_list.append(file)

with open(PATH + "corrupted_file.txt", 'x') as file :
    file.writelines(error_filenames)


print("done.")