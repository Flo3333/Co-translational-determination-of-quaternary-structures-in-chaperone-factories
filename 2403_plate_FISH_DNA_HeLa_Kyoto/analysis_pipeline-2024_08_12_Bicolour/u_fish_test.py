import torch
import numpy as np
import imageio.v3 as io
import bigfish.stack as stack
import bigfish.detection as detc
import napari
from ufish.api import UFish
from czifile import imread



path_im_test = "/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_seq/Soha quantification/b-cat quantif 230905 n1 b-cat bac APC/230905 n1 b-cat bac apc if fitc ires neo smfish cy3 with puromycin-06.czi"
ufish = UFish()
ufish.load_weights()

# image_3d = io.imread(path_im_test)
image_3d = imread(path_im_test)

print(image_3d.shape)
TEST_2D = False
TEST_3D = True

#### 2D test ####
if TEST_2D :
    voxel_size = (103,103)
    image_2d = np.max(image_3d, axis=0)[1] #max proj
    print(image_2d.shape)

    bigfish_spot = detc.detect_spots(
        images= image_2d,
        threshold=None,
        voxel_size=voxel_size,
        spot_radius=(100,100)
    )

    ufish_spot, enh_im = ufish.predict(image_2d)

    ufish_spot = np.array(list(zip(ufish_spot['axis-0'], ufish_spot['axis-1'])), dtype=int)

    viewer = napari.Viewer(title='2D test')
    viewer.add_image(image_2d)
    viewer.add_image(enh_im)
    viewer.add_points(bigfish_spot, face_color='blue', size= 5)
    viewer.add_points(ufish_spot,  face_color='green', size= 5)

    viewer.show(block=True)

if TEST_3D :
    voxel_size = (300,103,103)
    image_3d = image_3d[1,:,:,:,0]
    print(image_3d.shape)

    bigfish_spot = detc.detect_spots(
        images= image_3d,
        threshold=None,
        voxel_size=voxel_size,
        spot_radius=(300, 100,100)
    )

    ufish_spot, enh_im = ufish.predict(image_3d, spots_calling_method= 'local_maxima')
    ufish_spot_cc, enh_im = ufish.predict(image_3d, enh_img=enh_im, spots_calling_method= 'cc_center')
    
    ufish_spot = np.array(list(zip(ufish_spot['axis-0'], ufish_spot['axis-1'], ufish_spot['axis-2'])), dtype=int)
    ufish_spot_cc = np.array(list(zip(ufish_spot_cc['axis-0'], ufish_spot_cc['axis-1'], ufish_spot_cc['axis-2'])), dtype=int)

    viewer = napari.Viewer(title='3D test')
    viewer.add_image(image_3d, )
    viewer.add_image(enh_im, )
    viewer.add_points(bigfish_spot, face_color='blue', size= 5)
    viewer.add_points(ufish_spot,  face_color='green', size= 5)
    viewer.add_points(ufish_spot_cc,  face_color='red', size= 5)

    viewer.show(block=True)
