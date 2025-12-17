# Notes on FishSeq developement roadmap

## Images info :

*DAPI z-stack* -> (c,z,y,x)
channel :   0 dapi,
            1 beads

*DAPI FISH* -> (c,z,y,x)
channel :   0 color1,
            1 color2,
            2 beads

*Note : One image-stack out of 2 is a wash acquisition*

*Beads-Field* -> (c,z,t,x)
channel :   0 wavelength color1,
            0 wavelength color2,
            0 wavelength beads,


## Pipeline structure

### 1. Segmentation *(FishSeq_pipeline_segmentation.py)*
--> translation correction between **dapi-stack** beads and **FISH-stack00** beads (reference).
--> nucleus segmentation on corrected dapi
--> cytoplasm segmentation on beads channel + corrected dapi?
--> output as **.npyz**

*Note : for dapi correction we ignore the chromatical shift*

### 2. Detection *(FishSeq_pipeline_detection.py)*
--> Perform regular spot detection on both color for all image-stack. (each location, each cycle) (loop + multi-thread)
--> output as np.arrays to optimize next steps. Or maybe store np.arrays in df to easy access.

### 3. Shift correction *(FishSeq_pipeline_quantification.py)*
--> Perform shift correction on spot coordinates from **FISH-stack** beads channel (2) with stack 00 as reference. **(Linear correction)**
--> Perform chromatical correction on spot coordinates from **Beads-Field** channels reference wavelength to be discussed. **(Polynomial correction)**
--> *(opt) output*

### 4. Quantification
--> colocalisation quantification. *idea* : use 4D np array to compute distance map with inf distance in 4th dimension ?
--> cell quantification
--> output