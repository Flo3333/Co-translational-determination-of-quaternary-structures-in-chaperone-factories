
#testing random spots from mask
import numpy as np
import CustomPandasFramework.FishBicolor.quantification as quant

seed= None
shape = (2,10,10)
spot_number = 5

random_gen = np.random.default_rng()
mask = random_gen.integers(0,2, size=shape, dtype= bool)
print("MASK\n", mask)

spots = quant._random_spot_list_from_mask(shape, spot_number, mask)
print(spots)

Z,Y,X = list(zip(*spots))


print(all(mask[Z,Y,X]))
print((mask[Z,Y,X]))