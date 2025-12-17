
#testing random spots from mask
import numpy as np
import CustomPandasFramework.Fish_seq.quantification as quant


if False :
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


#Testing _harmonise_sym_data
if True :
    from CustomPandasFramework.Fish_seq.plot import _harmonize_sym_data
    import pandas as pd
    index1 = pd.Index(['pomme','poire','pêche'])
    index2 = pd.Index(['poire','pomme','pêche'])
    index3 = pd.Index(['poire','pomme','pêche'])
    index4 = pd.Index(['poire','pêche','mela'])

    print("index1 == index2", index1.equals(index2))
    print("index2 == index3", index2.equals(index3))

    df1 = pd.DataFrame({'1' : ['verte', 'jaune', 'rose']},index=index1)
    df4 = pd.DataFrame({'4' : ['verte', 'blanche', 'rossa']},index=index4)

    print("df1 : ", df1)
    print("df4 : ", df4)

    df1, df4 = _harmonize_sym_data(df1, df4)

    print("df : ", df1)
    print("df : ", df4)
    print()
    print(df1.loc[('pomme')])
