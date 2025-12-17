import pbwrap.quantification.statistical_test as stest
import pandas as pd
import numpy as np


def compute_Ttest(
        df : pd.DataFrame,
        group_keys : 'list[str]',
        measure_keys : 'list[str]',
        new_key = 'pvalue',
) :
    data_grouped = df.groupby(group_keys)
    idx = data_grouped.first().index
    level_number = idx.nlevels
    sample_names = idx.get_level_values(level_number-1).unique()
    sample_number = len(sample_names)
    header = [measure + "_{0}".format(new_key) for measure in measure_keys]
    
    #sample_number check : T_test are valid only to 1vs1 
    if sample_number != 2 : raise ValueError("Incorrect number of sample to try : {0}. Expected 2 for t-test, check group_keys argument is correct or use different test.".format(sample_number))

    if level_number > 1 :
        groups = idx.droplevel(level_number-1).unique()
        pvalue_df = pd.DataFrame(index= groups,columns=header)

        for key in measure_keys :
            
            data = data_grouped[key].apply(list).sort_index(level=1)
            data = data.reset_index(level=level_number-1)
            sample1_idx = data[data[group_keys[-1]] == sample_names[0]].index
            sample2_idx = data[data[group_keys[-1]] == sample_names[1]].index
            data = data.loc[data.index.isin(sample1_idx) & data.index.isin(sample2_idx)].set_index(group_keys[-1], append=True)
            data = data.sort_index(level=level_number-1)

            data_len = int(len(data))
            assert data_len %2 == 0
            samples1 = list(data.iloc[range(0,data_len//2)][key])
            samples2 = list(data.iloc[range(data_len//2,data_len)][key])

            pvalue = stest.t_test(
                samples1= samples1,
                samples2= samples2,
                equal_variance= False
            )
            pvalue_df.loc[data.index.droplevel(level_number-1).unique(),["{0}_{1}".format(key,new_key)]] = np.reshape(pvalue, (len(pvalue),1))
    else :
        pvalue_df = pd.Series(index=header)

        for key in measure_keys :
            data = data_grouped[key].apply(list)
            pvalue = stest.t_test(
                samples1= data.iloc[0],
                samples2= data.iloc[1],
                equal_variance= False
            )[0]

            pvalue_df.loc[key + "_{0}".format(new_key)] = pvalue

    return pvalue_df