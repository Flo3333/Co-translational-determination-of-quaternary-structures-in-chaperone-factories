import pbwrap.quantification.statistical_test as stest
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import SymLogNorm
from itertools import product


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
            samples1 = list(data.iloc[range(0,data_len//2)][key].dropna())
            samples2 = list(data.iloc[range(data_len//2,data_len)][key].dropna())

            pvalue = stest.t_test(
                samples1= samples1,
                samples2= samples2,
                equal_variance= False,
                nan_policy='omit',
            )
            pvalue_df.loc[data.index.droplevel(level_number-1).unique(),["{0}_{1}".format(key,new_key)]] = np.reshape(pvalue, (len(pvalue),1))
    else :
        pvalue_df = pd.Series(index=header)

        for key in measure_keys :
            data = data_grouped[key].apply(list)
            pvalue = stest.t_test(
                samples1= data.iloc[0].dropna(),
                samples2= data.iloc[1].dropna(),
                equal_variance= False,
                nan_policy='omit',
            )[0]

            pvalue_df.loc[key + "_{0}".format(new_key)] = pvalue

    return pvalue_df

def compute_N_treatment_pvalue(
        df : pd.DataFrame,
        group_keys : 'list[str]',
        measure_keys : 'list[str]',
        new_key = '',
        significance = 0.01

) :
    """
    Computes p-values for data with N (N>2) samples and unequal variance.

    The p-value is computed from succesive Alexander-govern (ANOVA for unequal variance samples) test and Games_howell test (Tukey honnestly significant difference test for unequal variance samples).
    The kept p-values are the one from GamesHowell test (which test all samples against one another) unless the Alexander-govern (one per group of sample)

    """
    
    data_grouped = df.groupby(group_keys)
    idx = data_grouped.first().index
    level_number = idx.nlevels
    sample_names = idx.get_level_values(level_number-1).unique()
    sample_number = len(sample_names)
    
    #sample_number check : test made for N > 2 samples.
    if sample_number <= 2 : raise ValueError("Incorrect number of sample to try : {0}. Expected >=2 for this test, check group_keys argument is correct or use different test.".format(sample_number))

    groups = idx.droplevel(level_number-1).unique()

    pvalue_alexander_df = pd.DataFrame(index= idx.droplevel(-1).unique())
    columns_multi_index = product(measure_keys,sample_names)
    pvalue_games_howell_df = pd.DataFrame(index= idx, columns=pd.MultiIndex.from_tuples(columns_multi_index), dtype= float)
    if level_number == 1 : 
        for key in measure_keys :
            data = data_grouped[key].apply(list).sort_index(level=1)
            pvalue = stest.alexandergovern(list(data))
            pvalue_alexander_df[key + "_alexandergovern"] = pvalue
            pvalue = stest.games_howell(list(data))
            pvalue_games_howell_df[key] = pvalue

        pvalue_alexander_df.columns = pvalue_alexander_df.columns.str.cat([new_key]*len(pvalue_alexander_df.columns))
    
    else :
        for key in measure_keys :
            data = data_grouped[key].apply(list)

            #AlexanderGovern            
            grouped_data = data.reset_index(drop=False).groupby(group_keys[:-1])[key].apply(list)
            pvalue =grouped_data.apply(stest.alexandergovern)
            pvalue_alexander_df[key] = pvalue
            
            #GamesHowell
            for group in pvalue_alexander_df.index :
                sub_data = data.loc[group]
                pvalue = stest.games_howell(list(sub_data))
                pvalue_games_howell_df.loc[group, (key,list(sub_data.index))] = pvalue

        for group in pvalue_alexander_df.index :
            fail_col = pvalue_alexander_df.columns[(pvalue_alexander_df.loc[group] > significance) & ~(pvalue_alexander_df.loc[group].isna())]
            if len(fail_col) == 0 : continue
            pvalue_games_howell_df.loc[group,fail_col] = float(pvalue_alexander_df.loc[group, fail_col])

    col_lvl0= pvalue_games_howell_df.columns.get_level_values(0).unique()
    pvalue_games_howell_df.columns = pd.MultiIndex.from_tuples(product(
        col_lvl0.str.cat([new_key]*len(col_lvl0)),
        sample_names
    ))

    pvalue_games_howell_df = pvalue_games_howell_df.sort_index(axis=0).sort_index(axis=1)

    return pvalue_games_howell_df


def pairwise_test_tile(
        ax : plt.Axes,
        pvalues: pd.DataFrame,
        measure,
        significance,
) :
    
    """
    Parameters
    ----------

        pvalues : pd.DataFrame
            Expected multi-index for axis 0  with samples on last level.  
            *EG (cell_line, rna, treatment) for a pvalue test on `treatments`.* 
            Expected multi index for axis 1 with samples on last level. Use or see `CustomPandasFramework.pvalues.compute_N_treatment_pvalue` to obtain such output.
        
        measure : str
            key on which test was performed.
        significance : float
            significance of the test, used for color scale.

    """


    if isinstance(pvalues.columns, pd.MultiIndex) :
        if measure not in pvalues.columns.get_level_values(0) : raise KeyError("{0} key not found in first level of column index. If your key is in another level please use a cross section to select correct p-values first.".format(measure))
        pvalues = pvalues[measure]
    if not pvalues.columns.equals(pvalues.index.get_level_values(-1)) : raise IndexError('pvalues dataframe : Last level of axis 0 must match axis1 index or axis index last level.\nindex0\n{0}\nindex1\n{1}'.format(pvalues.index.get_level_values((-1)), pvalues.columns))
    if len(pvalues.index.get_level_values(0)) != len(pvalues.index.get_level_values(-1).unique()) : raise ValueError("Not supported for  'multi-tiles'")
    pvalues = pvalues.sort_index(axis=0).sort_index(axis=1)

    #Colormesh
    datamesh = pvalues.to_numpy()
    colormesh = ax.pcolormesh(
        datamesh,
        edgecolors= 'black',
        norm= SymLogNorm(vmin= 10e-19, vmax= 1, linthresh= significance, linscale= -np.log10(significance)),
        cmap= 'bwr',
        )

    #Colorbar
    cbar = plt.colorbar(colormesh, ax=ax, location= 'right', ticks = [1, significance, 1e-1, 1e-3, 1e-18])
    
    #ticks
    ticks = np.arange(0.5,0.5 + len(datamesh),1)
    ticks_label = list(pvalues.columns)
    ax.set_xticks(ticks, ticks_label)
    ax.set_yticks(ticks, ticks_label)
    ax.set_title("Pairwise p-values({0})".format(measure))

    return ax