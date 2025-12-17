import numpy as np
import pandas as pd
import CustomPandasFramework.pvalues as test


"""
Sample creation - Normal distribution
"""
SEED = 1
N = 10000
random_gen = np.random.default_rng(seed=SEED)

#SAMPLE1
sample1_mean = 10.00
sample1_std = 0.1

#SAMPLE2
sample2_mean = 10.00
sample2_std = 0.1

#SAMPLE3
sample3_mean = 30
sample3_std = 0.1

#SAMPLE4
sample4_mean = 40
sample4_std = 1


sample1 = random_gen.normal(loc=sample1_mean, scale=sample1_std, size=N)
sample2 = random_gen.normal(loc=sample2_mean, scale=sample2_std, size=N)
sample3 = random_gen.normal(loc=sample3_mean, scale=sample3_std, size=N)
sample4 = random_gen.normal(loc=sample4_mean, scale=sample4_std, size=N)

# gender
sample_number = 2
gender = random_gen.choice(['male', 'female'], sample_number * N)

# DATAFRAME

df = pd.DataFrame({
    "id" : np.arange(N*2),
    "size" : np.concatenate([sample1 , sample2]),
    "weight" : np.concatenate([sample3 , sample4]),
    "population" : [1] * N + [2] * N,
    "gender" : gender
})

# RESULTS


print("DATA\n",df)
pvalue = test.compute_Ttest(
    df = df,
    group_keys= ["gender","population"],
    measure_keys= ["size","weight"],
)


print(pvalue)