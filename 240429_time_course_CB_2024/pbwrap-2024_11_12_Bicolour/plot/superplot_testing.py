import pandas as pd
import numpy as np
from pbwrap.plot import distribution_super_plot
import matplotlib.pyplot as plt

ARN_NUMBER = 4
TREATMENT_NUMBER = 3
EXPERIMENT_NUMBER = 3
TOTAL_POINT = 100000
SEED = 1


arn_list = np.array(["ARN {0}".format(i) for i in range(ARN_NUMBER)])
treatment_list = np.array(["Treatment {0}".format(i) for i in range(TREATMENT_NUMBER)])
experiment_list = np.array(["exp {0}".format(i) for i in range(EXPERIMENT_NUMBER)])

random = np.random.default_rng(seed=SEED)


Df = pd.DataFrame({
    'id' : np.arange(TOTAL_POINT),
    'arn' : arn_list[random.integers(0,ARN_NUMBER,size=TOTAL_POINT)],
    'treatment' : treatment_list[random.integers(0,TREATMENT_NUMBER,size=TOTAL_POINT)],
    'exp' : experiment_list[random.integers(0,EXPERIMENT_NUMBER,size=TOTAL_POINT)],
    'measure' : random.random(size= TOTAL_POINT),
})

index =Df.loc[Df['exp'] == 'exp 2'].index
Df.loc[index, ['measure']] += 1.5
index =Df.loc[Df['exp'] != 'exp 2']['measure'].index
Df.loc[index, ['measure']] -= 1.5

data = Df.groupby(['arn', 'treatment','exp'])['measure'].apply(list)

fig = plt.figure(figsize=(10,10))
ax = fig.gca()
ax = distribution_super_plot(
    data=data,
    ax=ax,
    xlabel= 'ARN',
    ylabel='measure'
)
plt.show()