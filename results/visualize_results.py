"""Visualize the results using Python
"""

import feather
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from datetime import datetime

sns.set("talk", style='ticks', palette='dark', font_scale=1.5)

sf6 = feather.read_dataframe("sf6_emissions.feather")
mod = feather.read_dataframe("min_model_results_mc.feather")
res = feather.read_dataframe("iterative_results.feather")
sw = feather.read_dataframe("sw-results.feather")

# Let's make a confusion matrix to show which variables actually matter!
confmat = res[['t.strat', 't.hemi.inter', 't.hemi.intra', 'strat.frac', 'min.box.all']]

corr = confmat.corr()

sns.heatmap(corr)

x = confmat['t.strat']
y = confmat['t.hemi.inter']

X, Y = np.meshgrid(x, y)
Z = confmat['min.box.all']

plt.figure()
plt.pcolormesh(confmat[['t.strat', 't.hemi.inter', 'min.box.all']].values, cmap='RdBu')
plt.colorbar()

res.info()

res.describe()


# Grab results where t.strat is 2

del mod["Group.1"]
del sf6["Group.1"]
del sw['Group.1']

# Convert the year to an actual timestamp
mod['year'] = mod['year'].apply(lambda x: datetime(math.floor(x), 1, 1))
sf6['year'] = sf6['year'].apply(lambda x: datetime(math.floor(x), 1, 1))
sw['year'] = sw['year'].apply(lambda x: datetime(math.floor(x), 1, 1))

# Convert the index to a datetime
mod.set_index("year", inplace=True)
sf6.set_index("year", inplace=True)
sw.set_index("year", inplace=True)

sf6.head()

sw.head()

mod.head()

cp = sns.color_palette()

fig, ax = plt.subplots(1, figsize=(14, 10))

# Plot the modeled results
ax.plot(mod['box.1'], lw=3, label='box.1', c=cp[0])
ax.plot(mod['box.2'], lw=3, label='box.2', c=cp[1])
ax.plot(mod['box.3'], lw=3, label='box.3', c=cp[2])
ax.plot(mod['box.4'], lw=3, label='box.4', c=cp[3])
ax.plot(mod['box.5'], lw=4, label='box.5', c=cp[4])

# Plot the known values
ax.plot(sf6['SF6.box.1'], '--', c=cp[0])
ax.plot(sf6['SF6.box.2'], '--', c=cp[1])
ax.plot(sf6['SF6.box.3'], '--', c=cp[2])
ax.plot(sf6['SF6.box.4'], '--', c=cp[3])

ax.plot(sw['box.1'], 'o', c=cp[0], label="SW.box.1")
#ax.plot(sw['box.2'], 'o', c=cp[1], label="SW.box.2")
#ax.plot(sw['box.3'], 'o', c=cp[2], label="SW.box.3")
#ax.plot(sw['box.4'], 'o', c=cp[3], label="SW.box.4")

ax.legend(loc='best')
ax.set_ylabel("SF6 (ppt)")
ax.set_xlabel("")
ax.set_title("Modeled SF6 Results", y=1.05)
sns.despine(offset=5)

# Timeseries of the residulas
fig, ax = plt.subplots(1, figsize=(14, 10))

# Plot the modeled results
ax.plot(sf6['SF6.box.1'] - mod['box.1'], lw=3, label='box.1', c=cp[0])
ax.plot(sf6['SF6.box.2'] - mod['box.2'], lw=3, label='box.2', c=cp[1])
ax.plot(sf6['SF6.box.3'] - mod['box.3'], lw=3, label='box.3', c=cp[2])
ax.plot(sf6['SF6.box.4'] - mod['box.4'], lw=3, label='box.4', c=cp[3])
#ax.plot(mod['box.5'], lw=4, label='box.5', c=cp[4])

ax.legend(loc='best')
ax.set_ylabel("SF6 Residuals (ppt)")
ax.set_xlabel("")
ax.set_title("Modeled SF6 Results", y=1.05)
sns.despine(offset=5)




# Histogram Distplot showing distribution of residuals
plt.figure()
ax = sns.distplot(a=(sf6["SF6.box.1"] - mod["box.1"]), hist=False, label="Box.1")
ax = sns.distplot(a=(sf6["SF6.box.2"] - mod["box.2"]), hist=False, label="Box.2", ax=ax)
ax = sns.distplot(a=(sf6["SF6.box.3"] - mod["box.3"]), hist=False, label="Box.3", ax=ax)
ax = sns.distplot(a=(sf6["SF6.box.4"] - mod["box.4"]), hist=False, label="Box.4", ax=ax)
