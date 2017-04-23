"""Visualize the results using Python
"""

import feather
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from datetime import datetime

sns.set("talk", style='ticks', palette='dark', font_scale=1.5, color_codes=True)

sf6 = feather.read_dataframe("sf6_emissions.feather")
mod = feather.read_dataframe("mc_results_1M_final.feather")
res = feather.read_dataframe("mc_results_1M_iters.feather")
sw = feather.read_dataframe("sw-results.feather")

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


cp = sns.color_palette()

fig, ax = plt.subplots(1, figsize=(14, 10))

# Plot the modeled results
ax.plot(mod['box.1'], lw=4, label='box.1', c=cp[0])
ax.plot(mod['box.2'], lw=4, label='box.2', c=cp[1])
ax.plot(mod['box.3'], lw=4, label='box.3', c=cp[2])
ax.plot(mod['box.4'], lw=4, label='box.4', c=cp[3])
ax.plot(mod['box.5'], lw=4, label='box.5', c=cp[4])

# Plot the known values
ax.plot(sf6['SF6.box.1'], '-.', c=cp[0])
ax.plot(sf6['SF6.box.2'], '-.', c=cp[1])
ax.plot(sf6['SF6.box.3'], '-.', c=cp[2])
ax.plot(sf6['SF6.box.4'], '-.', c=cp[3])

ax.legend(loc='best')
ax.set_ylabel("SF6 (ppt)")
ax.set_xlabel("")
ax.set_title("Monte Carlo SF6 Results", y=1.05)
sns.despine(offset=5)



# Timeseries of the residulas
fig, ax = plt.subplots(1, figsize=(14, 10))

# Plot the modeled results
ax.plot(sf6['SF6.box.1'] - mod['box.1'], 'o-', lw=3, label='box.1', c=cp[0])
ax.plot(sf6['SF6.box.2'] - mod['box.2'], 'o-', lw=3, label='box.2', c=cp[1])
ax.plot(sf6['SF6.box.3'] - mod['box.3'], 'o-', lw=3, label='box.3', c=cp[2])
ax.plot(sf6['SF6.box.4'] - mod['box.4'], 'o-', lw=3, label='box.4', c=cp[3])

yloc = plt.MaxNLocator(3)
ax.yaxis.set_major_locator(yloc)
plt.ylim([-0.1, 0.1])

ax.legend(loc='best')
ax.set_ylabel("SF6 Residuals (ppt)")
ax.set_xlabel("")
ax.set_title("Monte Carlo SF6 Residuals", y=1.05)
sns.despine(offset=5)

# Plot the KDE of various Parameters
best = res.sort_values('min.box.all').head(1000)
best.tail()


with sns.axes_style("white"):
    f, ax = plt.subplots(2, 2, figsize=(14, 10))

    sns.despine(left=True)

    # Plot the t.strat distribution
    sns.distplot(best['t.strat'], kde=True, hist=False, color='b', ax=ax[0, 0], kde_kws=dict(shade=True))

    # Plot the t.hemi.inter distribution
    sns.distplot(best['t.hemi.inter'], kde=True, hist=False, color='g', ax=ax[0, 1], kde_kws=dict(shade=True))

    # Plot the t.hemi.inter distribution
    sns.distplot(best['t.hemi.intra'], kde=True, hist=False, color='r', ax=ax[1, 0], kde_kws=dict(shade=True))

    # Plot the t.hemi.inter distribution
    sns.distplot(best['strat.frac'], kde=True, hist=False, color='m', ax=ax[1, 1], kde_kws=dict(shade=True))

    ax[0,0].set_ylim([0, ax[0,0].get_ylim()[-1]])
    ax[0,1].set_ylim([0, ax[0,1].get_ylim()[-1]])
    ax[1,0].set_ylim([0, ax[1,0].get_ylim()[-1]])
    ax[1,1].set_ylim([0, ax[1,1].get_ylim()[-1]])

    plt.setp(ax, yticks=[])

    plt.tight_layout()
