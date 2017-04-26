"""Visualize the results using Python
"""

# Make imports
import feather
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from datetime import datetime

# Set default seaborn parameters
sns.set("talk", style='ticks', palette='dark', font_scale=1.5, color_codes=True)

# Define the color palette as a list of RGB values
cp = sns.color_palette()

# Simple function to plot the timeseries
def tsplot(modeled, observed, max_y_ticks=None, alpha=1, title=''):
    """Plot the observed sf6 timeseries against modeled results.
    """
    fig, ax = plt.subplots(1, figsize=(14, 10))

    # Plot the modeled results
    ax.plot(modeled['box.1'], lw=4, label='North-High Lat.', c=cp[0], alpha=alpha)
    ax.plot(modeled['box.2'], lw=4, label='North Tropics', c=cp[1], alpha=alpha)
    ax.plot(modeled['box.3'], lw=4, label='South Tropics', c=cp[2], alpha=alpha)
    ax.plot(modeled['box.4'], lw=4, label='South-High Lat.', c=cp[3], alpha=alpha)
    ax.plot(modeled['box.5'], lw=4, label='Stratosphere', c=cp[4], alpha=alpha)

    ax.legend(loc='best')

    ax.plot(observed['SF6.box.1'], 'o', c=cp[0])
    ax.plot(observed['SF6.box.2'], 'o', c=cp[1])
    ax.plot(observed['SF6.box.3'], 'o', c=cp[2])
    ax.plot(observed['SF6.box.4'], 'o', c=cp[3])

    if max_y_ticks is not None:
        yloc = plt.MaxNLocator(max_y_ticks)
        ax.yaxis.set_major_locator(yloc)

    ax.set_ylabel("SF6 (ppt)", fontsize=24)
    ax.set_xlabel("")
    ax.set_title(title, y=1.05, fontsize=36)

    sns.despine(offset=5)

    return ax

# Build an array of each individual pulsed response to plot
pulses = []
pulses.append(("North-High Lat.", "results/pulse.box1.feather"))
pulses.append(("North Tropics", "results/pulse.box2.feather"))
pulses.append(("South Tropics", "results/pulse.box3.feather"))
pulses.append(("South-High Lat.", "results/pulse.box4.feather"))

for pulse in pulses:
    df = feather.read_dataframe(pulse[1])

    # Plot the data
    fig, ax = plt.subplots(1, figsize=(14, 10))

    ax.plot(df['time'], df['box.1'], c=cp[0], lw=6, label="North-High Lat.")
    ax.plot(df['time'], df['box.2'], c=cp[1], lw=6, label="North Tropics")
    ax.plot(df['time'], df['box.3'], c=cp[2], lw=6, label="South Tropics")
    ax.plot(df['time'], df['box.4'], c=cp[3], lw=6, label="South-High Lat.")
    ax.plot(df['time'], df['box.5'], c=cp[4], lw=6, label="Stratosphere")

    sns.despine(offset=5)

    ax.set_ylabel("SF6 (ppt)", fontsize=24)
    ax.legend(loc='best')

    yloc = plt.MaxNLocator(3)
    ax.yaxis.set_major_locator(yloc)


# Next, we want to plot the grid-searched results
# Set up a dataframe with the SF6 Actual Measurements
sf6 = feather.read_dataframe("results/sf6_emissions.feather")
gs = feather.read_dataframe("results/gridsearch_optimal.feather")
sw = feather.read_dataframe("results/sw-results.feather")

# Convert the year to an actual datetime
sf6['year'] = sf6['year'].apply(lambda x: datetime(math.floor(x), 1, 1))
gs['year'] = gs['year'].apply(lambda x: datetime(math.floor(x), 1, 1))
sw['year'] = sw['year'].apply(lambda x: datetime(math.floor(x), 1, 1))

# Set the datetime as the index
sf6.set_index("year", inplace=True)
gs.set_index("year", inplace=True)
sw.set_index("year", inplace=True)

# Make the Grid Search Plot
ax = tsplot(gs, sf6, max_y_ticks=3)

# Make the Grid Search Plot
ax = tsplot(sw, sf6, max_y_ticks=3)

# Make the Monte Carlo Results
mod = feather.read_dataframe("results/mc_results_1M_final.feather")
res = feather.read_dataframe("results/mc_results_1M_iters.feather")

# Convert the year to an actual timestamp
mod['year'] = mod['year'].apply(lambda x: datetime(math.floor(x), 1, 1))

# Convert the index to a datetime
mod.set_index("year", inplace=True)

ax = tsplot(mod, sf6)

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
