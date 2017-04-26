# EPS236 Final Class Project
EPS236 project source for David Hagan (MIT), Chris Lim (MIT), and Sidhant Pai (MIT)

## Objective

Set up a 5-box model and minimize the residuals between experimental results from EDGAR and that of our model output.

## Requirements

The following R packages are required:

  * feather
  * MonteCarlo
  * dplyr
  * progress
  
The following Python packages are required to make the plots:

  * seaborn
  * matplotlib
  * pandas
  * feather

## Getting Started

We have attempted to minimize the error in the model in 3 different ways:
 
  1. Iterative grid-search
  2. Monte Carlo using the MonteCarlo package
  3. Monte Carlo as written by us
  
Each method is contained within it's own file. Below, we show how to run each file.

### Iterative Grid Search

The run file for this method is `iterative_gridsearch.R`. To run the file, you should be able to just run from within RStudio or from the R command line. The output is three files, all in feather format.

  1. **results/iterative_gridsearch_results.feather**
  
      Contains results for each individual run in a table containing columns for each of the optimized   parameters, as well as the sum of the squared residuals.
      
  2. **results/model_results_minimized.feather**
  
      Contains the final model results from the run with the minimized residuals as a timeseries.
      
  3. **results/sf6_emissions.feather**
  
      Contains the observed SF6 emissions as a timeseries (for plotting purposes)

### Monte Carlo using the MonteCarlo package

We found a MonteCarlo package, however it appears that it works for as an iterative grid-search than a true Monte Carlo simulation. We tried parsing the source, but it was written pretty horribly...Regardless, the results output the files from running `montecarlo.R`:

Not updated...

### Monte Carlo from Scratch

We weren't totally sure the MonteCarlo package was working the way we thought it should, so we took matters into our own hands and wrote the file `monte_carlo_mit.R`. This version is a Monte Carlo algorithm we wrote to pick from input distributions for each parameter. You can change the number of iterations through the `num.iterations` variable. There are two output files:

  1. **results/mc_results_by_iter.feather**
  
      Contains the run-by-run results of the Monte Carlo simulation with parameters and results.
  
  2. **results/mc_results_final.feather**
  
      Contains the final model results as found by optimizing the parameters.


## Final Model Parameters

After running a Monte Carlo for 1M iterations, we arrived at our optimal solution:

| Variable | Result |
|:--------:|:------:|
| `t.strat` | 5.21 yr |
| `t.hemi.inter` | 0.774 |
| `t.hemi.intra` | 0.150 |
| `strat.frac` | 0.484 |


## Variable Definitions

| Variable Name | Definition |
|:-------------:|:-----------|
|`tau.global.lifetimes.years`| Global mean photochemical lifetime for several species|
|`ghg.observations`| NASA GHG observations for several species by box for many years (>1977)|
|`loss.monthly`|Monthly loss rate for several compounds by box per month|
|`sf6.sources`|Yearly SF6 emissions by box from EDGAR|
|`time.range.sf6.sources`|Vector with the first and last year of EDGAR sources |
|`mass.global`|Total mass of the atmosphere in units of kg|
|`mass.stratosphere`|Total mass of the stratosphere in kg|
|`molecular.wt.all`|Molecular weight of various species in terms of just terms of Carbon and Nitrogen|

## Description of Data Files

| Filename | Description |
|:---------|:------------|
| data/tglobal_all.txt | Table of global lifetimes in years |
| data/ghg.GMD.conc.tbl | NASA Greenhouse Gas emissions |
| data/loss.monthly.tbl | Box-by-box losses as a function of month |
| data/SF6.emissions.bybox.tbl | SF6 Emissions by box from 1970-2008 |

