# eps236-project
EPS236 project source for DH, CL, and SP

## Objective

Set up a 5-box model

## To Do

  1.

## Final Model Parameters

After running a Monte Carlo for 1M iterations, we arrived at our optimal solution:

  * `t.strat`: 6.03 yr
  * `t.hemi.inter`: 1.026
  * `t.hemi.intra`: 0.129
  * `strat.frac`: 0.025

Updated with more resolved value range

* `t.strat`: 5.21 yr
* `t.hemi.inter`: 0.774
* `t.hemi.intra`: 0.150
* `strat.frac`: 0.484

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

## Description of Data

| Filename | Description |
|:--------:|:------------|
| "SF6.emissions.bybox.tbl" | SF6 Emissions by box from 1970-2008 |
