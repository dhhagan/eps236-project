# Initialize the Model and define utility functions
# Inititalize the Model by reading in all experimental data

# Import libraries to use
library(dplyr)
library(feather)
library(progress)

# Source the model initialization file (import all data, etc.)
source("model_initialize.r")

# Import our run.model function
source("utils.r")

# Define the species of interest (NOTE: SF6 is the only one that works...too many things are hard-coded)
species = "SF6"

# Bring in our observations for SF6
# Again, we're lazy and R is dumb. Normally, we would never want to hardcode the index, but it's easier 
# since working with real datetimes in R is unbelievably complex for it being a "good" language to do 
# timeseries analysis in...
idx.low <- 218
idx.high <- 383

# Grab the SF6 observations from 1995.125-2008.875
sf6.observations <- data.frame(
                      ghg.observations[,paste(species, "box", 1:4, sep=".")], 
                      row.names=ghg.observations[,"Year"])[idx.low:idx.high,]

# Get the year from the index
sf6.observations <- cbind(sf6.observations, year=floor(as.numeric(row.names.data.frame(sf6.observations))))

# Average the observations by year
sf6.observations.boxed.annual.means <- aggregate(sf6.observations, list(sf6.observations$year), mean)

# stevewofsy results (sw.res) using the default params listed on the slides
# tau.stratosphere <- 4
# tau.hemisphere.inter <- 2
# tau.hemisphere.intra <- 0.5
# strat.nh.fraction <- 0.55
sw.res <- run.model(species, 4, 2, 0.5, 0.55)

# Perform a grid-search to run all possible combinations as defined by the following grids
tau.stratosphere.list = seq(1,20,0.5)
tau.hemisphere.inter.list = seq(1,5,0.5)
tau.hemisphere.intra.list = seq(0.1, 1, 0.1)
strat.nh.fraction.list = seq(0.4, 0.6, 0.05)

# Define an empty matrix that we will use to store results
res.matrix = matrix(ncol=9, nrow=0)

# determine the number of total runs to initiate the sweet progress bar
n_runs = length(tau.stratosphere.list) * length(tau.hemisphere.inter.list) * length(tau.hemisphere.intra.list) * 
            length(strat.nh.fraction.list)

# Initialize the progress bar so you know how many years you will have to wait...
progress.bar <- progress::progress_bar$new(total=n_runs, format=" running the codez [:bar] :percent eta: :eta")

# Perform the grid-search
for (t.strat in tau.stratosphere.list) {
  for (t.hemi.inter in tau.hemisphere.inter.list) {
    for (t.hemi.intra in tau.hemisphere.intra.list) {
      for (strat.frac in strat.nh.fraction.list) {
        # Run the model
        res.ind <- run.model(species, tau.stratosphere = t.strat, tau.hemisphere.inter = t.hemi.inter, 
                  tau.hemisphere.intra = t.hemi.intra, strat.nh.fraction = strat.frac)
        
        # Return the residuals for each box
        # NOTE: Ended up only using the min.box.all results
        min.box.1 <- min.cost(res.ind, sf6.observations.boxed.annual.means, box.no = 1)
        min.box.2 <- min.cost(res.ind, sf6.observations.boxed.annual.means, box.no = 2)
        min.box.3 <- min.cost(res.ind, sf6.observations.boxed.annual.means, box.no = 3)
        min.box.4 <- min.cost(res.ind, sf6.observations.boxed.annual.means, box.no = 4)
        min.box.all <- min.cost(res.ind, sf6.observations.boxed.annual.means, box.no = NaN)
        
        res.matrix <- rbind(res.matrix, c(t.strat, t.hemi.inter, t.hemi.intra, strat.frac, min.box.1, min.box.2, 
                                         min.box.3, min.box.4, min.box.all))
        
        # Update the progress bar
        progress.bar$tick()
      }
    }
  }
}

# Store the results as a data.frame rather than a matrix
res.df <- data.frame(res.matrix)

# rename the columns in the data.frame
colnames(res.df) <- c("t.strat", "t.hemi.inter", "t.hemi.intra", "strat.frac", "min.box.1", 
                      "min.box.2", "min.box.3", "min.box.4", "min.box.all")

# Find the row with the lowest min.box.all value
min.params <- arrange(res.df, min.box.all)[1,]

# Obtain the actual model results with the minimum params
min.results <- run.model(species, min.params$t.strat, min.params$t.hemi.inter, min.params$t.hemi.intra,
                         min.params$strat.frac)

# Save the data in feather format so we can plot them later in Python
feather::write_feather(res.df, "results/iterative_gridsearch_results.feather")
feather::write_feather(min.results, "results/model_results_minimized.feather")
feather::write_feather(sf6.observations.boxed.annual.means, "results/sf6_emissions.feather")

