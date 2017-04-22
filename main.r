# Initialize the Model and define utility functions
# Inititalize the Model by reading in all experimental data ()
# Set constants:
#
#
#
#

library(dplyr)
library(feather)
library(progress)

source("model_initialize.r")

# Import our run.model function that takes various params and returns the results (box-by-box concentration timeseries)
source("utils.r")


# Calculate params for SF6
species = "SF6"

# Bring in our observations for SF6
idx.low <- 218
idx.high <- 383

# Observations go from 1995.125-2008.875
sf6.observations <- data.frame(
                      ghg.observations[,paste(species, "box", 1:4, sep=".")], 
                      row.names=ghg.observations[,"Year"])[idx.low:idx.high,]

# Get the year from the index, and then average the observations by year
sf6.observations <- cbind(sf6.observations, year=floor(as.numeric(row.names.data.frame(sf6.observations))))

sf6.observations.boxed.annual.means <- aggregate(sf6.observations, list(sf6.observations$year), mean)

# Define default params
#tau.stratosphere <- 4
#tau.hemisphere.inter <- 2
#tau.hemisphere.intra <- 0.5
#strat.nh.fraction <- 0.55

sw.res <- run.model(species, 4, 2, 0.5, 0.55)

# Iterate over all combinations
tau.stratosphere.list = seq(1,20,0.5)
tau.hemisphere.inter.list = seq(1,5,0.5)
tau.hemisphere.intra.list = seq(0.1, 1, 0.1)
strat.nh.fraction.list = seq(0.4, 0.6, 0.05)

res.matrix = matrix(ncol=9, nrow=0)

n_runs = length(tau.stratosphere.list) * length(tau.hemisphere.inter.list) * length(tau.hemisphere.intra.list) * 
            length(strat.nh.fraction.list)

progress.bar <- progress::progress_bar$new(total=n_runs, format=" running the codez [:bar] :percent eta: :eta")

for (t.strat in tau.stratosphere.list) {
  for (t.hemi.inter in tau.hemisphere.inter.list) {
    for (t.hemi.intra in tau.hemisphere.intra.list) {
      for (strat.frac in strat.nh.fraction.list) {
        # Run the model
        res.ind <- run.model(species, tau.stratosphere = t.strat, tau.hemisphere.inter = t.hemi.inter, 
                  tau.hemisphere.intra = t.hemi.intra, strat.nh.fraction = strat.frac)
        
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


res.df <- data.frame(res.matrix)
colnames(res.df) <- c("t.strat", "t.hemi.inter", "t.hemi.intra", "strat.frac", "min.box.1", 
                      "min.box.2", "min.box.3", "min.box.4", "min.box.all")

# Find the row with the lowest min.box.all value
min.params <- arrange(res.df, min.box.all)[1,]

# Obtain the actual model results with the minimum params and write to a feather file
min.results <- run.model(species, min.params$t.strat, min.params$t.hemi.inter, min.params$t.hemi.intra,
                         min.params$strat.frac)

# Save the data in feather format
feather::write_feather(res.df, "results/iterative_results.feather")
feather::write_feather(min.results, "results/model_results_min2.feather")
feather::write_feather(sf6.observations.boxed.annual.means, "results/sf6_emissions.feather")

