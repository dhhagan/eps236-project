# Our attempt at a Monte Carlo Simulation
# Run MonteCarlo Simulations to determine the optimal settings for the model

library(dplyr)
library(feather)
library(progress)

# Initialize the model params
source("model_initialize.r")

# Import our run.model function that takes various params and returns the results (box-by-box concentration timeseries)
source("utils.r")

# Calculate params for SF6
species = "SF6"

# Define the number of iterations to use in the Monte Carlo simulation
num.iterations <- 100000

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

# Set up our version of a monte carlo where we choose values from a random sequence thousands of time.


# Set up the progress bar
progress.bar <- progress::progress_bar$new(total=num.iterations, format=" running the codez [:bar] :percent eta: :eta")

# Set up an empty matrix for results
res.matrix = matrix(ncol=5, nrow=0)

# Min and Max values
t.strat.min <- 1
t.strat.max <- 15
t.hemi.inter.min <- 0.01
t.hemi.inter.max <- 2.0
t.hemi.intra.min <- 0.01
t.hemi.intra.max <- 0.5
strat.frac.min <- 0.2
strat.frac.max <- 0.8

for (i in 1:num.iterations) {
  # Get random values between min and max for each param
  t.strat <- runif(1, min=t.strat.min, max=t.strat.max)
  t.hemi.inter <- runif(1, min=t.hemi.inter.min, max=t.hemi.inter.max)
  t.hemi.intra <- runif(1, min=t.hemi.intra.min, max=t.hemi.intra.max)
  strat.frac <- runif(1, min=strat.frac.min, max=strat.frac.max)
  
  # Run the model
  res.ind <- run.model(species, t.strat, t.hemi.inter, t.hemi.intra, strat.frac)
  
  min.box.all <- min.cost(res.ind, sf6.observations.boxed.annual.means, box.no = NaN)
  
  res.matrix <- rbind(res.matrix, c(t.strat, t.hemi.inter, t.hemi.intra, strat.frac, min.box.all))
  
  # Update the progress bar
  progress.bar$tick()
}

res.df <- data.frame(res.matrix)
colnames(res.df) <- c("t.strat", "t.hemi.inter", "t.hemi.intra", "strat.frac", "min.box.all")

# Find the row with the lowest min.box.all value
min.params <- arrange(res.df, min.box.all)[1,]

# Obtain the actual model results with the minimum params and write to a feather file
min.results <- run.model(species, min.params$t.strat, min.params$t.hemi.inter, min.params$t.hemi.intra,
                         min.params$strat.frac)

# Write results to file
feather::write_feather(res.df, "results/mc_results_by_iter.feather")
feather::write_feather(min.results, "results/mc_results_final.feather")
