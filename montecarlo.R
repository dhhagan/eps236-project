# Run MonteCarlo Simulations to determine the optimal settings for the model

library(feather)
library(MonteCarlo)

# Initialize the model params
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


# Define the output function to send to the Monte Carlo Simulation
# Essentially, we need to give just the params and it will return the minimized values
monte_carlo_model <- function(tau.stratosphere, tau.hemisphere.inter, tau.hemisphere.intra, 
                              strat.nh.fraction) {
  
  species <- "SF6"
  
  # Run the model
  model.results <- run.model(species, tau.stratosphere, tau.hemisphere.inter, 
                             tau.hemisphere.intra, strat.nh.fraction)
  
  
  # Return the min values for all boxes
  min <- min.cost(model.results, sf6.observations.boxed.annual.means, box.no = NaN)
  #list("min"=min)
  return (min)
}

# Set up the param grid to use
tau.strat.grid <- seq(1, 10, .1)
tau.hemi.inter.grid <- seq(1, 5, .1)
tau.hemi.intra.grid <- seq(0.1, 1, 0.05)
strat.nh.grid <- seq(0.4, 0.6, 0.05)

grid_search = list(
                "tau.stratosphere"=tau.strat.grid,
                "tau.hemisphere.inter"=tau.hemi.inter.grid,
                "tau.hemisphere.intra"=tau.hemi.intra.grid,
                "strat.nh.fraction"=strat.nh.grid)

# Run the Monte Carlo Simulation
mc.results <- MonteCarlo::MonteCarlo(
                            func=monte_carlo_model, 
                            nrep=10, 
                            param_list=grid_search, 
                            max_grid=1000000,
                            time_n_test=TRUE,
                            save_res=TRUE,
                            ncpus = 4)

