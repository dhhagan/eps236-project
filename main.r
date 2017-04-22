# Initialize the Model and define utility functions
# Inititalize the Model by reading in all experimental data ()
# Set constants:
#
#
#
#
library(dplyr)
library(feather)

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


# Define params
tau.stratosphere <- 4
tau.hemisphere.inter <- 2
tau.hemisphere.intra <- 0.5
strat.nh.fraction <- 0.55

model.results <- run.model(species, tau.stratosphere = tau.stratosphere, tau.hemisphere.inter = tau.hemisphere.inter, 
                           tau.hemisphere.intra = tau.hemisphere.intra, strat.nh.fraction = strat.nh.fraction)



min <- min.cost(model.results, sf6.observations.boxed.annual.means, box.no = 4)

# Iterate over all combinations
tau.stratosphere.list = seq(1,10,0.5)
tau.hemisphere.inter.list = seq(1,5,0.5)
tau.hemisphere.intra.list = seq(0.1, 1, 0.1)
strat.nh.fraction.list = seq(0.4, 0.6, 0.1)

res.matrix = matrix(ncol=9, nrow=0)

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
      }
    }
  }
}


res.df <- data.frame(res.matrix)
colnames(res.df) <- c("t.strat", "t.hemi.inter", "t.hemi.intra", "strat.frac", "min.box.1", 
                      "min.box.2", "min.box.3", "min.box.4", "min.box.all")

head(arrange(res.df, min.box.all))

# Save the data in feather format
feather::write_feather(res.df, "results/iterative_results.feather")

# Plot the model results
#matplot(c(timestamps[1] - delta, timestamps), t(magic.matrix), type='l', lwd=4)
#legend("topleft", legend=paste("Box",1:5), col=1:5, lty=1:5, text.col=1:5, text.font=2)

# Add the observed points to the list
#yy <- ghg.observations[, "Year"]
#matpoints(yy, ghg.observations[,paste("SF6","box",1:4,sep=".")],pch=16,cex=.5)


# Generate a matrix to compate to (index 218-452)
# To compare we need to eliminate a row from magic.matrix
