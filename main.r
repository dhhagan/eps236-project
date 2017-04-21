# Initialize the Model and define utility functions
source("model_initialize.r")
source("utils.r")


# Calculate params for SF6
species = "SF6"
vals <- calc.flux.and.lifetimes(species)

##################################Start definition of K matrix #########################################
k.matrix <- matrix(nrow=5, ncol=5, data=0)
  
k.matrix[1, 1] = -1*(vals$flux.adv.strat.to.nh + vals$flux.midlat)
k.matrix[1, 2] = vals$flux.midlat
k.matrix[1, 5] = vals$flux.adv.strat.to.nh

k.matrix[2, 1] = vals$flux.adv.strat.to.nh + vals$flux.midlat
k.matrix[2, 2] = -1 * (vals$flux.midlat + vals$flux.adv.ind.trop.to.strat + vals$flux.adv.ntrop.to.strop + vals$flux.ns)
k.matrix[2, 3] = vals$flux.ns

k.matrix[3, 2] = vals$flux.adv.ntrop.to.strop + vals$flux.ns
k.matrix[3, 3] = -1 * (vals$flux.midlat + vals$flux.adv.ind.trop.to.strat + vals$flux.ns)
k.matrix[3, 4] = vals$flux.adv.strat.to.sh + vals$flux.midlat

k.matrix[4, 3] = vals$flux.midlat
k.matrix[4, 4] = -1 * (vals$flux.adv.strat.to.sh + vals$flux.midlat)
k.matrix[4, 5] = vals$flux.adv.strat.to.sh

k.matrix[5, 2] = vals$flux.adv.ind.trop.to.strat
k.matrix[5, 3] = vals$flux.adv.ind.trop.to.strat
k.matrix[5, 5] = -1 * (vals$flux.adv.strat.to.nh + vals$flux.adv.strat.to.sh + (vals$mass.stratosphere / vals$tau.strat.eff))

# Normalize the k-matrix to the mass in each box of our 5-box model
# We divide by the mass in each box to convert from kg/yr to conc./yr
k.matrix.norm <- k.matrix / c(vals$mass.troposphere, vals$mass.stratosphere)

#print("Rowsum Check")
#print(apply(k.matrix.norm, 1, sum) / apply(k.matrix.norm, 1, max)) 

#print("Colsum Checks")
#print(apply(k.matrix.norm, 2, sum) / apply(k.matrix.norm, 2, max)) 
#print(apply(k.matrix, 2, sum) / apply(k.matrix, 2, max)) 


## compute the eigenvalues and eigenvectors
eigen.out <- eigen(k.matrix.norm)

eigen.vals <- eigen.out$values
eigen.vecs <- eigen.out$vectors
lifetimes <- -1 * round(1/eigen.vals, 3)

# Calculate the max tropospheric lifetime
max.tropo.lifetime <- max(lifetimes[1:4])

print ("time constants (yr-1): ")
print (lifetimes)
print (round(eigen.vecs, 3))

#################################################################
# Run the forward model with default values for the parameters
#################################################################
# global.conv is the conversion factor between Gg to mole fraction in ppt (0.029 is the mol. wt of air)
gg.to.ppt <- (1e12 * 1e9 / vals$mol.wt) / (vals$mass.global / 0.029)

# First, add a fifth row (stratosphere) with zero emissions
# Then, convert from Gg to ppt using the above-defined global.conv
sources.sf6 <- t(cbind(sf6.sources[, 2:5], rep(0, nrow(sf6.sources))))*gg.to.ppt

# Set the names of the rows and columns
dimnames(sources.sf6) <- list(1:5, sf6.sources[,1])




#note EDGAR adds 2.9 ppt bet 1995 and end 2008, but atm adds 3.23, so we should scale up total sources by 10% +-
## initial conditions, time range, and time step (.025 years, apprx 9 days)
delta <- .025;
time.range.years <- c(1995.125,2009 - delta)

# Build a vector of timestamps
timestamps <- seq(time.range.years[1],time.range.years[2], delta)  ## time period to be considered for this run

# Query the ghg.observations file for measurements between the start and end times
mod.start <- ghg.observations[,"Year"] >= time.range.years[1]
mod.end <- ghg.observations[,"Year"] <= time.range.years[2]


# Initial conditions
# Grab the index in ghg.observations which corresponds to the beginning of our measurements
obs.start.idx <- min(which(mod.start))

# Grab the Species box emissions observations for the 4 boxes from ghg emissions
# Find the time delay for 1 e-fold in the troposphere
delay <- floor(max.tropo.lifetime * 12)

obs.initial <- as.numeric(ghg.observations[obs.start.idx, paste(species, "box", 1:4, sep=".")])
obs.delayed <- as.numeric(ghg.observations[obs.start.idx + delay, paste(species, "box", 1:4, sep=".")])

# Find our initial guess for the forward solver
initial.guess <- c(obs.initial, 2*mean(obs.initial) - mean(obs.delayed))

# Create a matrix of the initial guesses
magic.matrix <- matrix(ncol=1, nrow=5, data=initial.guess)

for (tstep in timestamps) {
  gradient <- (k.matrix.norm %*% magic.matrix[, ncol(magic.matrix)] + 
                 sources.sf6[, as.character(trunc(tstep))] / c(mass.troposphere, mass.stratosphere) * mass.global) * delta
  
  # Forward step
  new.vals <- magic.matrix[, ncol(magic.matrix)] + gradient
  
  magic.matrix <- cbind(magic.matrix, new.vals)
}

dimnames(magic.matrix) <- list(paste("box", 1:5, sep='.'), c(timestamps[1] - delta, timestamps))

# Plot the model results
matplot(c(timestamps[1] - delta, timestamps), t(magic.matrix), type='l', lwd=4)
legend("topleft", legend=paste("Box",1:5), col=1:5, lty=1:5, text.col=1:5, text.font=2)

# Add the observed points to the list
yy <- ghg.observations[, "Year"]
matpoints(yy, ghg.observations[,paste("SF6","box",1:4,sep=".")],pch=16,cex=.5)