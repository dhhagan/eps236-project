# Function definitions file

calc.flux.and.lifetimes = function(species, tau.stratosphere=4, tau.hemisphere.inter=2, tau.hemisphere.intra=0.5,
                                   strat.nh.fraction=0.55, print_header=F) {
  # data is molecular.wt.all
  # Source data
  #source("model_initialize.r")
  
  # Grab the mol. weight of 'species'
  mol.wt.ind.species <- molecular.wt.all[species]
  
  # Define the photochemical lifetime (years) of 'species'
  tau.global.ind.species <- tau.global.lifetimes.years[species]
  
  # Get data for 'species
  conc.observations <- ghg.observations[, c("Year", paste(species, "box", 1:4, sep='.'))]
  
  # Get species-specific loss rates from the greater loss file (by month)
  loss.freq.ind.species <- loss.monthly[, paste(species, "box", 1:5, sep=".")]
  
  # Change the rownames to be 1:12
  rownames(loss.freq.ind.species) <- 1:12
  
  # Calculate the annual frequencies of photochemical loss
  loss.freq.annual <- apply(loss.freq.ind.species, MARGIN=2, FUN=mean)
  
  # Define a vector with initial guesses for lifetimes
  initial.lifetime.guesses <- c(tau.hemisphere.intra, tau.hemisphere.inter, strat.nh.fraction, tau.stratosphere)
  
  names(initial.lifetime.guesses) <- c(
    deparse(substitute(tau.hemisphere.intra)),
    deparse(substitute(tau.hemisphere.inter)),
    deparse(substitute(strat.nh.fraction)),
    deparse(substitute(tau.stratosphere)))
  
  # definition of exch rates kg/yr
  flux.midlat <- mass.troposphere[1] / tau.hemisphere.intra   # Midlat-Tropics exchange flux kg/yr (box 1-2 advection)
  flux.ns <- sum(mass.troposphere[1:2]) / tau.hemisphere.inter # north-source hemisphere flux
  
  # Define advection flux and stuff
  flux.trop.to.strat <- mass.stratosphere / tau.stratosphere  # mass flux of air from the troposphere to the stratosphere
  
  # Define advection
  flux.adv.ind.trop.to.strat <- flux.trop.to.strat / 2.                           # box 2, 3 to stratosphere
  flux.adv.strat.to.nh <- flux.trop.to.strat * strat.nh.fraction                  # stratosphere to box 1 
  flux.adv.strat.to.sh <- flux.trop.to.strat * (1 - strat.nh.fraction)            # strat to box 4
  flux.adv.ntrop.to.strop <- flux.trop.to.strat * (2*strat.nh.fraction - 1) / 2   # box 2 to box 3 (watch for divide by 2)
  
  # compute the photochemical lifetimes to match the global lifetime and model box structure
  mass.frac.trop <- 0.842
  mass.frac.strat <- 1 - mass.frac.trop
  
  # Fudge/Estimate parameter
  alpha <- (1/tau.global.ind.species + 1/tau.stratosphere) / 
    (1/tau.stratosphere - (1/tau.global.ind.species)*(mass.frac.trop/mass.frac.strat))
  
  # Effective checical-loss lifetime in the stratosphere
  tau.strat.eff <- tau.global.ind.species * mass.frac.strat / (mass.frac.trop * alpha + mass.frac.strat)
  
  header <- sprintf(
    "
    Species: %s
    Mol Wt. %.3f
    tglobal: %.1f
    taustrat: %.1f
    strat.nh.fract: %.2f
    tau.intrahemis: %.3f
    tau.interhemis: %.3f
    alpha: %.3f
    tstrat: %.3f
    ",
    species, mol.wt.ind.species, tau.global.ind.species, tau.stratosphere, strat.nh.fraction, 
    tau.hemisphere.intra, tau.hemisphere.inter, alpha, tau.strat.eff)
  
  if (print_header == T){
    cat(header)
  }
  
  vars.to.return <- list(
    species=species, 
    mol.wt=mol.wt.ind.species,
    flux.midlat=flux.midlat, 
    flux.ns=flux.ns,
    flux.trop.to.strat=flux.trop.to.strat,
    flux.adv.ind.trop.to.strat=flux.adv.ind.trop.to.strat,
    flux.adv.ntrop.to.strop=flux.adv.ntrop.to.strop,
    flux.adv.strat.to.nh=flux.adv.strat.to.nh,
    flux.adv.strat.to.sh=flux.adv.strat.to.sh,
    tau.strat.eff=tau.strat.eff,
    mass.global=mass.global,
    mass.stratosphere=mass.stratosphere,
    mass.troposphere=mass.troposphere)
  
  return (vars.to.return)
}

run.model <- function(species, tau.stratosphere=4, tau.hemisphere.inter=2, tau.hemisphere.intra=0.5, strat.nh.fraction=0.55,
                      DEBUG = FALSE, AVG.1.YR=TRUE) {
  # Run the model and return the matrix of box-by-box results
  # First, calculate the fluxes
  # We start in 1995 because that's where the observations begin and go to 2008 because that's where EDGAR sources end
  vals <- calc.flux.and.lifetimes(species = species, tau.stratosphere = tau.stratosphere, tau.hemisphere.inter = tau.hemisphere.inter,
                                  tau.hemisphere.intra = tau.hemisphere.intra, strat.nh.fraction = strat.nh.fraction)
  
  # Set up the matrix and run  the model
  
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
  
  if (DEBUG == TRUE) {
    print("Rowsum Check")
    print(apply(k.matrix.norm, 1, sum) / apply(k.matrix.norm, 1, max)) 
    
    print("Colsum Checks")
    print(apply(k.matrix.norm, 2, sum) / apply(k.matrix.norm, 2, max)) 
    print(apply(k.matrix, 2, sum) / apply(k.matrix, 2, max)) 
  }
  
  ## compute the eigenvalues and eigenvectors
  eigen.out <- eigen(k.matrix.norm)
  
  eigen.vals <- eigen.out$values
  eigen.vecs <- eigen.out$vectors
  lifetimes <- -1 * round(1/eigen.vals, 3)
  
  # Calculate the max tropospheric lifetime
  max.tropo.lifetime <- max(lifetimes[1:4])
  
  if (DEBUG == TRUE) {
    print ("time constants (yr-1): ")
    print (lifetimes)
    print (round(eigen.vecs, 3))
  }
  
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
  
  results <- data.frame(t(magic.matrix)[-1,])
  results <- cbind(results, year=floor(as.numeric(timestamps)))
  
  if (AVG.1.YR == TRUE) {
    results <- aggregate(results, list(results$year), mean)
  }
  #year=floor(as.numeric(row.names.data.frame(sf6.observations)))
  
  return (results)
}


min.cost <- function(model.output, obs.output, box.no=1) {
  # if box.no == NA, return total sum
  # min.cost(model.results, sf6.observations.boxed.annual.means, box.no = 4)
  x1 <- obs.output[, c("SF6.box.1", "SF6.box.2", "SF6.box.3", "SF6.box.4")]
  x2 <- model.output[, c("box.1", "box.2", "box.3", "box.4")]
  
  colsums <- colSums(abs(x1 - x2))
  
  if (is.nan(box.no)) {
    res <- sum(colsums)
  }
  else {
    res <- colsums[box.no]
  }
  
  return (res)
}