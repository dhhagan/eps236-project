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
