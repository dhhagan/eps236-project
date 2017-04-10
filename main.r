# main.r


# Initialize the Model
#  1) Read in data sources
source("model_initialize.r")

# Initialize Some Variables for Semi Hemi?
species <- "SF6"
tau.stratosphere <- 4  # units in years
tau.hemisphere.intra <- 0.5
tau.hemisphere.inter <- 2
strat.nh.fraction <- 0.55

# Grab the mol. weight of 'species'
mol.wt.ind.species <- molecular.wt.all[species]

# Define the photochemical lifetime (years) of 'species'
tau.global.ind.species <- tau.global[species]

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
flux.adv.ntrop.to.strop <- flux.trop.to.strat * (2*strat.nh.fraction - 1) #/ 2   # box 2 to box 3 (watch for divide by 2)

# compute the photochemical lifetimes to match the global lifetime and model box structure
alpha <- 

#   next lines compute stratospheric photochemical lifetime to match global lifetime and model box structure
alpha <- (1/tglobal + 1/taustrat)/(1/taustrat -1/tglobal*.842/.158)
tstrat <- tglobal*.158/(.842*alpha+.158) # effective chemical-loss-lifetime-in-strat
header <- (paste(species,"Mol.wt",Mol.wt,"; tglobal=",tglobal,"taustrat=",taustrat,"strat.nh.fraction=",
                 strat.nh.fraction,"tau.intrahemis=",tau.intrahemis,"tau.interhemis=",tau.interhemis,"alpha=",alpha,"tstrat=",tstrat))
print(header)