# Model Initialization
# Read in experimental data for emissions and losses
# Initalize global variables for mass of atmosphere, troposphere, and the individual compounds

##   GLOBAL LIFETIMES in YEARS  (REQUIRED file):  vvvvvvvvvvvvvv
tau.global.lifetimes.years <- scan("data/tglobal_all.txt",skip=1)

# Name the columns in the tglobal.all array
names(tau.global.lifetimes.years) <- scan("data/tglobal_all.txt", nlines=1, what=character())

#  Read in NASA Greenhouse Gas Global Monitoring Devision values
ghg.observations <- read.table("data/ghg.GMD.conc.tbl", header = TRUE)

# Read in losses as a function of month
loss.monthly <- read.table("data/loss.monthly.tbl", header=TRUE)

# EDGAR 2015 emissions of SF6 by box, 1970 to 2008; for 2008 to 2014, we will use our derived model!
# Read in SF6 Emissions Table from EDGAR
sf6.sources <- read.table("data/SF6.emissions.bybox.tbl", header=TRUE)

# Select the time range that is covered by the EDGAR SF6 Data (min, max)
time.range.sf6.sources <- range(sf6.sources[,"yr"])

# Calculate the masses in each of the boxes
mass.global <- ( 5.1e14 * 0.984e5)/9.8   #5.27e18 kg; M*g g=9.8  M=0.984/9.8 kg/m2  surface area of earth=5.1e14
mass.stratosphere <- 0.1577 * mass.global

#semi-hemi masses         N-temp      N-trop       S-trop    S-temp (subtr strat mass from each)
mass.troposphere = (1 - c(exp(-12/7), exp(-14/7), exp(-14/7), exp(-12/7) ) ) * mass.global/4

# Make molecular weights in terms of just Carbon and Nitrogen
molecular.wt.all <- c(28, 137.37, 32.065 + 6*18.998, 120.91, 86.47, 12, 102, 117)

# Name the columns for the molecular weights vector
names(molecular.wt.all) <- c("N2O", "CFC11", "SF6", "CFC12", "HCFC22", "CO2", "HFC134A", "HCFC141B")