# *****************************************************************************
#          GENERATING DIRECTLY Cmax or AUC FROM A CROSSOVER DESIGN
# *****************************************************************************

# My own package to simulate data from a crossover design (at log scale)
library(simcrossover)

# Four sequences, four periods, two treatments crossover design in
# Jones, B. and Kenward M.G. Design and Analysis of Cross-Over Trials
# CRC Press, 2015, p.127:
myDesign <- matrix(c(
  "A", "A", "B", "B",
  "B", "B", "A", "A",
  "A", "B", "B", "A",
  "B", "A", "A", "B"
), nrow = 4, byrow = TRUE)
# Just to observe fixed factors codification:
simCrossover(10, design = myDesign, treatEff = c(3,7), ssize = c(3,2,5,4), as.factors = TRUE)
# Introducing variability (an extremely heteroscedastic case:)
dat <- simCrossover(2, design = myDesign,
                    treatEff = c(3,7),
                    ssize = c(3,2,5,4), sigmaSubj = 0.294, sigma = c(0.2, 2),
                    as.factors = TRUE)
dat

# Typical bioequivalence crossover design: 2x2, TR/RT
set.seed(12345)
nsim <- 10
sim.dat <- simCrossover(nsim, ssize = c(8,8), treatEff = log(0.95), 
                        sigmaSubj = 0.1980422, sigma = 0.09975135, as.factors = TRUE)
sim.dat
attr(sim.dat, "design")
attr(sim.dat, "ssize")

# From Cmax or AUC to basic statistics (estimate of the treatment -or formulation
# in bioequivalence- effect and estimate of the "residual" or "within subjects" variability)

options(contrasts=c("contr.treatment", "contr.poly"))

# Basic computations for the first simulation replicate (column V1):
dat1.lm <- lm(sim.dat[,6+1] ~ sequence + period + treatment + subj.in.seq %in% sequence,
              data = sim.dat)
hat.treatEff <- coef(dat1.lm)["treatmentT"]
names(hat.treatEff) <-NULL
hat.treatEff
anova(dat1.lm)
hat.varWithin <- anova(dat1.lm)[5,3]
hat.varWithin

# Package which provides many bioequivalence computations
library(PowerTOST)

# Basic statistics at original scale:
gmr <- exp(hat.treatEff)
gmr
cv <- se2CV(sqrt(hat.varWithin))
cv
ci.origScal <- CI.BE(pe = gmr, CV = cv, n = attr(sim.dat, "ssize"), design = "2x2")
ci.origScal
log(ci.origScal)
if ((0.80 < ci.origScal[1]) && (ci.origScal[2] < 1.25)) {
  be <- 1
} else {
  be <- 0
}
be

# These computations must be repeated for each simulation replicate,
# i.e., for each column of 'sim.dat', from column 7 to column 16.

# Encapsultating all these computations in a function:
beComput <- function(iSim, dat) {
  dat.lm <- lm(dat[,6+iSim] ~ sequence + period + treatment + subj.in.seq %in% sequence,
               data = dat)
  hat.treatEff <- coef(dat.lm)["treatmentT"]
  names(hat.treatEff) <- NULL
  gmr <- exp(hat.treatEff)
  hat.varWithin <- anova(dat.lm)[5,3]
  cv <- se2CV(sqrt(hat.varWithin))
  ci.origScal <- CI.BE(pe = gmr, CV = cv, 
                       n = attr(dat, "ssize"), design = "2x2")
  if ((0.80 < ci.origScal[1]) && (ci.origScal[2] < 1.25)) {
    be <- 1
  } else {
    be <- 0
  }
  result <- c(formEff = hat.treatEff, varResid = hat.varWithin, GMR = gmr, CV = cv, 
              ci.origScal[1], ci.origScal[2], Bioequivalence = be)
  return(result)
}

simSummary <- vapply(1:nsim, beComput, dat = sim.dat, FUN.VALUE = rep(0.0,7))
simSummary

# Simulation estimate of the probability of declaring bioequivalence
# (very bad estimate, only 10 simulation replicates!):
sum(simSummary["Bioequivalence",]) / nsim

nsim <- 1000
set.seed(12345)
sim.dat <- simCrossover(nsim, ssize = c(8,8), treatEff = log(1.25), 
                        sigmaSubj = 0.09975135, sigma = 0.09975135, as.factors = TRUE)
simSummary <- vapply(1:nsim, beComput, dat = sim.dat, FUN.VALUE = rep(0.0,7))
simSummary[,1:10]

# Simulation estimate of the probability of declaring bioequivalence
sum(simSummary["Bioequivalence",]) / nsim

# Function that performs a "complete" simulation process (subject simulation):
sim <- function(nsim, ssize, treatEff, sigmaSubj = 0.0, sigma, seed = 12345) {
  set.seed(seed)
  sim.dat <- simCrossover(nsim, ssize = ssize, treatEff = treatEff, 
                          sigmaSubj = sigmaSubj, sigma = sigma, as.factors = TRUE)
  simSummary <- vapply(1:nsim, beComput, dat = sim.dat, FUN.VALUE = rep(0.0,7))
  result <- apply(simSummary, 1, mean)
  names(result) <- c("E(hat.formEff)", "E(hat.varWithin)", "E(hat.GMR)", "E(hat.CV)", 
                     "E(CI[1])", "E(CI[2])", "Prob{declare BE}")
  return(result)
}

sim(nsim = 10000, ssize = c(8,8), treatEff = log(1.25), sigma = 0.09975135)
sim(nsim = 10000, ssize = c(8,8), treatEff = log(0.95), sigma = 0.09975135)
sim(nsim = 10000, ssize = c(8,8), treatEff = log(0.95), sigma = 0.2935604)
sim(nsim = 10000, ssize = c(16,16), treatEff = log(0.95), sigma = 0.2935604)

# Direct simulation of the formulation effect estimate and the residual variance:
# According to the linear model and its normality assumpions:
# 1) The formulation effect estimate follows a normal distribution centered at its
#    "true" simulated population value and with 
#    variance = (1/2) * (1/n1 + 1/n2) * sigma^2
#    where sigma^2 stands for the "true" simulated residual variance and n1, n2 stand
#    for the sample sizes in each crossover sequence.
# 2) The residual variance estimate follows a 
#    (chi-square with ddf degrees of freedom) * (true residual variance) / ddf
#    where ddf = n1 + n2 - 2

nsim <- 10000
set.seed(12345)
ssize <- c(8,8)
hat.treatEffs <- rnorm(nsim, mean = log(1.25), sd = 0.09975135 * sqrt(sum(1 / (2*ssize))))
ddf <- sum(ssize) - 2
hat.varWithins <- rchisq(nsim, df = ddf) * (0.09975135^2) / ddf

beComput2 <- function(iSim, hat.treatEffs, hat.varWithins, ssize) {
  hat.treatEff <- hat.treatEffs[iSim]
  gmr <- exp(hat.treatEff)
  hat.varWithin <- hat.varWithins[iSim]
  cv <- se2CV(sqrt(hat.varWithin))
  ci.origScal <- CI.BE(pe = gmr, CV = cv, n = ssize, design = "2x2")
  if ((0.80 < ci.origScal[1]) && (ci.origScal[2] < 1.25)) {
    be <- 1
  } else {
    be <- 0
  }
  result <- c(formEff = hat.treatEff, varResid = hat.varWithin, GMR = gmr, CV = cv, 
              ci.origScal[1], ci.origScal[2], Bioequivalence = be)
  return(result)
}

simSummary2 <- vapply(1:nsim, beComput2, 
                      hat.treatEffs = hat.treatEffs, 
                      hat.varWithins = hat.varWithins, ssize = ssize, 
                      FUN.VALUE = rep(0.0,7))
simSummary2[,1:10]

# Simulation estimate of the probability of declaring bioequivalence
apply(simSummary2, 1, mean)

# Function that performs a "complete" simulation process (statistics direct generation):
sim2 <- function(nsim, ssize, treatEff, sigmaSubj = 0.0, sigma, seed = 12345) {
  set.seed(seed)
  hat.treatEffs <- rnorm(nsim, mean = treatEff, sd = sigma * sqrt(sum(1 / (2*ssize))))
  ddf <- sum(ssize) - 2
  hat.varWithins <- rchisq(nsim, df = ddf) * (sigma^2) / ddf
  simSummary <- vapply(1:nsim, beComput2, 
                       hat.treatEffs = hat.treatEffs, 
                       hat.varWithins = hat.varWithins, ssize = ssize,
                       FUN.VALUE = rep(0.0,7))
  result <- apply(simSummary, 1, mean)
  names(result) <- c("E(hat.formEff)", "E(hat.varWithin)", "E(hat.GMR)", "E(hat.CV)", 
                     "E(CI[1])", "E(CI[2])", "Prob{declare BE}")
  return(result)
}

sim2(nsim = 10000, ssize = c(8,8), treatEff = log(1.25), sigma = 0.09975135)
sim2(nsim = 10000, ssize = c(8,8), treatEff = log(0.95), sigma = 0.09975135)
sim2(nsim = 10000, ssize = c(8,8), treatEff = log(0.95), sigma = 0.2935604)
sim2(nsim = 10000, ssize = c(16,16), treatEff = log(0.95), sigma = 0.2935604)
