#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RETREAT GRBIO 2023
# SYNTHETIC DATA
# 7/7/2023
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##-- Remove objects
rm(list=ls())

################################################################################
# Load packages
install.packages(brms)
install.packages(rstudioapi)
install.packages(readr)
install.packages(dplyr)
install.packages(ggplot2)
install.packages(reclin)

library(brms)
library(rstudioapi)
library(readr)
library(dplyr)
library(ggplot2)
library(reclin)



################################################################################
# Example 1: Sequential synthesis
# CEData
################################################################################


################################################################################
# Load data
setwd(dirname(getActiveDocumentContext()$path))
source('functions.R')
CEdata <- readr::read_csv(file = "../Datasets/CEdata.csv")
head(CEdata,3)

################################################################################
##-- Transform variables
CEdata$LogIncome      <- log(CEdata$Income)
CEdata$LogExpenditure <- log(CEdata$Expenditure)

################################################################################
# Inspect data
##-- Normality -----------------------------------------------------------------
par(mfrow=c(2,2))
hist(CEdata$Income)
hist(CEdata$Expenditure)
hist(CEdata$LogIncome)
hist(CEdata$LogExpenditure)

##-- Bivariate relationship ----------------------------------------------------
par(mfrow=c(1,1))
plot(LogExpenditure~LogIncome,CEdata)

################################################################################
# Bayesian model

##-- Get priors ----------------------------------------------------------------
brms::get_prior(data    = CEdata,
                family  = gaussian,
                formula = LogExpenditure ~ 1 + LogIncome)

##-- Prior for intercept
mean(CEdata$LogExpenditure)
sd(CEdata$LogExpenditure) * 3/(3-2)
curve(dstudent_t(x, 3, 8.9, 2.5),0,20)

##-- Fit model -----------------------------------------------------------------
SLR_fit <- brms::brm(data    = CEdata,
                     family  = gaussian,
                     formula = LogExpenditure ~ 1 + LogIncome,
                     iter    = 5000, # iterations per chain (including warmup)
                     warmup  = 3000, # warmup iterations
                     thin    = 1,    # thin rate --> Set thin > 1 to save memory and computation time
                     chains  = 1,    # Number of Markov chains
                     seed    = 539)
SLR_fit 

##-- Posteriors ----------------------------------------------------------------
post_SLR <- brms::posterior_samples(x = SLR_fit)
# post_SLR <- brms::as_draws(x = SLR_fit)
dim(post_SLR)
post_SLR[1:3, ]


##-- Diagnostics ---------------------------------------------------------------
bayesplot::mcmc_trace(x = SLR_fit,
                      pars = c("b_Intercept", "b_LogIncome", "sigma"))
bayesplot::mcmc_acf  (x = SLR_fit,
                      pars = c("b_Intercept", "b_LogIncome", "sigma"))


#-- To streamline our synthesis process, we create the design matrix X
# based on the chosen model
SLR_ff    <- stats::as.formula(LogExpenditure ~ 1 + LogIncome)
SLR_model <- stats::model.frame(SLR_ff, CEdata)
SLR_X     <- data.frame(stats::model.matrix(SLR_ff, SLR_model))



##-- Syntethic dataset ---------------------------------------------------------
n <- nrow(CEdata)
SLR_synthetic_one <- SLR_synthesize(X          = SLR_X,
                                    post_draws = post_SLR,
                                    index      = 1,
                                    n          = nrow(SLR_X),
                                    seed       = 246)
names(SLR_synthetic_one) <- c("Intercept", "LogIncome", "LogExpenditure")

##-- Comparison ----------------------------------------------------------------
d_width<- tibble(conf   = CEdata$LogExpenditure,
                 synt   = SLR_synthetic_one$LogExpenditure)
d_long <- tibble(type   = c(rep('confidential',n),rep('synthetic',n)),
                 logExp = c(CEdata$LogExpenditure,SLR_synthetic_one$LogExpenditure))

##-- Scatterplot. Relationship confidential and synthetic expenditure
ggplot(d_width, aes(x=conf, y=synt)) + 
  geom_point(size=1,alpha=0.5) +
  xlab('Confidential') + ylab('Synthetic')

##-- Density. Distribution of logExpenditure
ggplot(d_long, aes(x=logExp, color=type, linetype=type)) + 
  geom_density(linewidth=1) +
  xlab('LogExpenditure')


##-- A Bayesian simple linear regression: m > 1 synthetic datasets -------------
n <- nrow(CEdata) # Sample size
m <- 20           # Synthetic datasets 
SLR_synthetic_m <- vector("list", m)

for (l in 1:m){
  SLR_synthetic_one <- SLR_synthesize(X          = SLR_X,
                                      post_draws = post_SLR,
                                      index      = 1980 + l,
                                      n          = nrow(SLR_X),
                                      seed       = m + l)
  names(SLR_synthetic_one) <- c("Intercept", "LogIncome", "LogExpenditure")
  SLR_synthetic_m[[l]] <- SLR_synthetic_one
}
SLR_synthetic_m[[1]]

##-- Generate Kidcounts --------------------------------------------------------
Poisson_fit <- brms::brm(data   = CEdata,
                         family = poisson(link = "log"),
                         formula= KidsCount ~ 1 + LogIncome + LogExpenditure, # <-- Real log expenditure
                         iter   = 5000,
                         warmup = 3000,
                         thin   = 1,
                         chains = 1,
                         seed   = 398)

##-- Posterior parameter draws -------------------------------------------------
post_Poisson <- brms::posterior_samples(x = Poisson_fit)
post_Poisson[1:3, ]

##-- Syntethic dataset ---------------------------------------------------------
CEdata$LogExpenditure_synthetic <- SLR_synthetic_one$LogExpenditure
Poisson_ff     <- stats::as.formula(KidsCount ~ 1 + LogIncome + LogExpenditure_synthetic) # <-- Synthetic LogExpenditure
Poisson_model  <- stats::model.frame(Poisson_ff, CEdata)
Poisson_X_star <- data.frame(stats::model.matrix(Poisson_ff,Poisson_model))
n <- nrow(CEdata)

Poisson_ss_synthetic_one <- Poisson_synthesize(X          = Poisson_X_star,
                                               post_draws = post_Poisson,
                                               index      = 1,
                                               n          = nrow(Poisson_X_star),
                                               seed       = 653)
names(Poisson_ss_synthetic_one) <- c("Intercept", "LogIncome",
                                     "LogExpenditure", "KidsCount")

##-- Comparison ----------------------------------------------------------------
d_long2 <- tibble(type = factor(c(rep('confidential',n),rep('synthetic',n))),
                  kids = factor(c(CEdata$KidsCount,Poisson_ss_synthetic_one$KidsCount)))

##-- Density. Distribution of logExpenditure
ggplot(d_long2, aes(x=kids, fill=type)) + 
  geom_bar(position='dodge') +
  xlab('Type')


################################################################################
# Example 2: Joint synthesis
# CEData
################################################################################


################################################################################
# Load data
ACSdata <- data.frame(readr::read_csv(file = "../datasets/ACSdata.csv"))
head(ACSdata,3)

##-- Transform to fact
dj <- c(2, 6, 5, 2, 7, 2, 2, 3, 3, 2) # number of categories
r  <- ncol(ACSdata)

# make sure all variables are categorical
for (j in 1:r){
  ACSdata[, j] <- factor(ACSdata[, j], levels = 1:dj[j])
}

##-- Create the DPMPM model ----------------------------------------------------
## DPMPM: Dirichlet Process mixture of products of multinomials
model <- NPBayesImputeCat::CreateModel(X      = ACSdata,
                                       MCZ    = NULL,
                                       K      = 80,
                                       Nmax   = 0,
                                       aalpha = 0.25,
                                       balpha = 0.25,
                                       seed   = 221)

##-- Use the DPMPM_nozeros_syn() function to create m = 1 partially synthetic datasets
m <- 1
ACSdata_syn <- NPBayesImputeCat::DPMPM_nozeros_syn(X      = ACSdata,
                                                   dj     = dj,
                                                   nrun   = 10000,
                                                   burn   = 5000,
                                                   thin   = 10,
                                                   K      = 80,
                                                   aalpha = 0.25,
                                                   balpha = 0.25,
                                                   m      = m,
                                                   vars   = c("DIS", "HICOV"), # Only two synthetized variables
                                                   seed   = 221,
                                                   silent = TRUE)

##-- Comparison ----------------------------------------------------------------
d_long3 <- tibble(type = factor(c(rep('confidential',nrow(ACSdata)),rep('synthetic',nrow(ACSdata)))),
                  dis  = factor(c(ACSdata$DIS,ACSdata_syn$syndata[[1]]$DIS)))

##-- Density. Distribution of logExpenditure
ggplot(d_long3, aes(x=dis, fill=type)) + 
  geom_bar(position='dodge') +
  xlab('Type')


################################################################################
# R example 2
# Global Utility: pMSE
################################################################################

##-- Global utility: pMSE calculation ------------------------------------------

##-- Merge data
merged_data   <- rbind(CEdata            %>% select(Income,Expenditure),
                       SLR_synthetic_one %>% mutate(Income      = exp(LogIncome),
                                                    Expenditure = exp(LogExpenditure)) %>% 
                                             select(Income,Expenditure))
merged_data$S <- c(rep(0, n), rep(1, n))

##-- Logistic regression
log_reg <- stats::glm(formula = S ~ Income * Expenditure,
                      family  = binomial,
                      data    = merged_data)

##-- Propensity scores
probs <- predict(log_reg, type='response')


##-- pMSE
pMSE <- 1 / (2 * n) * sum((probs - 1 / 2)^2)
pMSE

#' The calculated propensity score utility measure pMSE is small and close to 0, 
#' meaning that the logistic regression model cannot really distinguish between 
#' the confidential and the synthetic datasets, indicating a high level of utility 
#' of our simulated synthetic data

################################################################################
# R example 3
# Global Utility: eCDF
################################################################################
##-- Visual comparative
d_long$Exp <- exp(d_long$logExp)
ggplot(d_long,aes(x=Exp,color=type)) + 
  stat_ecdf() + xlim(0,1000)

##-- Calculate measures  
exp_C <- d_long$Exp[d_long$type=='confidential']
exp_S <- d_long$Exp[d_long$type=='synthetic']

ecdf_C <- ecdf(exp_C)
ecdf_S <- ecdf(exp_S)

Um <- max  (ecdf_C(d_long$Exp) - ecdf_S(d_long$Exp)) 
Ua <- mean((ecdf_C(d_long$Exp) - ecdf_S(d_long$Exp))^2)
Um 
Ua

################################################################################
# R example 4
# Specific Utility: mean inference
################################################################################

##-- Specific utility: mean inference ------------------------------------------

##-- Calculate q(l) and v(l)
m <- 20
q <- numeric(m)
v <- numeric(m)
for (l in 1:m){
  SLR_synthetic_one_partial <- SLR_synthetic_m[[l]]
  q[l] <- mean(exp(SLR_synthetic_one_partial$LogExpenditure))
  v[l] <- var (exp(SLR_synthetic_one_partial$LogExpenditure))/n
}

##-- Calculate q(m) and b(m) and v(m)
q_bar_m <- mean(q)
b_m     <- var (q)
v_bar_m <- mean(v)

##-- Calculate T_p and degrees of freedom
T_p <- b_m / m + v_bar_m
v_p <- (m - 1) * (1 + v_bar_m / (b_m / m))^2

##-- Point and interval estimate
t_score_syn <- qt(p = 0.975, df = v_p)
c(q_bar_m,
  q_bar_m - t_score_syn * sqrt(T_p),
  q_bar_m + t_score_syn * sqrt(T_p))

##-- Real mean in confidential data
mean_con <- mean(CEdata$Expenditure)
sd_con <- sd(CEdata$Expenditure)
t_score_con <- qt(p = 0.975, df = n - 1)
c(mean_con,
  mean_con - t_score_con * sd_con / sqrt(n),
  mean_con + t_score_con * sd_con / sqrt(n))


################################################################################
# R example 5
# Specific Utility: Interval Overlap
################################################################################

syn_interval   <- c(q_bar_m - t_score_syn * sqrt(T_p),
                    q_bar_m + t_score_syn * sqrt(T_p))
confi_interval <- c(mean_con - t_score_con * sd_con / sqrt(n),
                    mean_con + t_score_con * sd_con / sqrt(n))

I <- CalculateIntervalOverlap_v1(confi_interval,
                                  syn_interval)
I

IO <- CalculateIntervalOverlap_v2(confi_interval,
                                   syn_interval)
IO

################################################################################
# R example 6
# Disclosure risk: Matching-based approaches
################################################################################

################################################################################
# Load data
ACSdata     <- data.frame(readr::read_csv(file = "../datasets/ACSdata.csv"))
ACSdata_syn <- data.frame(readr::read_csv(file = "../datasets/ACSdata_syn.csv"))
ACSdata_syn <- ACSdata_syn[, names(ACSdata)] ## make sure variables are in the same ordering
ACSdata_con <- ACSdata

n <- dim(ACSdata)[1]

##-- Matching-based measures
KeyQuantities_ACS <- CalculateKeyQuantities_cat(condata    = ACSdata_con,
                                                syndata    = ACSdata_syn,
                                                known.vars = c("SEX","RACE","MAR"),
                                                syn.vars   = c("DIS","HICOV"),
                                                n          = dim(ACSdata)[1])
KeyQuantities_ACS

##-- Matching-based measures
ACS_res <- IdentificationRiskCal(c_vector = KeyQuantities_ACS[["c_vector"]],
                                 T_vector = KeyQuantities_ACS[["T_vector"]],
                                 K_vector = KeyQuantities_ACS[["K_vector"]],
                                 F_vector = KeyQuantities_ACS[["F_vector"]],
                                 s        = KeyQuantities_ACS[["s"]],
                                 N        = n)

################################################################################
# R example 7
# Disclosure risk: Record Linkage approaches
################################################################################
#' In this example, YAus includes SEX, RACE, MAR and YAs includes DIS, HICOV. 
#' The intruder will perform step 1 using SEX, RACE, MAR and step 2 using DIS, 
#' HICOV.

##-- Prepare data: add index for each record
ACSdata$id <- 1:n
ACSdata_syn$id <- 1:n

#' We first generate pairs given the available variables SEX, RACE, MAR. This is 
#' done by using the pair_blocking() function in the reclin package, where the 
#' inputs include ACSdata_syn, ACSdata, and the set of available variables
ACS_pairs <- reclin::pair_blocking(ACSdata_syn, ACSdata, c("SEX", "RACE", "MAR"))

#' Next, for each generated pair in ACS_pairs, we create a set of similarity score 
#' over all the synthesized variables: the synthesized DIS, HICOV and the 
#' unsythesized LANX, WAOB, MIG, SCH, HISP
ACS_pairs_keys <- reclin::compare_pairs(ACS_pairs,
                                        by                 = c("DIS", "HICOV"),
                                        default_comparator = jaro_winkler(0.9))
ACS_pairs_keys[1:3, ]

#'We use the probabilistic record linkage approach with an EM algorithm to produce 
#'a weight for each scored pair. The higher the weight is, the more likely the pair 
#'of records belong to the same record. The EM algorithm is implemented using 
#'the problink_em() function and the weights are then calculated using the 
#'score_problink() function in the reclin package
m <- reclin::problink_em(ACS_pairs_keys)
ACS_pairs_keys_pRL <- reclin::score_problink(ACS_pairs_keys,
                                             model = m,
                                             var = "weight")

#' Now with calculated weight for each pair, we perform a one-to-one linkage by
#' comparing weights for all pairs while making sure that one record from the 
#' synthetic ACSdata_syn can be linked to at most one record from the confidential 
#' ACSdata and vice versa. There are a few choices provided by the reclin package, 
#' and for our data size, we choose to use the greedy algorithm with the 
#' select_greedy() function
ACS_pairs_keys_pRL <- reclin::select_greedy(ACS_pairs_keys_pRL,
                                            "weight",
                                            var       = "greedy",
                                            threshold = 0)
ACS_pairs_keys_pRL[1:3, ] # TRUE indicates a link and FALSE indicates no link

#'Lastly, we need to evaluate among all the TRUE’s in greedy, how many of them 
#'are true links and how many are false links
ACS_pairs_keys_pRL      <- add_from_x(ACS_pairs_keys_pRL, id_x = "id")
ACS_pairs_keys_pRL      <- add_from_y(ACS_pairs_keys_pRL, id_y = "id")
ACS_pairs_keys_pRL$true <- ACS_pairs_keys_pRL$id_x == ACS_pairs_keys_pRL$id_y

table(ACS_pairs_keys_pRL[c("true", "greedy")])
# The true linkage percentage is 6766/10000 = 67.66%, and the false
# linkage percentage is therefore 3234/10000 = 32.34%