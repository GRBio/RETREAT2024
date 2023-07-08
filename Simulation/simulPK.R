# library(mlxR)
library(RsSimulx)

# #############################################################################
# DETERMINISTIC MODEL (MEAN VALUES OR JUST ONE SUBJECT WITH KNOWN PARAMETERS)
# #############################################################################
# Given a character string codifying a Mlxtran model,
# function 'inlineModel' converts it to R code (or more specifically, C++ code
# runable from R) ready to be simulated:

model1 <- inlineModel("
[LONGITUDINAL]
  input = {ka, ke, V}                  
EQUATION:
  t0 = 0
  V0 = 0
  ddt_AGI = -ka * AGI
  ddt_A1 = ka * AGI - ke * A1
  C = A1 / V
  
PK:
depot(adm = 1, target = AGI)
")

tremnt <- list(time = 0, amount = 50, adm = 1)
output <- list(name = "C", time = c(seq(0.0, 4.0, by = 0.25), 
                                    c(6.0, 8.0, 12.0, 24.0)))
params <- c(ka = 1.22, ke = 0.150, V = 58.8)

sim1 <- simulx(model = model1, treatment = tremnt, parameter = params, output = output)
head(sim1)

Cc <- sim1$C
plot1 <- plot(ggplot(data = Cc) + geom_line(aes(x=time, y=C)) 
              + xlab("hores") + ylab("[Pr. actiu]")
              + theme_bw())

# Simulation returning metrics like AUC or Cmax:
sim1.metr <- exposure(model = model1, treatment = tremnt, parameter = params, output = output)
head(sim1.metr)

# To simulate a model, 'shinymlx' directly creates a Shiny application:
library(shiny)
library(shinydashboard)
shiny1 <- shinymlx(model = model1,
                   treatment = tremnt, parameter = params, output = output,
                   style = "dashboard", appname = "shiny1")
runApp("shiny1")

# #############################################################################
#                   INTRODUCING INTER-SUBJECT VARIABILITY
# #############################################################################
# For each subject "i", the PK parameters are values ka(i) and ke(i) generated
# from a probabilistic model, possibly centered in the previously specified
# ka and ke values. More specifically, in this simulation:
# log(ka(i)) = log(ka) + Ya(i), where Ya(i) is generated according to N(0, sigma_a)
# log(ke(i)) = log(ke) + Ye(i), where Ye(i) is generated according to N(0, sigma_e).
# (In other words, parameters ka(i) i Ke(i) come from a log-normal, 
# ka(i) = ka * exp(Ya(i)),  ke(i) = ke * exp(Ye(i))

model2 <- inlineModel("
[LONGITUDINAL]
  input = {ka, ke, V, Ya, Ye}                  
EQUATION:
  t0 = 0
  V0 = 0
  ka_i = ka * exp(Ya)
  ke_i = ke * exp(Ye)
  ddt_AGI = -ka_i * AGI
  ddt_A1 = ka_i * AGI - ke_i * A1
  C = A1 / V
  
PK:
depot(adm = 1, target = AGI)

[INDIVIDUAL]

input = {sigma_a, sigma_e}

DEFINITION:
Ya = { distribution = normal, prediction = 0, sd = sigma_a}
Ye = { distribution = normal, prediction = 0, sd = sigma_e}
")

tremnt <- list(time = 0, amount = 50, adm = 1)
output <- list(name = "C", time = c(seq(0.0, 4.0, by = 0.25), 
                                    c(6.0, 8.0, 12.0, 24.0)))
params <- c(ka = 1.22, ke = 0.150, V = 58.8, 
            sigma_a = 0.1980422, sigma_e = 0.09975135)

sim2 <- simulx(model = model2, 
               treatment = tremnt, parameter = params, output = output,
               group = list(size = 12, level = "individual"),
               settings = list(seed = 12345))
head(sim2)

Cc <- sim2$C
plot2 <- plot(ggplot(data = Cc) + geom_line(aes(x=time, y=C)) 
              + xlab("hores") + ylab("[Pr. actiu]")
              + theme_bw())

prctilemlx(Cc)

# Simulation returning metrics like AUC or Cmax:
sim2.metr <- exposure(model = model2, 
                      treatment = tremnt, parameter = params, 
                      output = output,
                      group = list(size = 12, level = "individual"),
                      settings = list(seed = 12345))
head(sim2.metr)


# #############################################################################
#     INTRODUCING INTRA-SUBJECT (or RESIDUAL or WITHIN-SUBJECT) VARIABILITY
# #############################################################################
# (No he trobat una forma de fer-ho dins aquest programari)
sim3 <- sim2
sim3$C[,"C"] <- sim3$C[,"C"] * exp(rnorm(nrow(sim3$C), sd = 0.09975135))

head(sim2$C[,"C"])
head(sim3$C[,"C"])

Cc <- sim3$C
plot2 <- plot(ggplot(data = Cc) + geom_line(aes(x=time, y=C)) 
              + xlab("hores") + ylab("[Pr. actiu]")
              + theme_bw())

prctilemlx(Cc)
