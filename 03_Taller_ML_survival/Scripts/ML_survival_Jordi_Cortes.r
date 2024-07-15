################################################################################
# Course: Machine learning in survival analysis
# RETREAT GRBIO 11th July 2024
# Author: Jordi Cortes Martinez
# Contact: jordi.cortes-martinez@upc.edu
################################################################################

rm(list=ls())

################################################################################
## Load libraries
################################################################################
library(survival)        # Survival analysis
library(caret)           # ML
library(dplyr)           # Tidy data
library(ggplot2)         # Tidy graphics
library(survminer)       # Plot survival curves
library(pec)             # Measures of predictive performance
library(glmnet)          # Measures of predictive performance
library(rpart)           # Survival tree
library(partykit)        # Survival tree
library(randomForestSRC) # Random Forest
library(survivalsvm)     # Survival SVM
library(survivalmodels)  # deepsurv
library(reticulate)      # import python packages
library(MTLR)            # MTLR
# remotes::install_github("IyarLin/survXgboost")
library(survXgboost)     # survXgboost
library(xgboost)         # xgboost
library(data.table)      # objects of xgboost
library(Matrix)          # objects of xgboost
 
################################################################################
## Descriptive analysis and data preparation
################################################################################

##-- Inspection ----------------------------------------------------------------
View(rotterdam)
?rotterdam
d0 <- as_tibble(rotterdam)

##-- Clean data ----------------------------------------------------------------
d <- d0 %>% 
  select(-c(pid,rtime,recur)) %>% 
  rename(time = dtime, status  = death)

##-- Some preliminary analysis -------------------------------------------------
surv_obj <- with(d,Surv(time,status))

survdiff(surv_obj ~ meno,  data = d)            # Log-rank test
survdiff(surv_obj ~ size,  data = d)            # Log-rank test

ggsurvplot(survfit(surv_obj ~ meno, data = d))  # Survival curves
ggsurvplot(survfit(surv_obj ~ size, data = d))  # Survival curves

##-- Train & Test data ---------------------------------------------------------
set.seed(123)
trainIndex <- createDataPartition(d$status, p = .7, list = FALSE, times = 1)[,1]
d_train    <- d %>% dplyr::slice( trainIndex) # Between 60% and 80%
d_test     <- d %>% dplyr::slice(-trainIndex)
dim(d_train);dim(d_test)

#-------------------------------------------------------------------------------
#
# Traditional methods
#
#-------------------------------------------------------------------------------

################################################################################
## Kaplan-Meier (naive)
################################################################################

##-- Model ---------------------------------------------------------------------
fit_km   <- survfit(Surv(time, status) ~ 1, data = d_train)            # Kaplan-Meier survival curve
plot(fit_km)
med_time <- summary(fit_km)$table['median']                            # median time
(pred_km <- 1 - predictSurvProb(fit_km, d_test, times = med_time)[,1]) # prob of death in the median time (=0.5)
y        <- cbind(time = d_test$time, status = d_test$status)          # time and status for test sample
(Cind_km_test <- Cindex(pred = pred_km, y = y))                        # Naive c-index

################################################################################
## Cox model
################################################################################
##-- Model  --------------------------------------------------------------------
fit_cox <- coxph(Surv(time,status) ~ ., data=d_train, x = TRUE) 
summary(fit_cox)

##-- C-Index in train data -----------------------------------------------------
(Cind_cox_train <- fit_cox$concordance['concordance'])

##-- C-Index in test data ------------------------------------------------------
pred_cox <- predict(fit_cox, d_test)                          # linear predictor (it does not depend on time). It uses centered continuous covariates
y        <- cbind(time = d_test$time, status = d_test$status) # time and status for test sample
(Cind_cox_test <- Cindex(pred = pred_cox, y = y))             # C-index for cox model

#-------------------------------------------------------------------------------
#
# ML methods
#
#-------------------------------------------------------------------------------

################################################################################
## Survival tree
################################################################################

##-- Model ---------------------------------------------------------------------
fit_stree <- as.party(rpart(Surv(time, status) ~ . , data = d_train,  
                            control = rpart.control(minsplit=30, cp=0.008)))

##-- Visualization -------------------------------------------------------------
plot(fit_stree, type = 'extended')

##-- C-Index in test data ------------------------------------------------------
pred_stree <- - predict(fit_stree, newdata = d_test, type='response')
y  <-  cbind(time = d_test$time, status = d_test$status)
(Cind_stree_test <- Cindex(pred = pred_stree, y = y)) # 0.66

################################################################################
## Survival Random forests
################################################################################
##-- Model ---------------------------------------------------------------------
d_train <- as.data.frame(d_train)       # It does not work as tibble for factors
fit_rf  <- rfsrc(Surv(time, status) ~ ., data = d_train, 
                 ntree = 500, splitrule="logrank")

##-- Values returned by rfsrc object --------------------------------
# chf	          --> Cumulative hazard function (CHF).
# chf.oob	      --> OOB CHF (NULL unless outcome="test").
# survival	    --> Survival function.
# survival.oob	--> OOB survival function (NULL unless outcome="test").
# time.interest	--> Ordered unique death times.
# ndead	        --> Number of deaths.

##-- CHF for each individual ---------------------------------------------------
time    <- fit_rf$time.interest
chf     <- fit_rf$chf
chf.oob <- fit_rf$chf.oob

par(mfrow=c(1,2),las=1)

##-- Ensemble (In bag)
plot(time,chf[1,],type='s',xlab='Time',ylab='CHF', 
     main='Ensemble',sub = "(10 individuals)")
apply(chf[2:10,],1,lines,x=time,type='s',col='grey80')

##-- OOB (Out of bag)
plot(time,chf.oob[1,],type='s',xlab='Time',ylab='CHF', 
     main='OOB',sub = "(10 individuals)")
apply(chf.oob[2:10,],1,lines,x=time,type='s',col='grey80')

##-- Survival functions for each individual ------------------------------------
time     <- fit_rf$time.interest
surv     <- fit_rf$survival
surv.oob <- fit_rf$survival.oob

par(mfrow=c(1,2),las=1)

##-- Ensemble (In bag)
plot(time,surv[1,],type='s',xlab='Time',ylab='S', 
     main='Ensemble',ylim=0:1,sub = "(10 individuals)")
apply(surv[2:10,],1,lines,x=time,type='s',col='grey80')

##-- OOB (Out of bag)
plot(time,surv.oob[1,],type='s',xlab='Time',ylab='S', 
     main='OOB',ylim=0:1,sub = "(10 individuals)")
apply(surv.oob[2:10,],1,lines,x=time,type='s',col='grey80')

##-- Variable importance  ------------------------------------------------------
fit_rf_imp1  <- vimp(fit_rf)$importance  # Do not run --> It takes a lot of time
fit_rf_imp2  <- subsample(fit_rf)        # Do not run --> It takes a lot of time
load(file = 'Rdata/delay_objects.Rdata')   # saved object with the executions
fit_rf_imp1       # w/o bootstrapping
par(mfrow=c(1,1))
plot(fit_rf_imp2) # by bootstrapping
# save(list=c('fit_rf_imp1','fit_rf_imp2'),file = 'Rdata/delay_objects.Rdata')

##-- C-Index in test data (in median time) -------------------------------------
med_time_pos <- which(fit_rf$time.interest>med_time)[1]  # closest higher value to median
pred_rf  <- predict(fit_rf,d_test)$chf[,med_time_pos]    # chf in median time
y  <-  cbind(time = d_test$time, status = d_test$status) # time and status in test sample
(Cind_rf_test <- Cindex(pred = pred_rf, y = y))          # 0.69

##-- More info -----------------------------------------------------------------
##-- Tune mtry and nodesize
opt_param   <- tune(Surv(time, status) ~ ., 
                    data     = d_train, 
                    ntreeTry = 50)
fit_rf_opt  <- rfsrc(Surv(time, status) ~ ., 
                     data      = d_train, 
                     mtry      = opt_param$optimal['mtry'],
                     nodesize  = opt_param$optimal['nodesize'],
                     ntree     = 500, 
                     splitrule = "logrank")

##-- C-Index in test data (with tune) 
pred_rf_opt       <- predict(fit_rf_opt,d_test)$chf[,med_time_pos] # chf in median time
(Cind_rf_test_opt <- Cindex(pred = pred_rf_opt, y = y))            # It does not improve

##-- Brier score
brier_rf_train <- get.brier.survival(fit_rf, cens.model="km")
plot(brier_rf_train$brier.score, type="s", col=2)

################################################################################
## Survival SVM
################################################################################

##-- Kernel SVM ranking model --------------------------------------------------
fit_svm_ra <- survivalsvm(Surv(time, status) ~ .,     # Formula
                          data     = d_train[1:500,], # Only 500 observations 
                          type     = "vanbelle1",     # Type of SVM
                          diff.meth= "makediff1",
                          gamma.mu = 1,               # Penalization parameter 
                          opt.meth = "quadprog",      # Optimization method 
                          kernel   = "lin_kernel")    # Kernel to map

##-- C-Index in test data ------------------------------------------------------
pred_svm_ra <- predict(fit_svm_ra, d_test)$predicted[1,]
y           <- cbind(time = d_test$time, status = d_test$status)
(Cind_svm_ra_test <- Cindex(pred = -pred_svm_ra, y = y)) # 0.58


################################################################################
## DeepSurv
################################################################################
# reticulate::install_python(version = '3.12.1')
# reticulate::py_install("pycox")
# reticulate::py_install("torch")
sys <- import("sys", convert = TRUE)

##-- Model Deep Surv -----------------------------------------------------------
set.seed(12345)
system.time(
fit_deep <- deepsurv(data           = d_train,           # data
                     frac           = 0.3,               # Fraction of data to validate dataset          
                     activation     = "relu",            # Activation function: REctified Linear Units. Piecewise linear function that will output the input directly if positive, otherwise, it will output zero 
                     num_nodes      = c(4, 8, 4, 2),     # number of nodes in each layer
                     dropout        = 0.1,               # dropout fraction tuned over [0, 1]
                     early_stopping = TRUE,              # Early stopping
                     batch_size     = 32,                # Elements in each batch
                     epochs         = 100)               # number of epochs
) # 31 seconds

##-- Interpretability ----------------------------------------------------------
# Influence of age and size on the risk
pred_deep_train    <- predict(fit_deep, type = "risk")
d_interpret <- tibble(age  = d_train$age,
                      size = d_train$size,
                      risk = pred_deep_train)
ggplot(d_interpret,aes(x = age, y = risk)) + 
  geom_smooth(se = TRUE)
ggplot(d_interpret,aes(x = size, y = risk)) + 
  geom_boxplot()

##-- Size influence on individual 1 (Personalized) -----------------------------
d_size      <- tibble(rbind(d_test[1,],d_test[1,],d_test[1,])) # Repeat individual 3 times
d_size$size <- levels(d_test$size)                             # Assign different sizes

##-- Risk for individual 1 according size
fit_deep_risk1 <- predict(fit_deep, type = "risk", d_size)     # Different risk values
fit_deep_risk1/fit_deep_risk1[1]                               # HRs (first is the reference)

##-- Survival for individual 1 according size 
d_size_int  <- predict(fit_deep, type = "surv", d_size)
times <- as.numeric(colnames(d_size_int))
matplot(times, t(d_size_int), type = "l", ylim=0:1,
        xlab='Time', ylab='S(t)', las=1)
legend('bottomleft',legend = paste0('Ind 1. Size: ',levels(d_test$size)),
       col = 1:3,lty=1:3, cex=0.8)

##-- C-Index in test data ------------------------------------------------------
pred_deep <- predict(fit_deep, type = "risk", d_test)
y         <- cbind(time = d_test$time, status = d_test$status)
(Cind_deep_test <- Cindex(pred = pred_deep, y = y)) # 0.69


################################################################################
## MTLR Model - machine learning model
################################################################################
##-- Model MTLR ----------------------------------------------------------------
fit_mtlr   <- mtlr(Surv(time, status)~. , data = d_train, nintervals = 9)

##-- Variable importance over time ---------------------------------------------
plot(fit_mtlr) # Absolute value indicates importance and the sign indicates the direction of association

##-- Survival curves -----------------------------------------------------------
pred_mtlr_surv <- predict(fit_mtlr)
plotcurves(pred_mtlr_surv, index = 1:10, remove_legend = FALSE) + 
  theme(legend.position = "none") # First 10 individuals
  
##-- C-Index in test data ------------------------------------------------------
pred_mtlr <- predict(fit_mtlr, d_test, type="mean_time")
y         <- cbind(time = d_test$time[sel], status = d_test$status[sel])
(Cind_mtlr_test <- Cindex(pred = -pred_mtlr, y = y)) # 0.69





################################################################################
## survXgboost
################################################################################

##-- Prepare data -------------------------------------------------------------------
label_train <- ifelse(d_train$status == 1, d_train$time, -d_train$time)
label_test  <- ifelse(d_test$status == 1,  d_test$time,  -d_test$time)

x_train <- as.matrix(d_train[, !names(d_train) %in% c("size","time", "status")]) # it does not work with factors
x_test0 <- as.matrix(d_test[, !names(d_train) %in% c("size","time", "status")])

x_test <- xgb.DMatrix(data  = as.matrix(x_test0), label = label_test)

##-- Model ---------------------------------------------------------------------
fit_xgboost <- xgb.train.surv(
  params = list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.05), # larger eta leads to algorithm not converging
  data      = x_train, 
  label     = label_train,
  watchlist = list(val2 = x_test),
  nrounds   = 1000, 
  early_stopping_rounds = 30
)

##-- Predict survival curves ---------------------------------------------------
times           <- seq(1, max(d_train$time), 50)
survival_curves <- predict(fit_xgboost, newdata = x_train, type = "surv", times = times)
matplot(times, t(survival_curves[1:5, ]), type = "l", ylim=0:1,
        xlab='Time', ylab='S(t)', las=1)

##-- Variable importance (only status) -----------------------------------------
sparse_matrix <- sparse.model.matrix(status~. - time, data = d_train)
dfxg          <- data.table(d_train, keep.rownames = FALSE)
output_vector <- dfxg$status
bst <- xgboost(data      = sparse_matrix, 
               label     = output_vector, 
               max.depth = 4, 
               eta       = 1, 
               nthread   = 2, 
               nrounds   = 10,
               objective = "binary:logistic")
importance <- xgb.importance(feature_names = sparse_matrix@Dimnames[[2]], model = bst)
xgb.plot.importance(importance_matrix = importance)

##-- C-Index in test data ------------------------------------------------------
pred_xgboost       <- predict(fit_xgboost, newdata = x_test, 
                              type = "surv", times = med_time)
(Cind_xgboost_test <- Cindex(pred = 1-pred_xgboost, y = y)) 


################################################################################
# Conclusion
################################################################################
Models <- c('KM','Cox','Stree','RandForest','SVM_ra','DeepSurv','MTLR','xgboost')
d_summary <- data.frame(Model  = factor(Models,levels=Models),
                        Cindex = c(Cind_km_test,
                                  Cind_cox_test,
                                  Cind_stree_test,
                                  Cind_rf_test,
                                  Cind_svm_ra_test,
                                  Cind_deep_test,
                                  Cind_mtlr_test,
                                  Cind_xgboost_test))

ggplot(d_summary, aes(x = Model, y = Cindex)) + geom_col()



