##-- Genarte synthetic data normal response
SLR_synthesize <- function(X, post_draws, index, n, seed){
  set.seed(seed)
  mean_Y <- as.matrix(X) %*%
    t(data.matrix(post_draws[index, c("b_Intercept", "b_LogIncome")]))
  synthetic_Y <- stats::rnorm(n, mean_Y, post_draws[index, "sigma"])
  data.frame(X, synthetic_Y)
}

##-- Genarte synthetic data poisson model
Poisson_synthesize <- function(X, post_draws, index, n, seed){
  set.seed(seed)
  lambda_log <- as.matrix(X) %*%
    t(data.matrix(post_draws[index, c("b_Intercept", "b_LogIncome",
                                      "b_LogExpenditure")]))
  synthetic_Y <- stats::rpois(n, exp(lambda_log))
  data.frame(X, synthetic_Y)
}

##-- Matching-based approaches
CalculateKeyQuantities_cat <- function(condata, syndata,
                                       known.vars, syn.vars, n){
  condata <- condata
  syndata <- syndata
  n <- n
  c_vector <- rep(NA, n)
  T_vector <- rep(NA, n)
  for (i in 1:n){
    match <- (eval(parse(text=paste("condata$",syn.vars,
                                    "[i]==syndata$",syn.vars,
                                    sep="",collapse="&")))&
                eval(parse(text=paste("condata$",known.vars,
                                      "[i]==syndata$",
                                      known.vars,sep="",
                                      collapse="&"))))
    match.prob <- ifelse(match, 1/sum(match), 0)
    if (max(match.prob) > 0){
      c_vector[i] <- length(match.prob[match.prob == max(match.prob)])
    }
    else
      c_vector[i] <- 0
    T_vector[i] <- is.element(i, rownames(condata)
                              [match.prob == max(match.prob)])
  }
  K_vector <- (c_vector * T_vector == 1)
  F_vector <- (c_vector * (1 - T_vector) == 1)
  s <- length(c_vector[c_vector == 1 & is.na(c_vector) == FALSE])
  res_r <- list(c_vector = c_vector,
                T_vector = T_vector,
                K_vector = K_vector,
                F_vector = F_vector,
                s = s
  )
  return(res_r)
}

#' We create the function IdentificationRiskCal() which takes the previously 
#' calculated key quantities, and N, the number of target inputs as the inputs
IdentificationRiskCal <- function(c_vector, T_vector,K_vector, F_vector,s, N){
  nonzero_c_index <- which(c_vector > 0)
  exp_match_risk <- sum(1/c_vector[nonzero_c_index] *
                          T_vector[nonzero_c_index])
  true_match_rate <- sum(na.omit(K_vector))/N
  false_match_rate <- sum(na.omit(F_vector))/s
  res_r <- list(exp_match_risk = exp_match_risk,
                true_match_rate = true_match_rate,
                false_match_rate = false_match_rate
  )
  return(res_r)
}

##-- Interval overlap definition 1 ---------------------------------------------
CalculateIntervalOverlap_v1 <- function(confi_interval, syn_interval){
  
  L_i <- max(confi_interval[1], syn_interval[1])
  U_i <- min(confi_interval[2], syn_interval[2])
  
  if (L_i <= U_i){
    overlap <- (U_i - L_i) / (2 * (confi_interval[2] - confi_interval[1])) + 
      (U_i - L_i) / (2 * (syn_interval  [2] - syn_interval  [1]))
  }else{overlap <- 0}
  
  return(overlap)
}

##-- Interval overlap definition 2 ---------------------------------------------
CalculateIntervalOverlap_v2 <- function(confi_interval, syn_interval){
  L_c <- confi_interval[1]
  U_c <- confi_interval[2]
  L_s <- syn_interval[1]
  U_s <- syn_interval[2]
  
  overlap <- 1 / 2 * ((min(U_c, U_s) - max(L_c, L_s)) / (U_c - L_c) +
                        (min(U_c, U_s) - max(L_c, L_s)) / (U_s - L_s))
  return(overlap)
}
