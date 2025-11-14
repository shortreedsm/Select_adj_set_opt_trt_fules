library(dplyr)
library(splines)
library(tidyr)
library(glmnet)
library(R.utils)
library(Ball)
library(grpreg)
library(glmpath)
library(np)
library(CBPSmod)
library(parallel)

sink("analysis_out.txt")

source("GLiDeR_Rfunctions_tailoring.R")
source("causal_ball_script_tailoring.R")
source("dips_tailoring.R")

#create expanded dataset
# rev - direction of rule, F for I(P0<theta), T for I(P0>=theta)
create.pseudo <- function(dat, t_ind, Theta_vec, rev) {
  n <- nrow(dat)
  dat.pseudo <- cbind(dat, 0)
  names(dat.pseudo) <- c(names(dat), "regime")
  dat.pseudo <- dat.pseudo[NULL, ]
  A <- dat$A
  
  P0 <- dat[, t_ind]
  
  for (ii in 1:(length(Theta_vec))) {
    if(rev){
      temp.index <-
        (A == 1 & P0 >= Theta_vec[ii]) | (A == 0 & P0 < Theta_vec[ii])
    } else{
      #regime of interest: A = 1(P0<theta)
      temp.index <-
        (A == 1 & P0 < Theta_vec[ii]) | (A == 0 & P0 >= Theta_vec[ii])
    }
    if (sum(temp.index) == 0)
      next
    dat.p <-
      cbind(dat[temp.index, ], rep(Theta_vec[ii], sum(temp.index)))
    names(dat.p) <- c(names(dat), "regime")
    dat.pseudo <- rbind(dat.pseudo, dat.p)
  }
  dat.pseudo <- dat.pseudo[order(dat.pseudo$regime), ]
  return(dat.pseudo)
}

 #propensity methods -----
#propensity model with no selection
no_selection <- function(dat, adjust_inds, t_ind) {
  
  var_num <- ncol(dat) - 2
  
  inds <- union(adjust_inds, t_ind)
  
  fit.treatment <- glm(dat$A ~ as.matrix(dat[, (inds)]),
                       family = binomial())
  
  
  dat$pred.treatment <- fit.treatment$fitted
  
  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment > 1),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > 1 - 10 ^ (-6))
  )
  
  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))
  
  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)
  
  return(list(dat = dat, selected = rep(NA, var_num),
              bad_props = bad_props))
  
}


#weighted absolute mean difference function to choose lamdba weight in OAL
#arguments
#    betas - vector of beta values from outcome model
#    X - covariate vector
#    A - treatment vector
#    weights - matrix of fitted weights, a column per lamdba
#returns
# wAMD
wAMD <- function(betas, X, A, weights) {
  d <- length(betas)
  
  exposed_denom <- A %*% weights
  
  unexposed_denom <- (1 - A) %*% weights
  
  exposed_num <- matrix(NA, nrow = d,
                        ncol = ncol (weights))
  
  unexposed_num <- exposed_num
  
  for (j in 1:d) {
    exposed_num[j, ] <- (X[, j] * A) %*% weights
    
    unexposed_num[j, ] <- (X[, j] * (1 - A)) %*% weights
    
  }
  
  wAMDs <- abs(betas) %*%
    abs(
      sweep(exposed_num, 2, exposed_denom, "/")
      - sweep(unexposed_num, 2, unexposed_denom, "/")
    )
  
  return(wAMDs)
  
}


oal <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 2
  
  inds <- c(setdiff(adjust_inds, t_ind), t_ind)
  
  X <- as.matrix(dat[, (inds)])
  
  d <- length(inds)
  
  n <- nrow(dat)
  
  #lambda <- n^c(-10,-5,-1,-0.75,-0.5,-0.25,0.25,0.49)
  
  orig_lambda <- n ^ c(-20, -15, -10, -5, -1, -0.75, 0, 0.49)
  
  temp <- as.matrix(cbind(dat$A, X))
  
  outcome_betas <- lm(dat$Y ~ temp)$coefficients[-(1:2)]
  
  gamma <- 6 - 2 * log(orig_lambda, base = n)
  
  b_weights <- matrix(NA, nrow = (d),
                      ncol = length(gamma))
  
  for (j in 1:(d - 1)) {
    b_weights[j, ] <- abs(outcome_betas[j]) ^ (-gamma)
  }
  
  #tailoring ind has 0 weight
  b_weights[d, ] <- 0
  
  scale_factors <- apply(b_weights, 2, sum) / d
  
  lambda <- orig_lambda * scale_factors

  weights <- matrix(NA, nrow = nrow(dat), ncol = length(lambda))
  
  
  for (r in 1:length(lambda)) {
    
    fit.treatment <- glmnet(
      x = X,
      y = dat$A,
      family = binomial(link = "logit"),
      alpha = 1,
      lambda = lambda[r],
      penalty.factor = b_weights[, r]
    )
    
    fits <- predict(fit.treatment,
                    X,
                    s = lambda[r],
                    type = "response")
    
    weights[, r] <- 1 / fits * dat$A +
      (1 - dat$A) * 1 / (1 - fits)
  }
  
  
  wAMD_vals <- wAMD(
    betas = outcome_betas,
    X = X,
    A = dat$A,
    weights = weights
  )
  
  
  best_lambda_ind <- which.min(wAMD_vals)
  

  final_mod <- glmnet(
    x = X,
    y = dat$A,
    family = binomial(link = "logit"),
    alpha = 1,
    lambda = lambda[best_lambda_ind],
    penalty.factor = b_weights[, best_lambda_ind]
  )
  
  fits <-   predict(final_mod,
                    X,
                    s = lambda[best_lambda_ind],
                    type = "response")
  
  coefs <-   coef(final_mod,
                  s = lambda[best_lambda_ind])[-1]
  
  
  #get selected variables in order of data 
  temp_selected <- as.vector(abs(coefs) > 10^(-5))

  
  
  selected <- rep(0, var_num)
  
  
  
  #get the true selected indices
  temp_selected <- inds[temp_selected]
  
  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1
  
  dat$pred.treatment <- fits
  
  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment > 1),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > 1 - 10 ^ (-6))
  )
  
  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))
  
  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)
  
  return(list(dat = dat, selected = selected, bad_props = bad_props))
  
}

glider <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 2
  
  #put tailoring index last, glider will set last index weight to 0
  inds <- c(setdiff(adjust_inds, t_ind), t_ind)
  
  mod.treatment <- GLiDeR(
    Xorig = dat[, (inds)],
    Yorig = dat$Y ,
    Aorig = dat$A,
    lambda = NULL
  )
 
  
  X <- as.matrix(cbind(1, dat[, (inds)]))
  
  gamma <- c(mod.treatment$gamma0, mod.treatment$gamma)
  
  dat$pred.treatment <- 1 / (1 + exp(X %*% gamma))
  
  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment > 1),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > 1 - 10 ^ (-6))
  )
  
  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))
  
  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)
  
  #get selected variables in order of input data 
  temp_selected <- abs(mod.treatment$gamma) > 10^(-5)

  
  selected <- rep(0, var_num)
  
  #get the true selected indices
  temp_selected <- inds[temp_selected]
  
  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1
  
  return(list(dat = dat, selected = selected, bad_props = bad_props))
  
}

causal_ball <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 2
  
  #set tailoring index as last
  inds <- c(setdiff(adjust_inds, t_ind), t_ind)
  
  #tell CBS that tailoring index is last
  
#use a subset of the rows so that memory wont be too much
 r_subset <- sample(nrow(dat), size = round(.4*nrow(dat)))
  
  temp <- CBS(
    X = dat[r_subset, (inds)],
    D = dat$A[r_subset],
    Y = dat$Y[r_subset],
    t_ind = length(inds)
  )
  
  temp_selected <- temp$selected
  
  selected <- rep(0, var_num)
  
  #get the true selected indices
  temp_selected <- inds[temp_selected == 1]

  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1
  
  dat$pred.treatment <- predict(temp$final_mod,
                                as.matrix(dat[,inds[temp$lasso_inds]]), 
                                type = "response")
  
  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment > 1),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > 1 - 10 ^ (-6))
  )
  
  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))
  
  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)
  
  return(list(dat = dat, selected = selected, bad_props = bad_props))
  
}

dip <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 2
  
  inds <- c(setdiff(adjust_inds, t_ind), t_ind)
  
  temp <-
    DiPS(
      Yi = dat$Y,
      Ti = dat$A,
      Xi = as.matrix(dat[, (inds)]),
      Wi = rep(1, nrow(dat)),
      bwd = "plug",
      q = 4,
      alp = 1 / (2 + 4),
      C = 1,
      constEff = FALSE,
      t_ind = length(inds)
    )
  
  bad_props <- c(mean(c(temp[[1]],temp[[2]]) < 0),
                 mean(c(temp[[1]],temp[[2]]) == 0),
                 mean(c(temp[[1]],temp[[2]]) > 0 & c(temp[[1]],temp[[2]]) < 10 ^ (-6)),
                 mean(c(temp[[1]],temp[[2]]) > 1),
                 mean(c(temp[[1]],temp[[2]]) == 1),
                 mean(c(temp[[1]],temp[[2]]) < 1 & c(temp[[1]],temp[[2]]) > 1 - 10 ^ (-6))
  )
  
  #put weights in [0,1] range
  probs_1 <- pmin(pmax(temp[[1]], 10 ^ (-6)), 1 - 10 ^ (-6))
  probs_0 <- pmin(pmax(temp[[2]], 10 ^ (-6)), 1 - 10 ^ (-6))
  
  dat$weights.trt <- dat$A * 1 / probs_1 +
    (1 - dat$A) * 1 / (1-probs_0)
  
  dat$pred.treatment <- probs_1
  
  #get selected variables in order of data in the DiPS function
  temp_selected <- temp[[3]]
  
  selected <- rep(0, var_num)
  
  #get the true selected indices
  temp_selected <- inds[temp_selected]
  
  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1
  
  return(list(dat = dat, selected = selected, bad_props = bad_props))
  
}

hd_balancing <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 2
  
  X <- as.matrix(dat[, (adjust_inds)])
  
  A <- dat$A
  
  dat$pred.treatment <- NA
  
  
 # withTimeout({
    dat$pred.treatment <- hdCBPS(
      formula = A ~ X,
      y = dat$Y,
      ATT = 0,
      #default iterations is 1000
      iterations = 100
    )$fitted.values
  #}, onTimeout = "error",
  #timeout = 3600)
    
    bad_props <- c(mean(dat$pred.treatment < 0),
                   mean(dat$pred.treatment == 0),
                   mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                   mean(dat$pred.treatment > 1),
                   mean(dat$pred.treatment == 1),
                   mean(dat$pred.treatment < 1 & dat$pred.treatment > 1 - 10 ^ (-6))
    )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))
    
  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)
  
  return(list(dat = dat, selected = rep(NA, var_num), bad_props = bad_props))
  
}


#treatment weights from logistic regression with baseline variables
#arguments:
#  dat - data frame with tx var A at least
#  adjust_inds - vector of columns of any confounders in dat
#  prop.selection - selection method
create.analysis <- function(dat, adjust_inds, t_ind, prop.selection) {
  var_num <- ncol(dat) - 2
  
  select_fun <- eval(as.name(prop.selection))
  
  #baseline treatment weights
  if (is.null(adjust_inds)) {
    fit.treatment <-
      glm(dat$A ~ as.matrix(dat[, (t_ind)]), family = binomial())
    
    dat$pred.treatment <- fit.treatment$fitted
    
    dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
      (1 - dat$A) * 1 / (1 - dat$pred.treatment)
    
    dat$weights.trt  <- pmin(dat$weights.trt, 10)
    
    return(list(dat = dat, selected  = rep(0, var_num)))
  }
  else {
    
    ret <- select_fun(dat, adjust_inds, t_ind)
    
    ret$dat$weights.trt <- pmin(ret$dat$weights.trt, 10)
    
    return(ret)
  }
  
}

# estimation function -----

#arguments:
# dat: dataframe of observations with IPW weights
# t_ind - the index of the tailoring variable under consideration
# prop_adjust - logical, whether to adjust for covariates in the propensity
# prop_adjust_inds - inds of vars to adjust for in propensity
# outcome_adjust - how to adjust for covariates in the marginal:
#    via "Lasso", via the "Oracle", or "None"
# outcome_adjust_inds - inds of vars to adjust for besides regime in marginal
# Theta_vec - vector of regime values under consideration
# knots - vector of spline knots
# prop.selection - propensity selection method
# rev - direction of rule, F for I(P0<theta), T for I(P0>=theta)
#returns:
# v - estimated value of regime correspondng to each entry of Theta_vec
# beta - coefficients of the model
estimate.regime <- function(dat,
                            t_ind,
                            prop_adjust,
                            prop_adjust_inds,
                            outcome_adjust,
                            outcome_adjust_inds,
                            Theta_vec,
                            knots,
                            prop.selection,
                            rev)
{
  #change from v6, IPW weights are specific to variable
  temp <- create.analysis(dat, prop_adjust_inds, t_ind, prop.selection)
  
  bad_props <- temp$bad_props

  
  #make sure tailoring variable isn't in margianl model
  outcome_adjust_inds <- setdiff(outcome_adjust_inds,t_ind)
  
  selected <- temp$selected
  
  dat <- temp$dat
  
  dat.pseudo <- create.pseudo(dat, (t_ind), Theta_vec, rev = rev)
  
  X_marg <- ns(x = dat.pseudo$regime, knots = knots)
  
  
  #add in adjustment variables
  if (outcome_adjust == "Lasso") {
    X <- cbind(X_marg, dat.pseudo[, (outcome_adjust_inds)])
    
    X <- as.matrix(X)
    
    penalty.factor <- c(rep(0, ncol(X_marg)),
                        rep(1, ncol(dat.pseudo[, (outcome_adjust_inds)])))
    
    mod <- glmnet(
      y = dat.pseudo$Y,
      x = X,
      alpha = 1,
      penalty.factor = penalty.factor,
      weights = dat.pseudo$weights.trt
    )
    
    
    pred <- predict(mod, newx = X)
    
    bic <- rep(NA, ncol(pred))
    
    for (j in 1:ncol(pred)) {
      #gets the weighted mse and number of knots, calculates criteria val
      wmse <- mean(dat.pseudo$weights.trt * (dat.pseudo$Y - pred[, j]) ^
                     2)
      
      p <- mod$df[j] - 1
      
      
      bic[j] <- log(wmse) * nrow(X) + log(nrow(X)) * p
    }
    
    best_s <- mod$lambda[which.min(bic)]
    
    v_selected <- as.vector(coef(mod, s = best_s))[-1] != 0
    
    res <-
      lm(dat.pseudo$Y ~ X[, v_selected], weights = dat.pseudo$weights.trt)
    
    beta <- as.vector(coef(mod, s = best_s))
    
    beta[c(TRUE, v_selected)] <- res$coefficients
    
    v <- rep(NA, length(Theta_vec))
    
    X_pred_all <- ns(x = Theta_vec, knots = knots)
    
    for (i in 1:length(Theta_vec)) {
      X_pred <- X_pred_all[i, ]
      
      X_pred <-
        data.frame(1, t(X_pred), dat[, (outcome_adjust_inds)])
      
      v[i] <- mean(as.matrix(X_pred) %*% beta)
    }
    
  } else {
    res <- lm(dat.pseudo$Y ~ X_marg, weights = dat.pseudo$weights.trt)
    
    beta <- res$coefficients
    
    v <- rep(NA, length(Theta_vec))
    
    X_pred <- ns(x = Theta_vec, knots = knots)
    
    X_pred <- as.matrix(cbind(1, X_pred))
    
    v <- X_pred %*% beta
    
    
  }
  
  
  return(
    list(
      beta = beta,
      v = v,
      props = dat$pred.treatment,
      weights = dat$weights.trt,
      selected = selected,
      bad_props = bad_props
    )
  )
  
}

# overall function ----

#helper function to run all lasso types within a given rep
# rev - direction of rule, F for I(P0<theta), T for I(P0>=theta)
get_all_estimates <- function(dat,
                              Theta_vec,
                              knots_vec,
                              tailor_search_inds,
                              prop_adjust,
                              prop_adjust_inds,
                              outcome_adjust,
                              outcome_adjust_inds,
                              rev) {
  t_num <- length(tailor_search_inds)
  
  ests <- list()
  
  for(prop.selection in c(
    "no_selection",
    "oal",
    "glider",
    "dip",
    #"causal_ball",
    "hd_balancing"
  )){
  
  s_time <- Sys.time()  
    
  value <- rep(NA, t_num)
  
  runtime <- rep(NA,t_num)
  
  rule <- rep(NA, t_num)
  
  coefs <-
    matrix(NA, length(c(knots_vec, outcome_adjust_inds)) + 3, t_num)
  
  bad_props <- matrix(NA, ncol = 6, nrow = t_num)
  
  props <- matrix(NA, nrow = nrow(dat), ncol = t_num)
  
  weights <- matrix(NA, nrow = nrow(dat), ncol = t_num)
  
  selected <- matrix(NA, nrow = (ncol(dat) - 2), ncol = t_num)
  
  for (r in 1:t_num) {
    t_ind <- tailor_search_inds[r]
    #get the msm model
    res <- estimate.regime(
      dat,
      t_ind = t_ind,
      prop_adjust = prop_adjust,
      prop_adjust_inds = prop_adjust_inds,
      outcome_adjust = outcome_adjust,
      outcome_adjust_inds = outcome_adjust_inds,
      Theta_vec = Theta_vec,
      knots = knots_vec,
      prop.selection = prop.selection,
      rev = rev
    )
    
    
    #record the maximum value function
    value[r] <- max(res$v)
    
    rule[r] <- Theta_vec[which.max(res$v)]
    
    coefs[1:length(res$beta), r] <- res$beta
    
    props[, r] <- res$props
    
    weights[, r] <- res$weights
    
    bad_props[r,] <- res$bad_props
    
    selected[, r] <- res$selected
    
    runtime[r] <- (Sys.time() - s_time)
    
  }
  
  ests[[prop.selection]] <- list(
      value = value,
      rule = rule,
      coefs = coefs,
      props = props,
      weights = weights,
      selected = selected,
      bad_props = bad_props,
      runtime = runtime
    )
  
  print(prop.selection)
  
  gc()
  
  }
  
  return(ests)
}

# data cleaning and running the method ----

for(i in 1:25){

string <- paste0("G:/CTRHS/IMATS/Data/MHRN data/Analytic/Imputations/Final datasets/ADstartImp",i,".rdata")
load(string)

full_dat <- ADstartImp

A <- full_dat$InitMed

#putting names back to medications
medlevels <- readRDS("G:/CTRHS/IMATS/Nina/Code for ITR with confounding data example/medlevels.rds")
A <- sapply(A, function(x) medlevels[x])

SSRI <- c("CITALOPRAM", 
          "ESCITALOPRAM", 
          "FLUOXETINE", 
          "PAROXETINE", 
          "SERTRALINE", 
          "FLUVOXAMINE", 
          "VILAZODONE", 
          "VORTIOXETINE")

SNRI <- c("DESVENLAFAXINE", 
          "DULOXETINE", 
          "VENLAFAXINE")

#treatment is either SSRI or SNRI, ignoring others
A <- case_when(
  A %in% SSRI ~ 1,
  A %in% SNRI ~ 0,
  .default = NA
)


#reversing difference so we are maxing negative diff
Y <-  -1 * (full_dat$phq8_m6_d +
              as.numeric(full_dat$phq9_catm6) -
              full_dat$phq8_m1_d -
              as.numeric(full_dat$phq9_catm1) )
 

#combined non_anxiety diagnoses
full_dat <- full_dat %>% mutate(non_anx_dx_priorYr =
                                aud_dx_priorYr +
                                asd_dx_priorYr +
                                ptsd_dx_priorYr +
                                ocd_dx_priorYr +
                                oud_dx_priorYr +
                                per_dx_priorYr +
                                sed_dx_priorYr > 0,
                                #any psychthreapy in prior 5 years
                                psychtherapy_prior5yr = numMH_prior5yr > 0,
                                #basline phq score
                                baseline_phq = IndexTrtPHQ8_score + 
                                  as.numeric(phq9_cat),
                                baseline_bmi = WeightInitMed*0.453592/(HeightInitMed*0.0254)^2)

#the covariates being used
covars <- 
    #Age in years at treatment initiation was calculated
  c("AgeAtIndex",
    #information on patient sex (male or female). 
    "gend",
    #Health records data on race and ethnicity
    "race",
    #insurance type (commercial, Medicaid, Medicare, or private) 
    "ins",
    #neighborhood educational attainment (less than 25\% college degrees)
    "educ",
    #income (median  lower than 40,000 USD)
    "income",
    #level of poverty (20\% of households below federal poverty level)
    "pov",
    #urban or rural area (1 to 6, with 1 the most urban and 6 the most rural)
    "urban_rural",
    #Charlson score at initiation
    "Charlson",
    #tobacco use in the year prior to init
    "tobacco_priorYr",
    #BMI
    "baseline_bmi",
    #past year anxiety
    "anx_dx_priorYr",
    #past year non-anxiety mental health or substance use disorder
    "non_anx_dx_priorYr",
    #number of suicide attempts 6 months prior
    "num_SH_prior6m",
    #number of psych hospitalizations 6 months prior
    "num_MHIP_prior6m",
    #number different antidepressants taken 5 years prior
    "num_AD_Prior5yr",
    #received psychotherapy in the 5 years prior
    "psychtherapy_prior5yr",
    #number of PHQ measurements recorded year prior
    "PriorPhqCount",
    #PHQ recorded closest to treatment initiation
    "baseline_phq",
    #calendar year of treatment initiation
    "epiYear")

X <- full_dat[,covars]

f <- paste0(c("~",names(X)),
            c("",rep(" + ",length(names(X))-1),""),
            collapse = "")

#creating indicator variables for categorical covars
X <- model.matrix(as.formula(f), data = X)[,-1]

analysis_dat <- na.omit(data.frame(X = X, A = A, Y = Y))

#potential tailoring var initial phq score
Theta_vec <- seq(from = quantile(analysis_dat$X.baseline_phq,0.05),
                to = quantile(analysis_dat$X.baseline_phq,0.95),
                length.out = 20)
knots_vec <- Theta_vec[-c(1,length(Theta_vec))]

results_phq <- get_all_estimates(dat = analysis_dat,
                  Theta_vec,
                  knots_vec,
                  #initial phq index
                  tailor_search_inds = c(30),
                  prop_adjust = T,
                  prop_adjust_inds = 1:31,
                  outcome_adjust = "Lasso",
                  outcome_adjust_inds = 1:31,
                  rev = F)

results_phq_rev <- get_all_estimates(dat = analysis_dat,
                                 Theta_vec,
                                 knots_vec,
                                 #initial phq index
                                 tailor_search_inds = c(30),
                                 prop_adjust = T,
                                 prop_adjust_inds = 1:31,
                                 outcome_adjust = "Lasso",
                                 outcome_adjust_inds = 1:31,
                                 rev = T)

saveRDS(results_phq,paste0("results_phq",i,".rds"))

saveRDS(results_phq_rev,paste0("results_phq_rev",i,".rds"))


print(paste0("Imputation number ",i))

}

sink()
