rm(list = ls())
library(rootSolve) 
library(mgcv)
library(tidyverse)
library(gee)
source("DataAnalysis_functions.R")

dat <- as.data.frame(readRDS("FINAL Dataset_A.rds"))
dat$Y <- dat$primary_outcome
dat$A <- dat$treatment
i_total <- length(unique(dat$ID))
prob_A <- 0.6
dat$prob_A <- prob_A

result_collect <- list()
## Marginal analysis ---------------------------------------------------------------
moderator <- c()
## Fit the model for r(s)
### r_model = r_t(S_t) = logit E(Y | S, A = 0) 
r_model_method <- "gam"
r_model_formula <- as.formula("primary_outcome ~ s(days_since_download)")

## Fit the model for m(s)
#### m_model = m(S) = E(A | S, Y = 0) 
m_model_method <- "gam"
m_model_formula <- as.formula("treatment ~ s(days_since_download)")

## Fit the model for mu
### mu(H, A) = E(Y | H, A)
mu_model_method <- "gam"
mu_model_formula <- as.formula("primary_outcome ~ s(days_since_download)")

## Fit the model for association model
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("primary_outcome ~ s(days_since_download)")


beta_hat <- beta_two_stage_nonparam(dat = dat, i_total = i_total, prob_A = prob_A,
                                    moderator = moderator,
                                    r_model_method = r_model_method, r_model_formula = r_model_formula,
                                    m_model_method = m_model_method, m_model_formula = m_model_formula,
                                    mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                    assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula)
beta_hat
## Moderated analysis (days since download) ---------------------------------------------------------------
moderator <- c("days_since_download")

## Fit the model for r(s)
### r_model = r_t(S_t) = logit E(Y | S, A = 0) 
r_model_method <- "gam"
r_model_formula <- as.formula("primary_outcome ~ s(days_since_download)")

## Fit the model for m(s)
#### m_model = m(S) = E(A | S, Y = 0) 
m_model_method <- "gam"
m_model_formula <- as.formula("treatment ~ s(days_since_download)")

## Fit the model for mu
### mu(H, A) = E(Y | H, A)
mu_model_method <- "gam"
mu_model_formula <- as.formula("primary_outcome ~ s(days_since_download)")

## Fit the model for association model
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("primary_outcome ~ s(days_since_download)")


beta_hat <- beta_two_stage_nonparam(dat = dat, i_total = i_total, prob_A = prob_A,
                                    moderator = moderator,
                                    r_model_method = r_model_method, r_model_formula = r_model_formula,
                                    m_model_method = m_model_method, m_model_formula = m_model_formula,
                                    mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                    assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula)
beta_hat
## Moderated analysis (before 8pm) ---------------------------------------------------------------
moderator <- c("before_8pm")

## Fit the model for r(s)
### r_model = r_t(S_t) = logit E(Y | S, A = 0) 
r_model_method <- "gam"
r_model_formula <- as.formula("primary_outcome ~ before_8pm + s(days_since_download)")

## Fit the model for m(s)
#### m_model = m(S) = E(A | S, Y = 0) 
m_model_method <- "gam"
m_model_formula <- as.formula("treatment ~ before_8pm + s(days_since_download)")

## Fit the model for mu
### mu(H, A) = E(Y | H, A)
mu_model_method <- "gam"
mu_model_formula <- as.formula("primary_outcome ~ before_8pm + s(days_since_download)")

## Fit the model for association model
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("primary_outcome ~ before_8pm + s(days_since_download)")


beta_hat <- beta_two_stage_nonparam(dat = dat, i_total = i_total, prob_A = prob_A,
                                    moderator = moderator,
                                    r_model_method = r_model_method, r_model_formula = r_model_formula,
                                    m_model_method = m_model_method, m_model_formula = m_model_formula,
                                    mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                    assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula)

beta_hat


## Moderated analysis (habituation) ---------------------------------------------------------------
moderator <- c("habituation")

## Fit the model for r(s)
### r_model = r_t(S_t) = logit E(Y | S, A = 0) 
r_model_method <- "gam"
r_model_formula <- as.formula("primary_outcome ~ habituation + s(days_since_download)")

## Fit the model for m(s)
#### m_model = m(S) = E(A | S, Y = 0) 
m_model_method <- "gam"
m_model_formula <- as.formula("treatment ~ habituation + s(days_since_download)")

## Fit the model for mu
### mu(H, A) = E(Y | H, A)
mu_model_method <- "gam"
mu_model_formula <- as.formula("primary_outcome ~ habituation + s(days_since_download)")

## Fit the model for association model
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("primary_outcome ~ habituation + s(days_since_download)")


beta_hat <- beta_two_stage_nonparam(dat = dat, i_total = i_total, prob_A = prob_A,
                                    moderator = moderator,
                                    r_model_method = r_model_method, r_model_formula = r_model_formula,
                                    m_model_method = m_model_method, m_model_formula = m_model_formula,
                                    mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                    assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula)
beta_hat