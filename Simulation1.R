rm(list = ls())
library(rootSolve) 
library(mgcv)
library(tidyverse)
library(gee)
source("Simulation1_functions.R")


beta0 <- 1
beta1 <- -0.9
beta_true <- c(beta0, beta1)
nsim <- 1000
t_total <- 20
i_total_list <- c(20, 30, 50, 100, 200)
set.seed(123)

# 1. correct + correct -------------------------------------------------------
setting_num <- 1
r_model_method <- "gam"
r_model_formula <- as.formula("Y ~ s(X) + s(t)")
m_model_method <- "gam"
m_model_formula <- as.formula("A ~ s(X) + s(t)")
mu_model_method <- "gam"
mu_model_formula <- as.formula("Y ~ s(X) + s(t)") #this is the model formula for A = 0 and A = 1.
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("Y ~ s(X) + s(t)")
moderator <- "X"
#comparator formulas:
gee_model_formula <- as.formula("Y ~ A*X + t")
gam_model_formula <- as.formula("Y ~ A + A:X + s(X) + s(t)")

setting_list <- list(list(setting = paste0(setting_num, "1"), control_pattern = "linear", 
                          h1_a = -1, h1_b = 1, h1_c = -0.1, 
                          h2_a = -0.6, h2_b = -0.1, h2_c = 0.2),
                     list(setting = paste0(setting_num, "2"), control_pattern = "beta", 
                          h1_a = -0.5, h1_b = 1.1, h1_c = -1.2,
                          h2_a = -0.6, h2_b = -0.4, h2_c = 2),
                     list(setting = paste0(setting_num, "3"), control_pattern = "periodic", 
                          h1_a = -0.5, h1_b = 0.8, h1_c = -0.8,
                          h2_a = -0.2, h2_b = -0.4, h2_c = 2))

for (functype_i in 1:length(setting_list)){
  
  control_pattern <- setting_list[[functype_i]]$control_pattern
  h1_a <- setting_list[[functype_i]]$h1_a
  h1_b <- setting_list[[functype_i]]$h1_b
  h1_c <- setting_list[[functype_i]]$h1_c
  h2_a <- setting_list[[functype_i]]$h2_a
  h2_b <- setting_list[[functype_i]]$h2_b
  h2_c <- setting_list[[functype_i]]$h2_c
  
  result_i_total <- logistic_cee_sim(i_total_list = i_total_list, t_total = t_total,
                                     control_pattern = control_pattern,
                                     h1_a = h1_a, h1_b = h1_b, h1_c = h1_c, 
                                     h2_a = h2_a, h2_b = h2_b, h2_c = h2_c,
                                     beta0 = beta0, beta1 = beta1,
                                     moderator = moderator,
                                     r_model_method = r_model_method, r_model_formula = r_model_formula,
                                     m_model_method = m_model_method, m_model_formula = m_model_formula,
                                     mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                     assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula,
                                     gee_model_formula = gee_model_formula, gam_model_formula = gam_model_formula)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result_i_total, file_name)
  cat(file_name, "")
}

# 2. incorrect + correct -------------------------------------------------------
setting_num <- 2
r_model_method <- "gam"
r_model_formula <- as.formula("Y ~ s(X)")
m_model_method <- "gam"
m_model_formula <- as.formula("A ~ s(X) + s(t)")
mu_model_method <- "gam"
mu_model_formula <- as.formula("Y ~ s(X)") #this is the model formula for A = 0 and A = 1.
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("Y ~ s(X)")
moderator <- "X"
#comparator formulas:
gee_model_formula <- as.formula("Y ~ A*X")
gam_model_formula <- as.formula("Y ~ A + A:X + s(X)")

setting_list <- list(list(setting = paste0(setting_num, "1"), control_pattern = "linear", 
                          h1_a = -1, h1_b = 1, h1_c = -0.1, 
                          h2_a = -0.6, h2_b = -0.1, h2_c = 0.2),
                     list(setting = paste0(setting_num, "2"), control_pattern = "beta", 
                          h1_a = -0.5, h1_b = 1.1, h1_c = -1.2,
                          h2_a = -0.6, h2_b = -0.4, h2_c = 2),
                     list(setting = paste0(setting_num, "3"), control_pattern = "periodic", 
                          h1_a = -0.5, h1_b = 0.8, h1_c = -0.8,
                          h2_a = -0.2, h2_b = -0.4, h2_c = 2))


for (functype_i in 1:length(setting_list)){
  
  control_pattern <- setting_list[[functype_i]]$control_pattern
  h1_a <- setting_list[[functype_i]]$h1_a
  h1_b <- setting_list[[functype_i]]$h1_b
  h1_c <- setting_list[[functype_i]]$h1_c
  h2_a <- setting_list[[functype_i]]$h2_a
  h2_b <- setting_list[[functype_i]]$h2_b
  h2_c <- setting_list[[functype_i]]$h2_c
  
  result_i_total <- logistic_cee_sim(i_total_list = i_total_list, t_total = t_total,
                                     control_pattern = control_pattern,
                                     h1_a = h1_a, h1_b = h1_b, h1_c = h1_c, 
                                     h2_a = h2_a, h2_b = h2_b, h2_c = h2_c,
                                     beta0 = beta0, beta1 = beta1,
                                     moderator = moderator,
                                     r_model_method = r_model_method, r_model_formula = r_model_formula,
                                     m_model_method = m_model_method, m_model_formula = m_model_formula,
                                     mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                     assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula,
                                     gee_model_formula = gee_model_formula, gam_model_formula = gam_model_formula)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result_i_total, file_name)
  cat(file_name, "")
}




# 3. correct + incorrect -------------------------------------------------------
setting_num <- 3
r_model_method <- "gam"
r_model_formula <- as.formula("Y ~ s(X) + s(t)")
m_model_method <- "gam"
m_model_formula <- as.formula("A ~ s(X)")
mu_model_method <- "gam"
mu_model_formula <- as.formula("Y ~ s(X) + s(t)") #this is the model formula for A = 0 and A = 1.
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("Y ~ s(X) + s(t)")
moderator <- "X"
#comparator formulas:
gee_model_formula <- as.formula("Y ~ A*X")
gam_model_formula <- as.formula("Y ~ A + A:X + s(X)")


setting_list <- list(list(setting = paste0(setting_num, "1"), control_pattern = "linear", 
                          h1_a = -1, h1_b = 1, h1_c = -0.1, 
                          h2_a = -0.6, h2_b = -0.1, h2_c = 0.2),
                     list(setting = paste0(setting_num, "2"), control_pattern = "beta", 
                          h1_a = -0.5, h1_b = 1.1, h1_c = -1.2,
                          h2_a = -0.6, h2_b = -0.4, h2_c = 2),
                     list(setting = paste0(setting_num, "3"), control_pattern = "periodic", 
                          h1_a = -0.5, h1_b = 0.8, h1_c = -0.8,
                          h2_a = -0.2, h2_b = -0.4, h2_c = 2))


for (functype_i in 1:length(setting_list)){
  
  control_pattern <- setting_list[[functype_i]]$control_pattern
  h1_a <- setting_list[[functype_i]]$h1_a
  h1_b <- setting_list[[functype_i]]$h1_b
  h1_c <- setting_list[[functype_i]]$h1_c
  h2_a <- setting_list[[functype_i]]$h2_a
  h2_b <- setting_list[[functype_i]]$h2_b
  h2_c <- setting_list[[functype_i]]$h2_c
  
  result_i_total <- logistic_cee_sim(i_total_list = i_total_list, t_total = t_total,
                                     control_pattern = control_pattern,
                                     h1_a = h1_a, h1_b = h1_b, h1_c = h1_c, 
                                     h2_a = h2_a, h2_b = h2_b, h2_c = h2_c,
                                     beta0 = beta0, beta1 = beta1,
                                     moderator = moderator,
                                     r_model_method = r_model_method, r_model_formula = r_model_formula,
                                     m_model_method = m_model_method, m_model_formula = m_model_formula,
                                     mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                     assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula,
                                     gee_model_formula = gee_model_formula, gam_model_formula = gam_model_formula)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result_i_total, file_name)
  cat(file_name, "")
}



# 4. incorrect + incorrect -------------------------------------------------------
setting_num <- 4
r_model_method <- "gam"
r_model_formula <- as.formula("Y ~ s(X)")
m_model_method <- "gam"
m_model_formula <- as.formula("A ~ s(X)")
mu_model_method <- "gam"
mu_model_formula <- as.formula("Y ~ s(X)")  #this is the model formula for A = 0 and A = 1.
assoc_model_method <- "gam"
assoc_model_formula <- as.formula("Y ~ s(X)")
moderator <- "X"
#comparator formulas:
gee_model_formula <- as.formula("Y ~ A*X")
gam_model_formula <- as.formula("Y ~ A + A:X + s(X)")


setting_list <- list(list(setting = paste0(setting_num, "1"), control_pattern = "linear", 
                          h1_a = -1, h1_b = 1, h1_c = -0.1, 
                          h2_a = -0.6, h2_b = -0.1, h2_c = 0.2),
                     list(setting = paste0(setting_num, "2"), control_pattern = "beta", 
                          h1_a = -0.5, h1_b = 1.1, h1_c = -1.2,
                          h2_a = -0.6, h2_b = -0.4, h2_c = 2),
                     list(setting = paste0(setting_num, "3"), control_pattern = "periodic", 
                          h1_a = -0.5, h1_b = 0.8, h1_c = -0.8,
                          h2_a = -0.2, h2_b = -0.4, h2_c = 2))


for (functype_i in 1:length(setting_list)){
  
  control_pattern <- setting_list[[functype_i]]$control_pattern
  h1_a <- setting_list[[functype_i]]$h1_a
  h1_b <- setting_list[[functype_i]]$h1_b
  h1_c <- setting_list[[functype_i]]$h1_c
  h2_a <- setting_list[[functype_i]]$h2_a
  h2_b <- setting_list[[functype_i]]$h2_b
  h2_c <- setting_list[[functype_i]]$h2_c
  
  result_i_total <- logistic_cee_sim(i_total_list = i_total_list, t_total = t_total,
                                     control_pattern = control_pattern,
                                     h1_a = h1_a, h1_b = h1_b, h1_c = h1_c, 
                                     h2_a = h2_a, h2_b = h2_b, h2_c = h2_c,
                                     beta0 = beta0, beta1 = beta1,
                                     moderator = moderator,
                                     r_model_method = r_model_method, r_model_formula = r_model_formula,
                                     m_model_method = m_model_method, m_model_formula = m_model_formula,
                                     mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                     assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula,
                                     gee_model_formula = gee_model_formula, gam_model_formula = gam_model_formula)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result_i_total, file_name)
  cat(file_name, "")
}






