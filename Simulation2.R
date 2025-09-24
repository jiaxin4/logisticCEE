rm(list = ls())
library(rootSolve) 
library(mgcv)
library(tidyverse)
library(gee)
source("Simulation2_function.R")


nsim <- 1000
t_total <- 20
i_total_list <- c(20, 30, 50, 100, 200)
set.seed(123)
setting_list <- list(list(setting = 1, control_pattern = "linear", 
                          a_a = 2, b_a = -2, 
                          a_y1 = 0.8, b_y1 = -0.3, c_y1 = 0.1,
                          a_y0 = 0.1, b_y0 = 0.3, c_y0 = 0.1),
                     list(setting = 2, control_pattern = "beta", 
                          a_a = 2, b_a = -2, 
                          a_y1 = 0.4, b_y1 = 0.3, c_y1 = -0.1,
                          a_y0 = 0.7, b_y0 = -0.4, c_y0 = 0.1),
                     list(setting = 3, control_pattern = "periodic", 
                          a_a = 2, b_a = -2, 
                          a_y1 = 0.6, b_y1 = 0.1, c_y1 = -0.1,
                          a_y0 = 0.45, b_y0 = 0.1, c_y0 = 0.05))
for (functype_i in 1:3){
  control_pattern <- setting_list[[functype_i]]$control_pattern
  a_a <- setting_list[[functype_i]]$a_a
  b_a <- setting_list[[functype_i]]$b_a
  a_y1 <- setting_list[[functype_i]]$a_y1
  b_y1 <- setting_list[[functype_i]]$b_y1
  c_y1 <- setting_list[[functype_i]]$c_y1
  a_y0 <- setting_list[[functype_i]]$a_y0
  b_y0 <- setting_list[[functype_i]]$b_y0
  c_y0 <- setting_list[[functype_i]]$c_y0
  sim_num <- setting_list[[functype_i]]$setting
  control_pattern <- setting_list[[functype_i]]$control_pattern
  
  # 1. S_t is emptyset -------------------------------------------------------
  r_model_method <- "gam"
  r_model_formula <- as.formula("Y ~ s(t)")
  m_model_method <- "gam"
  m_model_formula <- as.formula("A ~ s(t)")
  mu_model_method <- "gam"
  mu_model_formula <- as.formula("Y ~ s(t)")
  assoc_model_method <- "gam"
  assoc_model_formula <- as.formula("Y ~ s(t)")
  moderator <- c()
  
  result_i_total <- logistic_cee_sim(i_total_list = i_total_list, t_total = t_total,
                                     control_pattern = control_pattern,
                                     moderator = moderator,
                                     r_model_method = r_model_method, r_model_formula = r_model_formula,
                                     m_model_method = m_model_method, m_model_formula = m_model_formula,
                                     mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                     assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula,
                                     a_a, b_a, a_y1, b_y1, c_y1, a_y0, b_y0, c_y0)
  file_name <- paste0("sim", sim_num, ".RDS")
  saveRDS(result_i_total, file_name)
  cat(file_name, "")
  
}

