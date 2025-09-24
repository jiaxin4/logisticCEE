# necessary functions
expit <- function(x){
  1/(1+exp(-x))
}

logit <- function(p){
  log(p/(1-p))
}


gen_linear_nonlinear_xt_t <- function(whichtype, xt, t, t_total, a, b, c){
  if (whichtype == "beta"){
    out <- a + b*dbeta(xt/2, 2, 2) + c*dbeta(t/t_total, 2, 2)
  } else if (whichtype == "periodic"){
    out <- a + b*sin(xt) + c*sin(t)
  } else if (whichtype == "linear"){
    out <- a + b*xt + c*t
  }
  out
}

gen_data_jointZY <- function(i_total, t_total, 
                             control_pattern, 
                             h1_a, h1_b, h1_c,
                             h2_a, h2_b, h2_c,
                             beta0, beta1){
  dat <- c()
  for (i in 1:i_total){
    t <- 1:t_total
    
    ## generate X
    Xi <- runif(t_total, 0, 2)
    
    ## generate h_1(X_t)
    h1_xt <- gen_linear_nonlinear_xt_t(whichtype = control_pattern, 
                                       xt = Xi, 
                                       t = t, t_total = t_total, 
                                       a = h1_a, b = h1_b, c = h1_c)
    ## generate h_2(X_t)
    h2_xt <- gen_linear_nonlinear_xt_t(whichtype = control_pattern, 
                                       xt = Xi, 
                                       t = t, t_total = t_total, 
                                       a = h2_a, b = h2_b, c = h2_c)
    ## generate A_t and Y_t jointly
    ### probability table
    alpha1 = 0.25
    alpha2 = -0.25
    prob_table <- matrix(NA, nrow = t_total, ncol = 4)
    prob_table[,1] <- 1                                       # A=0, Y=0
    prob_table[,2] <- exp(alpha2 + h2_xt)                      # A=0, Y=1
    prob_table[,3] <- exp(alpha1 + h1_xt)                    # A=1, Y=0
    prob_table[,4] <- exp(beta0 + beta1*Xi + alpha1 + alpha2 + h1_xt + h2_xt) # A=1, Y=1
    ### Normalize probabilities row-wise to be (0,1), or we can see this as using bayes rule to get conditinal probability
    prob_table <- prob_table / rowSums(prob_table)
    ### Sample (A, Y) for each observation
    labels <- matrix(c(0,0, 0,1, 1,0, 1,1), ncol = 2, byrow = TRUE) # create labels: (A,Y) = (0,0), (0,1), (1,0), (1,1)
    #### Sample (A, Y) pairs
    AY <- t(apply(prob_table, 1, function(p) labels[sample(1:4, 1, prob = p), ]))
    Ai <- AY[, 1]
    Yi <- AY[, 2]
    
    ## calculate P(A_t = 1|S_t) = P(At = 1, Yt = 0|St) + P(At = 1, Yt = 1|St)
    prob_A <- prob_table[,3] + prob_table[,4]
    ## calculate P(Y_t = 1|S_t, A_t)
    expect_Y <- expit((beta0 + beta1*Xi)*Ai + alpha2 + h2_xt)
    
    # The data
    dat <- rbind(dat, cbind(i, t, Ai, Xi, Yi, prob_A, expect_Y))
    
  }
  colnames(dat) <- c("id", "t", "A", "X", "Y", "prob_A", "expect_Y")
  dat <- as.data.frame(dat)
  dat
}


beta_two_stage_nonparam <- function(dat, 
                                    moderator,
                                    r_model_method, r_model_formula,
                                    m_model_method, m_model_formula,
                                    mu_model_method, mu_model_formula,
                                    assoc_model_method, assoc_model_formula,
                                    gee_model_formula, gam_model_formula){
  
  i_total <- length(unique(dat$id))
  t_total <- nrow(dat)/i_total
  Y <- dat$Y
  A <- dat$A
  X <- dat$X
  
  ######################## Stage 1
  ## Fit the model for r(s)
  ### r_model = r(S) = logit E(Y | S, A = 0) = X_r%*%gamma_r
  if (r_model_method == "glm") {
    r_model_fit <- glm(r_model_formula, family = binomial(link="logit"), data = dat, weights = 1 - A)
    r_predicted <- predict(r_model_fit, newdata = dat, type = "link")
  } else if (r_model_method == "gam") {
    r_model_fit <- gam(r_model_formula, data = dat, family = binomial(link = "logit"), weights = 1 - A)
    r_predicted <- predict(r_model_fit, newdata = dat, type = "link")
  } else if (r_model_method == "true") {
    r_predicted <- dat$r0_given_S
  }
  ## Fit the model for m(s)
  #### m_model = m(S) = E(A | S, Y = 0) 
  if (m_model_method == "glm") {
    m_model_fit <- glm(m_model_formula, family = binomial(link = "logit"), data = dat, weights = 1 - Y)
    m_predicted <- as.numeric(predict(m_model_fit, newdata = dat, type = "response"))
  } else if (m_model_method == "gam") {
    m_model_fit <- gam(m_model_formula, data = dat, family = binomial(link = "logit"), weights = 1 - Y)
    m_predicted <- predict(m_model_fit, newdata = dat, type = "response")
  } else if (m_model_method == "true") {
    m_predicted <- dat$m0_given_S
  } 
  
  ## Fit the model for mu
  ###additional nuisance regression model for efficiency improvement: mu(H, A) = E(Y | H, A)
  if (mu_model_method == "zero") {
    mu_a1 <- mu_a0 <- mu <- rep(0, length(A))
  } else if (mu_model_method == "gam") {
    dat_a0 <- dat[dat$A == 0, ]
    mu0_model_fit <- gam(mu_model_formula, data = dat_a0, family = binomial(link = "logit"))
    mu_a0 <- predict(mu0_model_fit, newdata = dat, type = "response")
    
    dat_a1 <- dat[dat$A == 1, ]
    mu1_model_fit <- gam(mu_model_formula, data = dat_a1, family = binomial(link = "logit"))
    mu_a1 <- predict(mu1_model_fit, newdata = dat, type = "response")
    
    mu <- ifelse(A, mu_a1, mu_a0)
  } else if (mu_model_method == "glm") {
    dat_a0 <- dat[dat$A == 0, ]
    mu0_model_fit <- glm(mu_model_formula, data = dat_a0, family = binomial(link = "logit"))
    mu_a0 <- predict(mu0_model_fit, newdata = dat, type = "response")
    
    dat_a1 <- dat[dat$A == 1, ]
    mu1_model_fit <- glm(mu_model_formula, data = dat_a1, family = binomial(link = "logit"))
    mu_a1 <- predict(mu1_model_fit, newdata = dat, type = "response")
    mu <- ifelse(A, mu_a1, mu_a0)
  }
  
  ## Fit the model for association model
  ### do not do the improved effciency stuff for the estimation equation of nuisance association parameter.
  if (assoc_model_method == "gam"){
    fit_assoc_model <- gam(assoc_model_formula, data = dat, family = binomial(link = "logit"), weights = A / prob_A)
    a_predicted <- predict(fit_assoc_model, newdata = dat, type = "link")
  } else if (assoc_model_method == "glm"){
    fit_assoc_model <- glm(assoc_model_formula, data = dat, family = binomial(link = "logit"), weights = A / prob_A)
    a_predicted <- predict(fit_assoc_model, newdata = dat, type = "link")
  }
  
  ### save the predicted nuisance parameters
  dat$r_predicted <- r_predicted
  dat$m_predicted <- m_predicted
  dat$mu <- mu
  dat$mu_a1 <- mu_a1
  dat$mu_a0 <- mu_a0
  dat$a_predicted <- a_predicted
  
  ########################## Stage 2
  #### Error capturing function
  run_method <- function(method, var_method, ...) {
    tryCatch(
      {
        beta_hat <- method(dat = dat, moderator = moderator)
        beta_var <- var_method(dat = dat, beta = beta_hat, moderator = moderator, ...)
        list(beta_hat = beta_hat, beta_var = beta_var)
      },
      error = function(cond) {
        message("\nCaught error in ", deparse(substitute(method)))
        message(cond)
        S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
        list(beta_hat = rep(NaN, ncol(S_mat)),
             beta_var = rep(NaN, ncol(S_mat)))
      }
    )
  }
  
  # Compute for Method1
  beta_hat_var_m1 <- run_method(beta_est_method1, beta_var_method1_nonparam)
  
  # Compute for Method2 (include optional arguments)
  beta_hat_var_m2 <- run_method(beta_est_method2, beta_var_method2_nonparam,
                                fit_assoc_model = fit_assoc_model)
  
  ## logisticGEE
  logisticGEE_model <- gee(gee_model_formula, data = dat, id = id, family = binomial, corstr = "independence")
  beta_hat_logisticGEE <- as.numeric(logisticGEE_model$coefficients[c("A", "A:X")])
  beta_var_logisticGEE <- as.numeric(diag(logisticGEE_model$robust.variance)[c("A", "A:X")])
  
  ## logistcGAM
  logisticGAM_model <- gam(gam_model_formula, data = dat, family = binomial(link = "logit"))
  beta_hat_logisticGAM <- as.numeric(summary(logisticGAM_model)$p.coeff[c("A", "A:X")])
  beta_var_logisticGAM <- (as.numeric(summary(logisticGAM_model)$se[c("A", "A:X")]))^2
  
  # Combine results into a list
  beta_est <- c(beta_hat_var_m1, 
                beta_hat_var_m2,
                list(beta_hat_logisticGEE, beta_var_logisticGEE),
                list(beta_hat_logisticGAM, beta_var_logisticGAM))
  names(beta_est) <- c("beta_hat_m1", "beta_var_m1",
                       "beta_hat_m2", "beta_var_m2",
                       "beta_hat_logisticGEE", "beta_var_logisticGEE",
                       "beta_hat_logisticGAM", "beta_var_logisticGAM")
  beta_est
}


## Method 1 beta_hat and var estimator -------------------------------------
beta_est_method1 <- function(dat, moderator){
  i_total <- length(unique(dat$id))
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))
  prob_A <- dat$prob_A
  dim_beta <- ncol(S_mat)

  ee <- function(beta) {
    Sbeta <- as.numeric(S_mat %*% beta)
    exp_r <- exp(dat$r_predicted)
    U1 <- (dat$Y - dat$mu) * (exp(- dat$A * Sbeta) + exp_r) * (dat$A - dat$m_predicted)
    U2 <- (dat$mu_a1 * exp(-Sbeta) - (1 - dat$mu_a1) * exp_r) * (1 - dat$m_predicted) * prob_A
    U3 <- - (dat$mu_a0 - (1 - dat$mu_a0) * exp_r) * dat$m_predicted * (1 - prob_A)
    
    U <- U1 + U2 + U3
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (ibeta in 1:dim_beta) {
      ef[ibeta] <- sum( U * S_mat[, ibeta])
    }
    ef <- ef / i_total
    return(ef)
  }
  
  solution <-   multiroot(ee, rep(0, dim_beta), useFortran = FALSE)
  beta_hat_solve <- solution$root
  
  beta_hat_solve
}



beta_var_method1_nonparam <- function(dat, beta, moderator){
  i_total <- length(unique(dat$id))
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))
  beta <- as.matrix(beta, nrow = ncol(S_mat))
  Sbeta <- as.numeric(S_mat %*%beta) #ft*beta
  prob_A <- dat$prob_A
  dim_beta <- ncol(S_mat)
  
  total_person_decisionpoint <- nrow(dat)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, dim_beta, dim_beta))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, dim_beta, dim_beta))
  
  for (it in 1:total_person_decisionpoint){
    dat_it <- dat[it, ]
    fx <- t(t(S_mat[it, ])) #p by 1
    
    #calculate derivative of Utilde wrt parameter of interest
    ##w.r.t beta
    dee_part1 <- (dat_it$Y - dat_it$mu)*exp(- dat_it$A*Sbeta[it])*(dat_it$A - dat_it$m_predicted)*(-dat_it$A)*fx%*%t(fx)
    dee_part2 <- - dat_it$mu_a1*exp(-Sbeta[it])*(1-dat_it$m_predicted)*dat_it$prob_A*fx%*%t(fx)
    dee_Utilde_beta <- dee_part1 + dee_part2
    dee_sum[it, , ] <- dee_Utilde_beta
    
    #calculate meat
    exp_r <- exp(dat_it$r_predicted)
    U1 <- (dat_it$Y - dat_it$mu) * (exp(- dat_it$A * Sbeta[it]) + exp_r) * (dat_it$A - dat_it$m_predicted)
    U2 <- (dat_it$mu_a1 * exp(-Sbeta[it]) - (1 - dat_it$mu_a1) * exp_r) * (1 - dat_it$m_predicted) * dat_it$prob_A
    U3 <- - (dat_it$mu_a0 - (1 - dat_it$mu_a0) * exp_r) * dat_it$m_predicted * (1 - dat_it$prob_A)
    U_it <- (U1 + U2 + U3)*fx
    meat_sum[it, ,] <- U_it%*%t(U_it)
  }
  dee <- solve(apply(dee_sum, c(2,3), sum)/i_total)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee%*%meat%*%t(dee)/i_total
  asy_var <- diag(var_cov)
  
  asy_var
}

## Method 2 beta_hat and var estimator -------------------------------------

beta_est_method2 <- function(dat, moderator){
  i_total <- length(unique(dat$id))
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))
  prob_A <- dat$prob_A
  dim_beta <- ncol(S_mat)
  
  ee <- function(beta) {
    tmp <- (expit(dat$a_predicted - as.numeric(S_mat %*% beta)) - dat$mu_a0) -
      (1 - dat$A) / (1 - dat$prob_A) * (dat$Y - dat$mu_a0)
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (ibeta in 1:dim_beta) {
      ef[ibeta] <- sum(tmp * S_mat[, ibeta])
    }
    
    ef <- ef / i_total
    return(ef)
  }
  
  solution <-   multiroot(ee, rep(0, dim_beta), useFortran = FALSE)
  beta_hat_solve <- solution$root
  
  beta_hat_solve
}


beta_var_method2_nonparam <- function(dat, beta, moderator, fit_assoc_model){
  i_total <- length(unique(dat$id))
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))
  beta <- as.matrix(beta, nrow = ncol(S_mat))
  Sbeta <- as.numeric(S_mat %*%beta) #ft*beta
  dim_beta <- ncol(S_mat)
  
  total_person_decisionpoint <- nrow(dat)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, dim_beta, dim_beta))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, dim_beta, dim_beta))
  
  #For alpha
  X_alpha <- predict(fit_assoc_model, newdata = dat, type = "lpmatrix")
  alpha_hat <- coef(fit_assoc_model)
  #max(as.numeric(X_alpha %*% t(t(alpha_fit))) - predict(fit_assoc_model, newdata = dat, type = "link"))
  
  p <- ncol(X_alpha) + length(beta)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  for (it in 1:total_person_decisionpoint){
    dat_it <- dat[it, ]
    fx <- t(t(S_mat[it, ])) 
    
    #calculate derivative of Utilde wrt nuisance parameter and parameter of interest
    ##w.r.t beta
    dee_part1 <- (-1)/((1 + exp(-(dat_it$a_predicted -  Sbeta[it])))^2)
    dee_part2 <- exp(-(dat_it$a_predicted -  Sbeta[it]))
    dee_Utilde_beta <- (dee_part1*dee_part2)*fx%*%t(fx)
    ##w.r.t alpha
    dee_Utilde_alpha <- (dee_part1*dee_part2)*fx%*%(-X_alpha[it, ])
    
    #calcualte derivative of nuisance parameter's own estimating equations
    ## dee of alpha
    dee_alpha <- -(dat_it$A/dat_it$prob_A)*(1/(1 + exp(-dat_it$a_predicted))^2)*exp(-dat_it$a_predicted)*X_alpha[it, ]%*%t(X_alpha[it, ])
    
    #combine them
    max_rows <- max(nrow(dee_Utilde_beta), nrow(dee_Utilde_alpha))
    pad_matrix <- function(mat, max_rows) {
      nrow_diff <- max_rows - nrow(mat)
      if (nrow_diff > 0) {
        zero_pad <- matrix(0, nrow = nrow_diff, ncol = ncol(mat))
        mat <- rbind(mat, zero_pad)
      }
      return(mat)
    }
    dee_Utilde_beta <- pad_matrix(dee_Utilde_beta, max_rows)
    dee_Utilde_alpha <- pad_matrix(dee_Utilde_alpha, max_rows)
    combined_matrix_bottom <- cbind(dee_Utilde_alpha, dee_Utilde_beta)
    combined_matrix_diag <- dee_alpha
    dee_sum[it, , ] <- rbind(cbind(combined_matrix_diag, matrix(0, 
                                                                ncol = ncol(combined_matrix_bottom)- ncol(combined_matrix_diag),
                                                                nrow = nrow(combined_matrix_diag))
    ), combined_matrix_bottom)
    
    #calculate meat
    ## EE of alpha
    psi_alpha <- (dat_it$A/dat_it$prob_A)*(dat_it$Y - expit(dat_it$a_predicted))*t(t(X_alpha[it,]))
    ## EE of beta
    U1 <- expit(dat_it$a_predicted - Sbeta[it]) - dat_it$mu_a0
    U2 <- (1 - dat_it$A)/(1 - dat_it$prob_A)*(dat_it$Y - dat_it$mu_a0)
    U_it <- (U1 - U2)*fx
    ## combine them
    meat_sum[it, ,] <- rbind(psi_alpha, U_it)%*%t(rbind(psi_alpha, U_it))
    
  }
  dee <- apply(dee_sum, c(2,3), sum)/i_total
  dee_inv <- solve(dee)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee_inv%*%meat%*%t(dee_inv)/i_total
  asy_var <- diag(var_cov)
  asy_var_beta <- asy_var[(ncol(X_alpha)+1):p]
  
  asy_var_beta
}



## Perform simulations -------------------------------------

eval_ci <- function(beta_hat, beta_var, beta_true){
  ci_low <- beta_hat - 1.96*sqrt(beta_var)
  ci_high <- beta_hat + 1.96*sqrt(beta_var)
  result <- c()
  for (beta_i in 1:length(beta_hat)){
    beta_i_in <- ifelse(beta_true[beta_i] <= ci_high[beta_i] & beta_true[beta_i] >= ci_low[beta_i], 1, 0)
    result <- c(result, beta_i_in)
  }
  result
}



logistic_cee_sim <- function(i_total_list, t_total, 
                             control_pattern,
                             h1_a, h1_b, h1_c,
                             h2_a, h2_b, h2_c,
                             beta0, beta1,
                             moderator,
                             r_model_method, r_model_formula,
                             m_model_method, m_model_formula,
                             mu_model_method, mu_model_formula,
                             assoc_model_method, assoc_model_formula,
                             gee_model_formula, gam_model_formula){
  beta_true <- c(beta0, beta1)
  result_allsim <- list()
  result_i_total <- list()
  for (i_total in i_total_list){
    
    for (isim in 1:nsim){
      if (isim %% 10 == 0) {
        cat(isim, "")
      }
      dat <- gen_data_jointZY(i_total = i_total, t_total = t_total, 
                      control_pattern = control_pattern, 
                      h1_a = h1_a, h1_b = h1_b, h1_c = h1_c,
                      h2_a = h2_a, h2_b = h2_b, h2_c = h2_c,
                      beta0 = beta0, beta1 = beta1)
      
      beta_hat <- beta_two_stage_nonparam(dat = dat, 
                                 moderator = moderator,
                                 r_model_method = r_model_method, r_model_formula = r_model_formula,
                                 m_model_method = m_model_method, m_model_formula = m_model_formula,
                                 mu_model_method = mu_model_method, mu_model_formula = mu_model_formula,
                                 assoc_model_method = assoc_model_method, assoc_model_formula = assoc_model_formula,
                                 gee_model_formula = gee_model_formula, gam_model_formula = gam_model_formula)
      
      ci_in_m1 <- eval_ci(beta_hat$beta_hat_m1, beta_hat$beta_var_m1, beta_true = beta_true)
      ci_in_m2 <- eval_ci(beta_hat$beta_hat_m2, beta_hat$beta_var_m2, beta_true = beta_true)
      ci_in_logisticGEE <- eval_ci(beta_hat$beta_hat_logisticGEE, beta_hat$beta_var_logisticGEE, beta_true = beta_true)
      ci_in_logisticGAM <- eval_ci(beta_hat$beta_hat_logisticGAM, beta_hat$beta_var_logisticGAM, beta_true = beta_true)
      
      ci_in <- list(ci_in_m1, ci_in_m2, ci_in_logisticGEE, ci_in_logisticGAM)
      names(ci_in) <- c("ci_in_m1", "ci_in_m2", "ci_in_logisticGEE", "ci_in_logisticGAM")
      
      result <- list(c(beta_hat, ci_in),
                     expect_Y = mean(dat$expect_Y),
                     i_total = i_total, t_total = t_total,
                     model_formulas = list(r_model_formula = r_model_formula, 
                                           m_model_formula = m_model_formula, 
                                           mu_model_formula = mu_model_formula, 
                                           assoc_model_formula = assoc_model_formula, 
                                           gee_model_formula = gee_model_formula, 
                                           gam_model_formula = gam_model_formula))
      result_allsim[[isim]] <- result
    }
    result_i_total <- c(result_i_total, list(result_allsim))
  }
  result_i_total
}


















