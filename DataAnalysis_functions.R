# necessary functions
expit <- function(x){
  1/(1+exp(-x))
}

logit <- function(p){
  log(p/(1-p))
}

beta_two_stage_nonparam <- function(dat, i_total, prob_A,
                                    moderator,
                                    r_model_method, r_model_formula,
                                    m_model_method, m_model_formula,
                                    mu_model_method, mu_model_formula,
                                    assoc_model_method, assoc_model_formula){
  
  Y <- dat$Y
  A <- dat$A
  
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
        #beta_hat <- c(1, -0.5)
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
  beta_hat_var_m1 <- run_method(beta_est_method1, beta_var_method1_nonparam,
                                r_model_fit, m_model_fit, mu0_model_fit, mu1_model_fit)
  
  # Compute for Method2 (include optional arguments)
  beta_hat_var_m2 <- run_method(beta_est_method2, beta_var_method2_nonparam,
                                fit_assoc_model, mu0_model_fit)
  
  # Combine results into a list
  beta_est <- c(beta_hat_var_m1, 
                beta_hat_var_m2)
  names(beta_est) <- c("beta_hat_m1", "beta_var_m1",
                       "beta_hat_m2", "beta_var_m2")
  beta_est
}


## Method 1 beta_hat and var estimator -------------------------------------
beta_est_method1 <- function(dat, moderator){
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))
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

MatrixDiag <- function(A, B){
  new_rows <- nrow(A) + nrow(B)
  new_cols <- ncol(A) + ncol(B)
  # Create an empty matrix filled with zeros
  result_matrix <- matrix(0, nrow = new_rows, ncol = new_cols)
  # Place A in the top-left corner
  result_matrix[1:nrow(A), 1:ncol(A)] <- A
  # Place B in the bottom-right corner
  result_matrix[(nrow(A) + 1):new_rows, (ncol(A) + 1):new_cols] <- B
  result_matrix
}
pad_matrix <- function(mat, max_rows) {
  nrow_diff <- max_rows - nrow(mat)
  if (nrow_diff > 0) {
    zero_pad <- matrix(0, nrow = nrow_diff, ncol = ncol(mat))
    mat <- rbind(mat, zero_pad)
  }
  return(mat)
}

beta_var_method1_nonparam <- function(dat, beta, moderator, 
                                      r_model_fit, m_model_fit, mu0_model_fit, mu1_model_fit){
  
  S_mat <- as.matrix(cbind( rep(1, nrow(dat)), dat[, moderator]))
  Sbeta <- as.numeric(S_mat %*%beta) #ft*beta
  dim_beta <- length(beta)
  
  # For mu1
  X_mu1 <- predict(mu1_model_fit, newdata = dat, type = "lpmatrix")
  gamma_mu1 <- coef(mu1_model_fit)
  # For mu0
  X_mu0 <- predict(mu0_model_fit, newdata = dat, type = "lpmatrix")
  gamma_mu0 <- coef(mu0_model_fit)
  # For r
  X_r <- predict(r_model_fit, newdata = dat, type = "lpmatrix")
  gamma_r <- coef(r_model_fit)
  # For m
  X_m <- predict(m_model_fit, newdata = dat, type = "lpmatrix")
  gamma_m <- coef(m_model_fit)
  
  total_person_decisionpoint <- nrow(dat)
  nuisance_p <- ncol(X_mu1) + ncol(X_mu0) + ncol(X_r) + ncol(X_m)
  p <- nuisance_p + length(beta)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  
  for (it in 1:total_person_decisionpoint){
    dat_it <- dat[it, ]
    fx <- t(t(S_mat[it, ])) 
    
    #calculate derivative of Utilde wrt nuisance parameter and parameter of interest
    ## w.r.t beta
    dee_part1 <- (dat_it$Y - dat_it$mu)*exp(- dat_it$A*Sbeta[it])*(dat_it$A - dat_it$m_predicted)*(-dat_it$A)*fx%*%t(fx)
    dee_part2 <- - dat_it$mu_a1*exp(-Sbeta[it])*(1-dat_it$m_predicted)*dat_it$prob_A*fx%*%t(fx)
    dee_Utilde_beta <- dee_part1 + dee_part2
    
    ## w.r.t gamma_mu1
    dee_mu1_gamma_mu1 <- as.numeric(exp(-gamma_mu1%*%X_mu1[it,])/(1 + exp(-gamma_mu1%*%X_mu1[it,]))^2)*X_mu1[it,]
    dee_Utilde_gamma_mu1 <- -dat_it$A*(exp(-dat_it$A*Sbeta[it]) + exp(dat_it$r_predicted))*(dat_it$A - dat_it$m_predicted) + 
      (exp(-Sbeta[it]) + exp(dat_it$r_predicted))*(1 - dat_it$m_predicted)*dat_it$prob_A
    dee_Utilde_gamma_mu1 <- dee_Utilde_gamma_mu1*fx%*%dee_mu1_gamma_mu1
    
    ## w.r.t gamma_mu0
    dee_mu0_gamma_mu0 <- as.numeric(exp(-gamma_mu0%*%X_mu0[it,])/(1 + exp(-gamma_mu0%*%X_mu0[it,]))^2)*X_mu0[it,]
    dee_Utilde_gamma_mu0 <- -(1 - dat_it$A)*(exp(-dat_it$A*Sbeta[it]) + exp(dat_it$r_predicted))*(dat_it$A - dat_it$m_predicted) -
      (1 + exp(dat_it$r_predicted))*dat_it$m_predicted*(1 - dat_it$prob_A)
    dee_Utilde_gamma_mu0 <- dee_Utilde_gamma_mu0*fx%*%dee_mu0_gamma_mu0
    
    ## w.r.t gamma_r
    dee_Utilde_gamma_r <- (dat_it$Y - dat_it$A*dat_it$mu_a1 - (1 - dat_it$A)*dat_it$mu_a0)*exp(dat_it$r_predicted)*(dat_it$A - dat_it$m_predicted) -
      (1 - dat_it$mu_a1)*exp(dat_it$r_predicted)*(1 - dat_it$m_predicted)*dat_it$prob_A +
      (1 - dat_it$mu_a0)*exp(dat_it$r_predicted)*dat_it$m_predicted*(1 - dat_it$prob_A)
    dee_Utilde_gamma_r <- dee_Utilde_gamma_r*fx%*%X_r[it,]
    
    ## w.r.t gamma_m
    dee_m_gamma_m <- as.numeric(exp(-gamma_m%*%X_m[it,])/(1 + exp(-gamma_m%*%X_m[it,]))^2)*X_m[it,]
    dee_Utilde_gamma_m <- (dat_it$Y - dat_it$mu)*(exp(-dat_it$A*Sbeta[it]) + exp(dat_it$r_predicted)) + 
      (dat_it$mu_a1*exp(-Sbeta[it]) - (1 - dat_it$mu_a1)*exp(dat_it$r_predicted))*dat_it$prob_A +
      (dat_it$mu_a0 - (1 - dat_it$mu_a0)*exp(dat_it$r_predicted))*(1 - dat_it$prob_A)
    dee_Utilde_gamma_m <- -dee_Utilde_gamma_m*fx%*%dee_m_gamma_m
    
    #calcualte derivative of nuisance parameter's own estimating equations
    ##dee of mu1
    dee_mu1 <- -as.numeric(exp(-gamma_mu1%*%X_mu1[it,])/(1 + exp(-gamma_mu1%*%X_mu1[it,]))^2)*X_mu1[it,]%*%t(X_mu1[it,])
    
    ##dee of mu0 
    dee_mu0 <- -as.numeric(exp(-gamma_mu0%*%X_mu0[it,])/(1 + exp(-gamma_mu0%*%X_mu0[it,]))^2)*X_mu0[it,]%*%t(X_mu0[it,])
    
    ##dee of r
    dee_r <- -as.numeric(exp(-gamma_r%*%X_r[it,])/(1 + exp(-gamma_r%*%X_r[it,]))^2)*X_r[it,]%*%t(X_r[it,])
    
    ##dee of m
    dee_m <- -as.numeric(exp(-gamma_m%*%X_m[it,])/(1 + exp(-gamma_m%*%X_m[it,]))^2)*X_m[it,]%*%t(X_m[it,])
    
    #combine them
    max_rows <- max(nrow(dee_Utilde_beta), nrow(dee_Utilde_gamma_mu1), nrow(dee_Utilde_gamma_mu0), nrow(dee_Utilde_gamma_r), nrow(dee_Utilde_gamma_m))
    dee_Utilde_beta <- pad_matrix(dee_Utilde_beta, max_rows)
    dee_Utilde_gamma_mu1 <- pad_matrix(dee_Utilde_gamma_mu1, max_rows)
    dee_Utilde_gamma_mu0 <- pad_matrix(dee_Utilde_gamma_mu0, max_rows)
    dee_Utilde_gamma_r <- pad_matrix(dee_Utilde_gamma_r, max_rows)
    dee_Utilde_gamma_m <- pad_matrix(dee_Utilde_gamma_m, max_rows)
    combined_matrix_bottom <- cbind(dee_Utilde_gamma_mu1, dee_Utilde_gamma_mu0, dee_Utilde_gamma_r, dee_Utilde_gamma_m, dee_Utilde_beta)
    combined_matrix_diag <- MatrixDiag(MatrixDiag(MatrixDiag(dee_mu1, dee_mu0), dee_r), dee_m)
    
    dee_sum[it, , ] <- rbind(cbind(combined_matrix_diag, matrix(0, 
                                                                ncol = ncol(combined_matrix_bottom)- ncol(combined_matrix_diag),
                                                                nrow = nrow(combined_matrix_diag))
    ), combined_matrix_bottom)
    
    #calculate meat
    ## EE of gamma_mu1
    psi_gamma_mu1 <-  (dat_it$Y - dat_it$mu_a1)*t(t(X_mu1[it,]))
    
    ## EE of gamma_mu0
    psi_gamma_mu0 <-  (dat_it$Y - dat_it$mu_a0)*t(t(X_mu0[it,]))
    
    ## EE of gamma_r
    psi_gamma_r <-  (dat_it$Y - expit(dat_it$r_predicted))*t(t(X_r[it,])) 
    
    ## EE of gamma_m
    psi_gamma_m <-  (dat_it$A - dat_it$m_predicted)*t(t(X_m[it,]))
    
    ## EE of beta
    exp_r <- exp(dat_it$r_predicted)
    U1 <- (dat_it$Y - dat_it$mu) * (exp(- dat_it$A * Sbeta[it]) + exp_r) * (dat_it$A - dat_it$m_predicted)
    U2 <- (dat_it$mu_a1 * exp(-Sbeta[it]) - (1 - dat_it$mu_a1) * exp_r) * (1 - dat_it$m_predicted) * dat_it$prob_A
    U3 <- - (dat_it$mu_a0 - (1 - dat_it$mu_a0) * exp_r) * dat_it$m_predicted * (1 - dat_it$prob_A)
    U_it <- (U1 + U2 + U3)*fx
    meat_sum[it, ,] <- rbind(psi_gamma_mu1, psi_gamma_mu0, psi_gamma_r, psi_gamma_m, U_it)%*%t(rbind(psi_gamma_mu1, psi_gamma_mu0, psi_gamma_r, psi_gamma_m, U_it))
  }
  
  dee <- solve(apply(dee_sum, c(2,3), sum)/i_total)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee%*%meat%*%t(dee)/i_total
  asy_var <- diag(var_cov)
  asy_var_beta <- asy_var[(nuisance_p + 1):p]
  
  asy_var_beta
}



## Method 2 beta_hat and var estimator -------------------------------------

beta_est_method2 <- function(dat, moderator){
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))
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


beta_var_method2_nonparam <- function(dat, beta, moderator, 
                                      fit_assoc_model, mu0_model_fit){
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
  
  #For mu0
  X_mu0 <- predict(mu0_model_fit, newdata = dat, type = "lpmatrix")
  gamma_mu0 <- coef(mu0_model_fit)
  
  nuisance_p <- ncol(X_alpha) + ncol(X_mu0)
  p <- nuisance_p + length(beta)
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
    ##w.r.t mu0
    dee_mu0_gamma_mu0 <- as.numeric(exp(-gamma_mu0%*%X_mu0[it,])/(1 + exp(-gamma_mu0%*%X_mu0[it,]))^2)*X_mu0[it,]
    dee_Utilde_gamma_mu0 <- (-1 + (1 - dat_it$A)/(1 - dat_it$prob_A))*fx%*%dee_mu0_gamma_mu0
    
    #calcualte derivative of nuisance parameter's own estimating equations
    ## dee of alpha
    dee_alpha <- -(dat_it$A/dat_it$prob_A)*(1/(1 + exp(-dat_it$a_predicted))^2)*exp(-dat_it$a_predicted)*X_alpha[it, ]%*%t(X_alpha[it, ])
    
    ## dee of mu0
    dee_mu0 <- -as.numeric(exp(-gamma_mu0%*%X_mu0[it,])/(1 + exp(-gamma_mu0%*%X_mu0[it,]))^2)*X_mu0[it,]%*%t(X_mu0[it,])
    
    #combine them
    max_rows <- max(nrow(dee_Utilde_beta), nrow(dee_Utilde_alpha), nrow(dee_Utilde_gamma_mu0))
    
    dee_Utilde_beta <- pad_matrix(dee_Utilde_beta, max_rows)
    dee_Utilde_alpha <- pad_matrix(dee_Utilde_alpha, max_rows)
    dee_Utilde_gamma_mu0 <- pad_matrix(dee_Utilde_gamma_mu0, max_rows)
    combined_matrix_bottom <- cbind(dee_Utilde_alpha, dee_Utilde_gamma_mu0, dee_Utilde_beta)
    combined_matrix_diag <- MatrixDiag(dee_alpha, dee_mu0)
    dee_sum[it, , ] <- rbind(cbind(combined_matrix_diag, matrix(0, 
                                                                ncol = ncol(combined_matrix_bottom)- ncol(combined_matrix_diag),
                                                                nrow = nrow(combined_matrix_diag))
    ), combined_matrix_bottom)
    
    #calculate meat
    ## EE of alpha
    psi_alpha <- (dat_it$A/dat_it$prob_A)*(dat_it$Y - expit(dat_it$a_predicted))*t(t(X_alpha[it,]))
    ## EE of gamma_mu0
    psi_gamma_mu0 <-  (dat_it$Y - dat_it$mu_a0)*t(t(X_mu0[it,]))
    ## EE of beta
    U1 <- expit(dat_it$a_predicted - Sbeta[it]) - dat_it$mu_a0
    U2 <- (1 - dat_it$A)/(1 - dat_it$prob_A)*(dat_it$Y - dat_it$mu_a0)
    U_it <- (U1 - U2)*fx
    ## combine them
    meat_sum[it, ,] <- rbind(psi_alpha, psi_gamma_mu0, U_it)%*%t(rbind(psi_alpha, psi_gamma_mu0, U_it))
    
  }
  dee <- apply(dee_sum, c(2,3), sum)/i_total
  dee_inv <- solve(dee)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee_inv%*%meat%*%t(dee_inv)/i_total
  asy_var <- diag(var_cov)
  asy_var_beta <- asy_var[(nuisance_p+1):p]
  
  asy_var_beta
}

