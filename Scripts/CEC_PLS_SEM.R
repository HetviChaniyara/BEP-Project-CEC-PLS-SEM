# Cardinality and Equality Constrained PlS-SEM
# Hetvi Chaniyara
# Bachelor End Project


current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# install.packages("MASS") 
library(MASS)


Initialize_parameters <- function(X, R) {
  J <- dim(X)[2] # number of columns
  I <- dim(X)[1] # number of rows
  svd_X <- svd(X)
  
  # SVD-based components
  W_svd <- svd_X$v[, 1:R]
  P_svd <- t(svd_X$u[, 1:R])
  
  # Random components
  W_rand <- matrix(rnorm(length(W_svd), mean = 0, sd = 1), nrow = nrow(W_svd))
  P_rand <- matrix(rnorm(length(P_svd), mean = 0, sd = 1), nrow = nrow(P_svd))
  
  # Weighted combination: 0.7 * SVD + 0.3 * random
  W0 <- 0.7 * W_svd + 0.3 * W_rand
  P0_T <- 0.7 * P_svd + 0.3 * P_rand
  
  U <- matrix(0, nrow = J, ncol = R) # Initialize to 0
  rho <- 1 # penalty parameter
  alpha <- max(eigen(t(X) %*% X)$values) # max eigenvalue of X^TX
  
  return(list(W0 = W0, P0_T = P0_T, U = U, rho = rho, alpha = alpha))
}


CEC_PLS_SEM <-function(X, R, epsilon, phi){
  J = dim(X)[2] # number of columns
  I = dim(X)[1] # number of rows
  iter <- 0
  convAO <- 0
  MaxIter <- 200
  
  # Get initialized parameters
  params <- Initialize_parameters(X,R)
  W <- params$W0
  P_T <- params$P0_T
  U <- params$U
  rho <- params$rho
  alpha <- params$alpha
  
  # Initialize matrices and lists
  T_scores <- matrix(nrow = I, ncol = R)
  Lossc <- 1
  Lossvec <- Lossc
  est_sim_v <- c()
  true_sim_v <- c()
  # Loop
  while (convAO == 0) {
    
    # Update component scores
    T_scores <- X%*%W

    # Update loadings
    P_T = compute_P_new_T(X,W,T_scores,U,rho)

    # Compute b
    b <- compute_b(X,W,P_T, alpha)

    # Update weights
    W <- compute_w_new(X, R, P_T, b, alpha, rho, U, phi)
    
    # Update scaled variable
    U <- compute_U(U, W, P_T, rho)
    
    # Calculate loss
    Lossu <- loss_function(X,W,P_T,rho,U)
    Lossvec <- c(Lossvec,Lossu)
    
    #Check for convergence or if maximum iterations are reached
    if (iter > MaxIter) {
      convAO <- 1
      cat("Maxiter")
    }
    if (abs(Lossc-Lossu) < epsilon) {
      convAO <- 1
      cat("convergence")
    }
    print(paste("Iteration completed:", iter))
    iter <- iter + 1
    Lossc <- Lossu
  }
  
  print(Lossvec)
  uslpca <- list('weights' = W, 'loadings' = P_T, 'Lossvec' = Lossvec, 'Residual' = Lossu, 'Scores'= T_scores)
  return(uslpca)
}

compute_P_new_T <- function(X, W, T_scores, U, rho) {
    # Calculate X^T XW
    XtXW <- t(X) %*% T_scores
    
    # Add regularization term rho * (W + U)
    regularization_term <- rho * (W + U)
    
    # Combine the terms
    term1 <- 2 * XtXW + regularization_term
    
    # Calculate (2 * W^T X^T X W + rho * I)
    I <- diag(ncol(W))  # Identity matrix with size equal to number of columns of W
    term2 <- 2 *((t(W) %*% t(X) %*% T_scores) + (rho * I))
    
    # Inverse of term2
    term2_inv <- ginv(term2)
    
    # Multiply term1 by the inverse of term2
    result <- term1 %*% term2_inv
    
    # Transpose the result to get P_new^T
    P_new_T <- t(result)
    
    return(P_new_T)
}  

compute_b <- function(X,W_old,P_T, alpha){
  # Vectorized form of W and X
  vec_W = as.vector(W_old)
  vec_X = as.vector(X)
  P = t(P_T)
  
  # P kronecker X
  PX_kron = kronecker(P, X)

  # Compute: PX_kron^T*PX_kron*vec(W)
  term1 = t(PX_kron) %*% PX_kron %*% vec_W
  
  # PX_kron^T *vec(X)
  term2 = t(PX_kron) %*% vec_X
  
  # Subtract term2 from term 1 and dividing by alpha
  term3 = term1 - term2
  term4 = term3/alpha
  
  # Subtract vec_W - term 4
  b = vec_W - term4
  
  return(b)
  
}

compute_w_new <- function(X, R, P_T, b, alpha, rho, U, phi_prop) {
  vec_P <- as.vector(t(P_T))  # Flatten P_T row-wise
  W_new_vec <- (2 * alpha * as.vector(b) + rho * (vec_P - as.vector(U))) / (2 * alpha + rho)
  J <- dim(X)[2]
  W_new_matrix <- matrix(W_new_vec, nrow = J, ncol = R)
  for (r in 1:R) {
    b_col <- matrix(b, nrow = J, ncol = R)[, r]
    U_col <- matrix(U, nrow = J, ncol = R)[, r]
    P_col <- t(P_T)[, r]
    importance_scores <- (b_col)^2 + (U_col - P_col)^2
    sorted_indices <- order(importance_scores, decreasing = TRUE)
    mask <- rep(0, J)
    mask[sorted_indices[1:phi_prop]] <- 1
    W_new_matrix[, r] <- W_new_matrix[, r] * mask
  }
  
  return(W_new_matrix)
}

compute_U <- function(U,W,P_T,rho){
  U_new <- U + rho*(W- t(P_T))
  return(U_new)
}

loss_function <-function(X,W,P_T,rho,U){
  # First term
  recon_error <- sum((X - X %*% W %*% P_T)^2)
  
  # # Second term
  # equality_penalty <- (rho / 2) * sum((W - t(P_T) + U)^2)
  # 
  # # Quardratic term 
  # quadratic <- (rho / 2) * sum((U)^2)
  
  # Adding all loss together
  total_loss <- recon_error
  return(total_loss)
}

num_correct <- function (TargetW, EstimatedW){
  total_vnumber <- dim(TargetW)[1] * dim(TargetW)[2]
  TargetW[which(TargetW != 0)] <- 1
  sum_select <- sum(TargetW)
  sum_zero <- total_vnumber - sum_select
  EstimatedW[which(EstimatedW != 0)] <- 1
  total_correct <- sum(TargetW == EstimatedW) # this is the total number of variables correctedly selected and zeros correctly retained
  prop_correct <- total_correct/total_vnumber
  return(prop_correct)
}

similarity_estimated <- function(EstimatedW, EstimatedP_T){
  # How similar the estimated weights is to the estimated factor loadings
  P <- t(EstimatedP_T)
  difference <- sum(abs(EstimatedW-P))
  return(difference)
}

similarity_true <- function(EstimatedW, TrueP){
  # How similar is the weights to the true factor loadings
  difference <- sum(abs(EstimatedW-TrueP))
  return(difference)
}

similarity_scores <-function(EstimatedT, TrueT){
  difference <- sum(abs(EstimatedT-TrueT))
  return(difference)
}

explained_variance <- function(X, T_scores) {
  total_var <- sum(apply(X, 2, var))
  comp_vars <- apply(T_scores, 2, var)
  return(comp_vars / total_var)
}

reconstruction_metrics <- function(X, W, P_T) {
  X_hat <- X %*% W %*% P_T
  error_matrix <- X - X_hat
  mse <- mean(error_matrix^2)
  var_explained <- 1 - (sum(error_matrix^2) / sum((X - mean(X))^2))
  return(list(mse = mse, R2 = var_explained))
}

score_metrics <- function(T_est, T_true) {
  mae <- mean(abs(T_est - T_true))
  rmse <- sqrt(mean((T_est - T_true)^2))
  corrs <- diag(cor(T_est, T_true))  # assumes same column order
  avg_corr <- mean(corrs)
  return(list(mae = mae, rmse = rmse, correlation = avg_corr))
}

evaluate_variable_selection <- function(W_true, W_estimated) {
  W_true_bin <- ifelse(W_true != 0, 1, 0)
  W_est_bin <- ifelse(W_estimated != 0, 1, 0)
  
  TP <- sum(W_true_bin == 1 & W_est_bin == 1)
  FP <- sum(W_true_bin == 0 & W_est_bin == 1)
  FN <- sum(W_true_bin == 1 & W_est_bin == 0)
  TN <- sum(W_true_bin == 0 & W_est_bin == 0)
  
  precision <- TP / (TP + FP + 1e-8)
  recall <- TP / (TP + FN + 1e-8)
  f1_score <- 2 * (precision * recall) / (precision + recall + 1e-8)
  accuracy <- (TP + TN) / (TP + FP + FN + TN)
  
  return(list(
    precision = precision,
    recall = recall,
    f1 = f1_score,
    accuracy = accuracy
  ))
}

compute_bias_variance_mse <- function(W_true, W_est) {
  # Flatten matrices to vectors
  W_true_vec <- as.vector(W_true)
  W_est_vec <- as.vector(W_est)
  
  # Bias: mean difference
  bias <- mean(W_est_vec - W_true_vec)
  
  # Variance: variance of estimation errors
  variance <- var(W_est_vec - W_true_vec)
  
  # MSE: mean squared error
  mse <- mean((W_est_vec - W_true_vec)^2)
  
  return(list(bias = bias, variance = variance, mse = mse))
}

compute_path_metrics <- function(Beta_true, Beta_est) {
  bias_abs <- mean(abs(Beta_est - Beta_true))
  mse <- mean((Beta_est - Beta_true)^2)
  return(list(absolute_bias = bias_abs, mse = mse))
}
compute_variance <- function(Beta_est, Beta_true) {
  variance <- mean((Beta_est - Beta_true)^2) - (mean(Beta_est - Beta_true))^2
  return(variance)
}

explained_variance <- function(X, T_scores) {
  total_var <- sum(apply(X, 2, var))
  comp_vars <- apply(T_scores, 2, var)
  return(comp_vars / total_var)
}

sparsity_level <- function(W) {
  total_elements <- length(W)
  zero_elements <- sum(W == 0)
  return(zero_elements / total_elements)
}

compute_AVE <- function(X, T_scores) {
  R <- ncol(T_scores)
  AVEs <- numeric(R)
  
  for (r in 1:R) {
    correlations <- cor(X, T_scores[, r])  # correlations of all X variables with latent score r
    loadings_squared <- correlations^2
    AVEs[r] <- mean(loadings_squared)
  }
  
  return(AVEs)
}

compute_CR <- function(X, T_scores) {
  R <- ncol(T_scores)
  CRs <- numeric(R)
  
  for (r in 1:R) {
    correlations <- cor(X, T_scores[, r])
    loadings_squared <- correlations^2
    sum_loadings <- sum(correlations)
    numerator <- sum_loadings^2
    denominator <- numerator + sum(1 - loadings_squared)
    CRs[r] <- numerator / denominator
  }
  
  return(CRs)
}

compute_fornell_larcker_matrix <- function(X, T_scores) {
  R <- ncol(T_scores)
  AVEs <- compute_AVE(X, T_scores)
  inter_corrs <- cor(T_scores)
  fl_matrix <- matrix(NA, nrow = R, ncol = R)
  
  for (i in 1:R) {
    for (j in 1:R) {
      if (i == j) {
        fl_matrix[i, j] <- sqrt(AVEs[i])
      } else {
        fl_matrix[i, j] <- inter_corrs[i, j]
      }
    }
  }
  
  rownames(fl_matrix) <- paste0("Comp", 1:R)
  colnames(fl_matrix) <- paste0("Comp", 1:R)
  return(fl_matrix)
}


n_multistarts <- 3
base_seed <- 100

# load("Data-R-W-Sparse/Info_simulation.RData") # Contains design information
Info_matrix <- Infor_simulation$design_matrix_replication
Ndatasets <- Infor_simulation$n_data_sets
results_list <- list()

for (i in 8:8) {
  # Load data
  filename <- paste0("DATA-R-W-Sparse/Wsparse", i, ".RData")
  load(filename)
  
  # Extract data
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  T_true <- out$Z
  R <- dim(W_true)[2]
  conditions <- Info_matrix[i, ]
  phi <- as.numeric((1 - conditions[3]) * dim(X)[2])
  
  # MULTISTART IMPLEMENTATION
  best_result <- NULL
  best_loss <- Inf
  
  for (m in 1:n_multistarts) {
    set.seed(base_seed + m)  # optional for reproducibility
    result <- CEC_PLS_SEM(X, R, epsilon = 1e-6, phi)
    
    if (result$Residual < best_loss) {
      best_loss <- result$Residual
      best_result <- result
    }
  }
  
  # Metrics
  recovery_rate <- num_correct(best_result$weights, W_true)
  estimate_sim <- similarity_estimated(best_result$weights, best_result$loadings)
  true_sim <- similarity_true(best_result$weights, P_true)
  score_sim <- similarity_scores(best_result$Scores, T_true)
  bias_var_mse <- compute_bias_variance_mse(W_true, best_result$weights)
  expl_var <- explained_variance(X, best_result$Scores)
  sparsity <- sparsity_level(best_result$weights)
  AVE_scores <- compute_AVE(X, best_result$Scores)
  CR_scores <- compute_CR(X, best_result$Scores)
  FL_matrix <- compute_fornell_larcker_matrix(X, best_result$Scores)
  
  
  # Reconstruction metrics
  recon_metrics <- reconstruction_metrics(X, best_result$weights, best_result$loadings)
  score_metrics_list <- score_metrics(best_result$Scores, T_true)
  selection_eval <- evaluate_variable_selection(W_true, best_result$weights)
  
  # Store full results
  results_list[[i]] <- data.frame(
    Dataset = i,
    n_variables = conditions$n_variables,
    s_size = conditions$s_size,
    p_sparse = conditions$p_sparse,
    n_components = conditions$n_components,
    VAFx = conditions$VAFx,
    Recovery_Rate = recovery_rate,
    final_loss = best_result$Residual,
    est_sim = estimate_sim,
    true_sim = true_sim,
    score_sim = score_sim,
    mse = recon_metrics$mse,
    R2 = recon_metrics$R2,
    MAE = score_metrics_list$mae,
    RMSE = score_metrics_list$rmse,
    Corr = score_metrics_list$correlation,
    Precision = selection_eval$precision,
    Recall = selection_eval$recall,
    F1 = selection_eval$f1,
    Accuracy = selection_eval$accuracy,
    Bias = bias_var_mse$bias,
    Variance = bias_var_mse$variance,
    MSE_weights = bias_var_mse$mse,
    Expl_Var1 = expl_var[1],
    Expl_Var2 = ifelse(length(expl_var) >= 2, expl_var[2], NA),
    Sparsity = sparsity,
    AVE1 = AVE_scores[1],
    CR1 = CR_scores[1],
    AVE2 = ifelse(length(AVE_scores) >= 2, AVE_scores[2], NA),
    CR2 = ifelse(length(CR_scores) >= 2, CR_scores[2], NA)
  )
  
  cat("Dataset", i, "completed with best loss =", best_result$Residual, "\n")
}


results_table <- do.call(rbind, results_list)
write.csv(results_table, "CEC_PLS_SEM_updatedmetrics_results_1_18.csv", row.names = FALSE)
results_table

