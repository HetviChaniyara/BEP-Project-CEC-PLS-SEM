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
  # W_rand <- matrix(rnorm(length(W_svd), mean = 0, sd = 1), nrow = nrow(W_svd))
  # P_rand <- matrix(rnorm(length(P_svd), mean = 0, sd = 1), nrow = nrow(P_svd))

  # Weighted combination: 0.7 * SVD + 0.3 * random
  # W0 <- 0.8*W_svd + 0.2*W_rand
  # P0_T <- 0.8*P_svd + 0.2*P_rand
  W0 <- W_svd
  P0_T <- P_svd
  
  
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
  MaxIter <- 150
  
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
  uslpca <- list('weights' = W, 'loadings' = P_T, 'Lossvec' = Lossvec, 'Residual' = Lossu, 'Scores'= T_scores, 'n_iterations'= iter)
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
  total_loss <- sum((X - X %*% W %*% P_T)^2)
  return(total_loss)
}

## Evaluation Metrics

evaluate_variable_selection <- function(W_true, W_estimated) {
  # Variable Selection Metrics
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
    recovery = accuracy
  ))
}

reconstruction_metrics <- function(X, W, P_T) {
  # Reconstruction Metrics
  X_hat <- X %*% W %*% P_T
  error_matrix <- X - X_hat
  mse <- mean(error_matrix^2)
  var_explained <- 1 - (sum(error_matrix^2) / sum((X - mean(X))^2))
  return(list(mse = mse, R2 = var_explained))
}

score_metrics <- function(est, true) {
  # General Function
  mae <- mean(abs(est - true))
  rmse <- sqrt(mean((est - true)^2))
  corrs <- diag(cor(est, true))  # assumes same column order
  avg_corr <- mean(corrs)
  return(list(mae = mae, rmse = rmse, correlation = avg_corr))
}

compute_bias_variance_mse <- function(W_true, W_est) {
  W_true_vec <- as.vector(W_true)
  W_est_vec <- as.vector(W_est)
  bias <- mean(W_est_vec - W_true_vec)
  variance <- var(W_est_vec - W_true_vec)
  mse <- mean((W_est_vec - W_true_vec)^2)
  
  return(list(bias = bias, variance = variance, mse = mse))
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
  
  # Identify zero variance columns
  zero_var_cols <- apply(X, 2, sd) == 0
  num_zero_cols <- sum(zero_var_cols)
  
  # Remove zero variance columns from X
  X_filtered <- X[, !zero_var_cols, drop = FALSE]
  
  for (r in 1:R) {
    correlations <- cor(X_filtered, T_scores[, r])
    loadings_squared <- correlations^2
    AVEs[r] <- mean(loadings_squared)
  }
  
  return(list(
    AVE = AVEs,
    num_zero_variance_columns = num_zero_cols
  ))
}


compute_CR <- function(X, T_scores) {
  R <- ncol(T_scores)
  CRs <- numeric(R)
  
  # Identify columns with non-zero variance
  nonzero_var_cols <- apply(X, 2, sd) != 0
  X_filtered <- X[, nonzero_var_cols, drop = FALSE]
  
  for (r in 1:R) {
    correlations <- cor(X_filtered, T_scores[, r])
    loadings_squared <- correlations^2
    sum_loadings <- sum(correlations)
    numerator <- sum_loadings^2
    denominator <- numerator + sum(1 - loadings_squared)
    CRs[r] <- numerator / denominator
  }
  
  return(CRs)
}

compute_fornell_larcker_values <- function(X, T_scores) {
  AVE_result <- compute_AVE(X, T_scores)
  AVEs <- AVE_result$AVE
  inter_corr <- cor(T_scores)[1, 2]  # correlation between Comp1 and Comp2
  
  values <- list(
    Comp1_sqrtAVE = sqrt(AVEs[1]),
    Comp2_sqrtAVE = sqrt(AVEs[2]),
    Correlation_Comp1_Comp2 = inter_corr
  )
  
  return(values)
}

n_multistarts <- 3
base_seed <- 100

# load("Data-R-W-Sparse/Info_simulation.RData") # Contains design information
Infor_matrix <- Infor_simulation$design_matrix_replication
Ndatasets <- Infor_simulation$n_data_sets
results_list <- list()

for (i in 1:Ndatasets) {
  if (i %% 3 == 0) {
    cat("Skipping dataset", i, "because it is divisible by 3\n")
    next  # skip to the next iteration
  }
  
  filename <- paste0("DATA-R-P-Sparse/Psparse", i, ".RData")
  load(filename)
  print("hello")
  # Extract data
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  T_true <- out$Z
  R <- dim(W_true)[2]
  conditions <- Infor_matrix[i, ]
  phi <- as.numeric((1 - conditions[3]) * dim(X)[2])
  print("indentified")
  # MULTISTART IMPLEMENTATION
  best_result <- NULL
  best_loss <- Inf
  
  FL_matrix_list <- list()
  
  for (m in 1:n_multistarts) {
    set.seed(base_seed + m)  # optional for reproducibility
    print("started")
    result <- CEC_PLS_SEM(X, R, epsilon = 1e-6, phi)
    
    if (result$Residual < best_loss) {
      best_loss <- result$Residual
      best_result <- result
    }
  }
  
  # Unified similarity metrics using score_metrics
  sim_weights_loadings <- score_metrics(best_result$weights, t(best_result$loadings))
  sim_weights_true <- score_metrics(best_result$weights, W_true)
  sim_scores_true <- score_metrics(best_result$Scores, T_true)
  sim_p_est_p_true <- score_metrics(t(best_result$loadings), P_true)
  
  # Reconstruction, AVE/CR, and selection metrics
  bias_var_mse <- compute_bias_variance_mse(W_true, best_result$weights)
  expl_var <- explained_variance(X, best_result$Scores)
  sparsity <- sparsity_level(best_result$weights)
  AVE_scores <- compute_AVE(X, best_result$Scores)
  CR_scores <- compute_CR(X, best_result$Scores)
  FL_values <- compute_fornell_larcker_values(X, best_result$Scores)
  
  recon_metrics <- reconstruction_metrics(X, best_result$weights, best_result$loadings)
  selection_eval <- evaluate_variable_selection(W_true, best_result$weights)
  
  # Store all results in consistent format
  results_list[[i]] <- data.frame(
    Dataset = i,
    n_variables = conditions$n_variables,
    s_size = conditions$s_size,
    p_sparse = conditions$p_sparse,
    n_components = conditions$n_components,
    VAFx = conditions$VAFx,
    Final_Loss = best_result$Residual,
    Num_Iterations = best_result$n_iterations,
    
    
    # Similarity metrics)
    P_vs_Ptrue_MAE = sim_p_est_p_true$mae,
    P_vs_Ptrue_RMSE = sim_p_est_p_true$rmse,
    P_vs_Ptrue_Corr = sim_p_est_p_true$correlation,
    
    W_vs_Loadings_MAE = sim_weights_loadings$mae,
    W_vs_Loadings_RMSE = sim_weights_loadings$rmse,
    W_vs_Loadings_Corr = sim_weights_loadings$correlation,
    
    W_vs_Wtrue_MAE = sim_weights_true$mae,
    W_vs_Wtrue_RMSE = sim_weights_true$rmse,
    W_vs_Wtrue_Corr = sim_weights_true$correlation,
    
    Score_vs_True_MAE = sim_scores_true$mae,
    Score_vs_True_RMSE = sim_scores_true$rmse,
    Score_vs_True_Corr = sim_scores_true$correlation,
    
    # Reconstruction
    MSE_Recon = recon_metrics$mse,
    R2_Recon = recon_metrics$R2,
    
    # Selection
    Precision = selection_eval$precision,
    Recall = selection_eval$recall,
    F1_Score = selection_eval$f1,
    Recovery = selection_eval$recovery,
    
    # Bias-variance-MSE
    Bias = bias_var_mse$bias,
    Variance = bias_var_mse$variance,
    MSE_Weights = bias_var_mse$mse,
    
    # Other
    Explained_Var1 = expl_var[1],
    Explained_Var2 = ifelse(length(expl_var) >= 2, expl_var[2], NA),
    Sparsity = sparsity,
    AVE1 = AVE_scores$AVE[1],
    CR1 = CR_scores[1],
    AVE2 = ifelse(length(AVE_scores) >= 2, AVE_scores$AVE[2], NA),
    AVE_zero = AVE_scores$num_zero_variance_columns,
    CR2 = ifelse(length(CR_scores) >= 2, CR_scores[2], NA),
    Comp1_sqrtAVE = FL_values$Comp1_sqrtAVE,
    Comp2_sqrtAVE = FL_values$Comp2_sqrtAVE,
    Correlation_Comp1_Comp2 = FL_values$Correlation_Comp1_Comp2
  )
  cat("Dataset", i, "completed with best loss =", best_result$Residual, "\n")
}

results_table <- do.call(rbind, results_list)
write.csv(results_table, "CEC_PLS_SEM_P_Spase_updated_all_except_3s_svd_only.csv", row.names = FALSE)
results_table