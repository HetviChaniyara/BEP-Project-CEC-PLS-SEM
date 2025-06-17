# Cardinality and Equality Constrained PlS-SEM
# Hetvi Chaniyara
# Bachelor End Project

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# install.packages("MASS") 
library(MASS)
library(gtools)

# CEC-PLS-SEM Full Update Procedure
CEC_PLS_SEM <-function(X, R, epsilon, phi){
  
  J = dim(X)[2] # number of columns
  I = dim(X)[1] # number of rows
  iter <- 0
  convAO <- 0
  MaxIter <- 150 # Change accordingly
  
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
  
  # Update Loop
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
  results <- list('weights' = W, 'loadings' = P_T, 'Lossvec' = Lossvec, 'Residual' = Lossu, 'Scores'= T_scores, 'n_iterations'= iter)
  return(results)
}

########################################################################################################################################
# Helper Functions

Initialize_parameters <- function(X, R) {
  
  J <- dim(X)[2] # number of columns
  I <- dim(X)[1] # number of rows
  svd_X <- svd(X)
  
  # SVD-based components
  W_svd <- svd_X$v[, 1:R]
  P_svd <- t(svd_X$u[, 1:R])
  W0 <- W_svd
  P0_T <- P_svd
  
  # Random components
  # W_rand <- matrix(rnorm(length(W_svd), mean = 0, sd = 1), nrow = nrow(W_svd))
  # P_rand <- matrix(rnorm(length(P_svd), mean = 0, sd = 1), nrow = nrow(P_svd))
  
  # Weighted combination: 0.7 * SVD + 0.3 * random
  # W0 <- 0.8*W_svd + 0.2*W_rand
  # P0_T <- 0.8*P_svd + 0.2*P_rand
  
  U <- matrix(0, nrow = J, ncol = R) # Initialize to 0
  rho <- 1 # penalty parameter
  alpha <- max(eigen(t(X) %*% X)$values) # max eigenvalue of X^TX
  
  return(list(W0 = W0, P0_T = P0_T, U = U, rho = rho, alpha = alpha))
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
  W_new_matrix <- matrix(W_new_vec, nrow = J, ncol = R) # Reconstruct into a matrix
  
  # Coefficients with smallest bjr^2 + (Ujr-Pjr)^2 set to 0
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
  
  # Update U
  U_new <- U + rho*(W- t(P_T))
  
  return(U_new)
}

loss_function <-function(X,W,P_T,rho,U){
  
  # Loss function
  total_loss <- sum((X - X %*% W %*% P_T)^2)
  
  return(total_loss)
}

###############################################################################################################################
# Evaluation Metrics Functions

evaluate_variable_selection <- function(W_true, W_estimated) {
  
  # Checking which and how many coefficients are exactly 0
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
  
  return(list(precision = precision,recall = recall,f1 = f1_score,recovery = accuracy))
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
  
  # General Function For MAE, RMSE and Corrleation
  mae <- mean(abs(est - true))
  rmse <- sqrt(mean((est - true)^2))
  corrs <- diag(cor(est, true))  # assumes same column order
  avg_corr <- mean(corrs)
  
  return(list(mae = mae, rmse = rmse, correlation = avg_corr))
}

align_components <- function(est, true) {
  
  # Try combinations to see which estimated composite is corresponding one in the true matrix
  n_comp <- ncol(true)
  perm <- permutations(n_comp, n_comp)
  
  best_perm <- NULL
  best_score <- Inf
  
  # Selects permutation with best score and returns that order of composites
  for (i in 1:nrow(perm)) {
    aligned_est <- est[, perm[i, ]]
    
    # Flip signs for best match
    for (j in 1:n_comp) {
      if (cor(aligned_est[, j], true[, j]) < 0) {
        aligned_est[, j] <- -aligned_est[, j]
      }
    }
    
    score <- sum((aligned_est - true)^2)  # total squared error
    if (score < best_score) {
      best_score <- score
      best_perm <- aligned_est
    }
  }
  
  return(best_perm)
}


compute_bias_variance_mse <- function(W_true, W_est) {
  
  # Compute bias, variance and MSE
  W_true_vec <- as.vector(W_true)
  W_est_vec <- as.vector(W_est)
  bias <- mean(W_est_vec - W_true_vec)
  variance <- var(W_est_vec - W_true_vec)
  mse <- mean((W_est_vec - W_true_vec)^2)
  
  return(list(bias = bias, variance = variance, mse = mse))
}

explained_variance <- function(X, T_scores) {
  
  # Calculates the explained variance of the factor scores in ratio to the total variation in X
  total_var <- sum(apply(X, 2, var))
  comp_vars <- apply(T_scores, 2, var)
  return(comp_vars / total_var)
}

sparsity_level <- function(W) {
  
  # Checks the sparsity of the parameter
  total_elements <- length(W)
  zero_elements <- sum(W == 0)
  return(zero_elements / total_elements)
}

compute_AVE <- function(X, T_scores) {
  
  # Computes AVE
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
  
  # Computes Composite Reliability
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
  
  # Calculates Fornell-larcker criteria and stores it as three variables instead of a matrix
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

