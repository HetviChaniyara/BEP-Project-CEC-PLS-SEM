set.seed(123)

# Set working directory
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Load libraries
library(RGCCA)
library(qgraph)
library(MASS) 

###########################################################################################
# Evaluation Metrics

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
  X_hat <- X %*% W %*% t(P_T)
  error_matrix <- X - X_hat
  mse <- mean(error_matrix^2)
  var_explained <- 1 - (sum(error_matrix^2) / sum((X - mean(X))^2))
  return(list(mse = mse, R2 = var_explained))
}

score_metrics <- function(est, true) {
  
  # Compute MAE and RMSE
  mae <- mean(abs(est - true))
  rmse <- sqrt(mean((est - true)^2))
  return(list(mae = mae, rmse = rmse))
}

compute_bias_variance_mse <- function(W_true, W_est) {
  
  # compute bias, variance and mse
  W_true_vec <- as.vector(W_true)
  W_est_vec <- as.vector(W_est)
  bias <- mean(W_est_vec - W_true_vec)
  variance <- var(W_est_vec - W_true_vec)
  mse <- mean((W_est_vec - W_true_vec)^2)
  
  return(list(bias = bias, variance = variance, mse = mse))
}

explained_variance <- function(X, T_scores) {
  
  # Compute variance explained by T scores in comparison to total variance in X
  total_var <- sum(apply(X, 2, var))
  comp_vars <- apply(T_scores, 2, var)
  return(comp_vars / total_var)
}

sparsity_level <- function(W) {
  
  # Compute sparsity in parameter
  total_elements <- length(W)
  zero_elements <- sum(W == 0)
  return(zero_elements / total_elements)
}

compute_AVE <- function(X, T_scores) {
  
  # Compute AVE
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
  
  # Compute Composite Reliability
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
  
  # Computes Fornell-Larcker
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


