set.seed(123)

# Set working directory
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Load libraries
library(RGCCA)
library(qgraph)
library(MASS)  # for ginv, used by your functions if needed


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

set.seed(123)

# Set working directory
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Load libraries
library(RGCCA)
library(qgraph)
library(MASS)  # for ginv, used by your functions if needed


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
  mae <- mean(abs(est - true))
  rmse <- sqrt(mean((est - true)^2))
  return(list(mae = mae, rmse = rmse))
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

# Load simulation info
# load("DATA-R-P-Sparse/Info_simulation.RData")  # Adjust if path differs
Info_matrix <- Infor_simulation$design_matrix_replication
Ndatasets <- Infor_simulation$n_data_sets

results_list <- list()

for (i in 1:Ndatasets) {
  print(i)
  cat("Processing Dataset", i, "\n")
  
  filename <- paste0("DATA-R-P-Sparse/Psparse", i, ".RData")
  load(filename)  # loads `out` list
  
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  T_true <- out$Z  # True scores, if available
  
  R <- 2
  I <- nrow(X)
  J <- ncol(X)
  conditions <- Info_matrix[i, ]
  blocks <- list(bl1 = X)
  
  n_multistarts <- 1
  best_loss <- Inf
  best_fit_sgcca <- NULL
  best_estimatedW <- NULL
  
  n_vars <- ncol(blocks[[1]])
  var_names <- colnames(blocks[[1]])
  
  # Fit model only once (no multistarts)
  if (as.numeric(conditions[3]) == 0) {
    # Regular RGCCA (no sparsity)
    fit_sgcca <- rgcca(
      blocks = blocks,
      method = "rgcca",
      superblock = FALSE,
      ncomp = R,
      scheme = "factorial",
      comp_orth = FALSE,
      verbose = FALSE
    )
  } else {
    # Sparse SGCCA
    fit_sgcca <- rgcca(
      blocks = blocks,
      method = "sgcca",
      sparsity = as.numeric(conditions[3]),
      superblock = FALSE,
      ncomp = R,
      scheme = "factorial",
      comp_orth = FALSE,
      verbose = FALSE
    )
  }
  
  # Extract weights
  estimatedW <- matrix(0, ncol = R, nrow = n_vars)
  block_weights <- fit_sgcca$a$bl1
  
  if (is.null(rownames(block_weights))) {
    if (is.vector(block_weights)) {
      estimatedW[, 1:length(block_weights)] <- block_weights
    } else if (is.matrix(block_weights)) {
      estimatedW[,] <- block_weights
    }
  } else {
    # Extract row indices from block_weights row names like "V1_10"
    row_indices <- as.integer(sub("V1_", "", rownames(block_weights)))
    
    # Number of columns in block_weights
    n_col_bw <- ncol(block_weights)
    
    # Assign block_weights into estimatedW by row index and column index
    for (l in seq_along(row_indices)) {
      row_idx <- row_indices[l]
      estimatedW[row_idx, 1:n_col_bw] <- block_weights[l, ]
    }
  }
  
  
  # Scores and loadings
  T_scores <- fit_sgcca$Y$bl1
  P_estimated_reg <- solve(t(T_scores) %*% T_scores) %*% t(T_scores) %*% X  # R x p
  P_estimated <- t(P_estimated_reg)  # p x R to match loadings
  
  
  # Reconstruction loss
  X_hat <- T_scores %*% t(P_estimated)
  recon_mse <- mean((X - X_hat)^2)
  
  # Store results
  best_fit_sgcca <- fit_sgcca
  best_estimatedW <- estimatedW
  best_T_scores <- T_scores
  best_P_estimated <- P_estimated
  best_loss <- recon_mse
  
  # After best multistart selected, compute metrics:
  
  # Similarity metrics (weights vs loadings)
  sim_weights_loadings <- score_metrics(best_estimatedW, best_P_estimated)
  # Weights vs true
  sim_weights_true <- score_metrics(best_estimatedW, W_true)
  # Scores vs true
  sim_scores_true <- score_metrics(best_T_scores, T_true)
  # Loadings vs true
  sim_p_est_p_true <- score_metrics(best_P_estimated, P_true)
  
  # Reconstruction metrics
  recon_metrics <- reconstruction_metrics(X, best_estimatedW, best_P_estimated)
  
  # Variable selection evaluation
  selection_eval <- evaluate_variable_selection(W_true, best_estimatedW)
  
  # Bias-variance-MSE
  bias_var_mse <- compute_bias_variance_mse(W_true, best_estimatedW)
  
  # Explained variance per component
  expl_var <- explained_variance(X, best_T_scores)
  
  # Sparsity
  sparsity <- sparsity_level(best_estimatedW)
  
  # AVE and CR scores
  AVE_scores <- compute_AVE(X, best_T_scores)
  CR_scores <- compute_CR(X, best_T_scores)
  
  # Fornell-Larcker
  FL_values <- compute_fornell_larcker_values(X, best_T_scores)
  
  # Save results for this dataset
  results_list[[i]] <- data.frame(
    Dataset = i,
    Sample_Size = I,
    Num_Items = J,
    Num_Latent = R,
    Condition1 = conditions[1],
    Condition2 = conditions[2],
    Condition3 = conditions[3],
    Final_Loss = best_loss,
    
    # Similarity metrics
    P_vs_Ptrue_MAE = sim_p_est_p_true$mae,
    P_vs_Ptrue_RMSE = sim_p_est_p_true$rmse,
    
    W_vs_Loadings_MAE = sim_weights_loadings$mae,
    W_vs_Loadings_RMSE = sim_weights_loadings$rmse,
    
    W_vs_Wtrue_MAE = sim_weights_true$mae,
    W_vs_Wtrue_RMSE = sim_weights_true$rmse,
    
    Score_vs_True_MAE = sim_scores_true$mae,
    Score_vs_True_RMSE = sim_scores_true$rmse,
    
    # Reconstruction
    MSE_Recon = recon_metrics$mse,
    R2_Recon = recon_metrics$R2,
    
    # Selection metrics
    Precision = selection_eval$precision,
    Recall = selection_eval$recall,
    F1_Score = selection_eval$f1,
    Recovery = selection_eval$recovery,
    
    # Bias, variance, mse weights
    Bias = bias_var_mse$bias,
    Variance = bias_var_mse$variance,
    MSE_Weights = bias_var_mse$mse,
    
    # Explained variance components
    Explained_Var1 = expl_var[1],
    Explained_Var2 = ifelse(length(expl_var) >= 2, expl_var[2], NA),
    
    Sparsity = sparsity,
    
    AVE1 = AVE_scores$AVE[1],
    CR1 = CR_scores[1],
    AVE2 = ifelse(length(AVE_scores$AVE) >= 2, AVE_scores$AVE[2], NA),
    AVE_zero = AVE_scores$num_zero_variance_columns,
    CR2 = ifelse(length(CR_scores) >= 2, CR_scores[2], NA),
    
    Comp1_sqrtAVE = FL_values$Comp1_sqrtAVE,
    Comp2_sqrtAVE = FL_values$Comp2_sqrtAVE,
    Correlation_Comp1_Comp2 = FL_values$Correlation_Comp1_Comp2
  )
  cat("Dataset", i, "completed with best loss =", best_loss, "\n")
}

# Combine all results
results_table <- do.call(rbind, results_list)

# Save results
write.csv(results_table, "SGCCA_simulation_results_Psparse_correct_metrics.csv", row.names = FALSE)

# Print summary
print(results_table)

