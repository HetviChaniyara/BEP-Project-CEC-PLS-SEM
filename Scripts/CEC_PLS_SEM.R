# Cardinality and Equality Constrained PlS-SEM
# Hetvi Chaniyara
# Bachelor End Project


current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(current_working_dir)
getwd()

# load("DATA-R/Wsparse1.RData") # Load one simulation
# # install.packages("MASS") # Install the package if you haven't already
library(MASS) # Load the package
# X <- out$X  # Observed data (what you'd use for PLS-SEM)
# W_true <- out$W  # Sparse weights (for study, if needed)
# P_true <- out$P  # Loadings
# T_true <- out$Z  # Component scores
# R <- dim(W_true)[2]# number of latent variables


Initialize_parameters<-function(X, R){
  J = dim(X)[2] # number of columns
  I = dim(X)[1] # number of rows
  W0 <- matrix(rnorm(J * R), nrow = J, ncol = R)
  P0_T <- matrix(rnorm(R * J), nrow = R, ncol = J)
  U <- matrix(0, nrow = J, ncol = R) # Initialize to 0
  rho <- 10 # penalty parameter
  alpha <- max(eigen(t(X) %*% X)$values) #max eigen value of X^TX
  return(list(W0 = W0, P0_T = P0_T, U = U, rho = rho, alpha = alpha))
}

CEC_PLS_SEM <-function(X, R, epsilon, phi){
  J = dim(X)[2] # number of columns
  I = dim(X)[1] # number of rows
  iter <- 0
  convAO <- 0
  MaxIter <- 100
  
  # Get initialized parameters
  params <- Initialize_parameters(X,R)
  W <- params$W0
  P_T <- params$P0_T
  U <- params$U
  rho <- params$rho
  alpha <- params$alpha
  # Initialize matrices and lists
  T <- matrix(nrow = I, ncol = R)
  Lossc <- 1
  Lossvec <- Lossc
  est_sim_v <- c()
  true_sim_v <- c()
  # Loop
  while (convAO == 0) {
    
    #1. Update component scores
    T <- X %*% W

    #2. Update loadings
    P_T = compute_P_new_T(X,W,U,rho)

    #3. Compute b
    b <- compute_b(X,W,P_T, alpha)

    #4. Update weights
    W <- compute_w_new(X, R, P_T, b, alpha, rho, U, phi)
    
    #5. Update scaled variable
    U <- compute_U(U, W, P_T, rho)
    
    #Calculate loss
    Lossu <- loss_function(X,W,P_T,rho,U)
    Lossvec <- c(Lossvec,Lossu)
    # estimate_sim <- similarity_estimated(result$weights, result$loadings)
    # cat("Estimate_sim", estimate_sim)
    # true_sim <- similarity_true(result$weights, P_true)
    # cat("true_sim", true_sim)
    # est_sim_v <- c(est_sim_v, estimate_sim)
    # true_sim_v <- c(true_sim_v, true_sim)
    #Check for convergence or if maximum iterations are reached
    if (iter > MaxIter) {
      convAO <- 1
      cat("Maxiter")
    }
    if (abs(Lossc-Lossu) < epsilon) {
      convAO <- 1
      cat("convergence")
    }
    iter <- iter + 1
    Lossc <- Lossu
    print(paste("Iteration completed:", iter))
  }
  uslpca <- list('weights' = W, 'loadings' = P_T, 'Lossvec' = Lossvec, 'Residual' = Lossu)
  return(uslpca)
}

compute_P_new_T <- function(X, W, U, rho) {
    #Calculate X^T XW
    XtXW <- t(X) %*% X %*% W
    
    #Add regularization term rho * (W + U)
    regularization_term <- rho * (W + U)
    
    # Combine the terms
    term1 <- 2 * XtXW + regularization_term
    
    # Calculate (2 * W^T X^T X W + rho * I)
    I <- diag(ncol(W))  # Identity matrix with size equal to number of columns of W
    term2 <- 2 *((t(W) %*% t(X) %*% X %*% W) + (rho * I))
    
    # inverse of term2
    term2_inv <- ginv(term2)
    
    #Multiply term1 by the inverse of term2
    result <- term1 %*% term2_inv
    
    #Transpose the result to get P_new^T
    P_new_T <- t(result)
    
    return(P_new_T)
}  

compute_b <- function(X,W_old,P_T, alpha){
  # vectorized form of W and X
  vec_W= as.vector(W_old)
  vec_X = as.vector(X)
  P = t(P_T)
  # P kronecker X
  PX_kron = kronecker(P, X)

  # compute: PX_kron^T*PX_kron*vec(W)
  term1 = t(PX_kron) %*% PX_kron %*% vec_W
  
  # PX_kron^T *vec(X)
  term2= t(PX_kron) %*% vec_X
  
  # Subtract term2 from term 1 and dividing by alpha
  term3 = term1 - term2
  term4 = term3/alpha
  
  # subtract vec_W - term 4
  b = vec_W - term4
  
  return(b)
  
}

compute_w_new <- function(X, R, P_T, b, alpha, rho, U, phi_prop) {
  vec_P = as.vector(t(P_T))
  W_new_vec <- (2 * alpha * b + rho * (vec_P - as.vector(U))) / (2 * alpha + rho)
  
  J = dim(X)[2]
  W_new_matrix <- matrix(W_new_vec, nrow = J, ncol = R)
  
  # Apply column-wise sparsity
  for (r in 1:R) {
    col_vec <- W_new_matrix[, r]
    phi_r <- floor(phi_prop * J)  # Number of non-zeros to keep in this column
    sorted_indices <- order(abs(col_vec), decreasing = TRUE)
    col_vec[sorted_indices[(phi_r + 1):J]] <- 0
    W_new_matrix[, r] <- col_vec
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
  
  # Second term
  equality_penalty <- (rho / 2) * sum((W - t(P_T) + U)^2)

  # Quardratic term 
  quadratic <- (rho / 2) * sum((U)^2)
  
  # Adding all loss together
  total_loss <- recon_error + equality_penalty + quadratic
  
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
  P <- t(EstimatedP_T)
  difference <- sum(abs(EstimatedW-P))
  return(difference)
}

similarity_true <- function(EstimatedW, TrueP){
  difference <- sum(abs(EstimatedW-TrueP))
  return(difference)
}

# result <- CEC_PLS_SEM(X,R, epsilon = 10^-6)
# result
# 
# eval <- num_correct(result$weights, W_true)
# eval

# 
# # load("Data-R/Info_simulation.RData") # Contains design information
Info_matrix <- Infor_simulation$design_matrix_replication
Ndatasets <- Infor_simulation$n_data_sets
results_list <- list()

for (i in 1:2) {

  # Load each data file
  filename <- paste0("DATA-R/Wsparse", i, ".RData")
  load(filename)

  # Extract data
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  R <- dim(W_true)[2]  # number of latent variables
  conditions <- Info_matrix[i, ]
  phi <- as.numeric((1-conditions[3])*dim(X)[2]*dim(X)[1])
  print(phi)
  # Run your CEC-PLS-SEM model
  result <- CEC_PLS_SEM(X, R, epsilon = 10^-6, phi)

  # Evaluate recovery rate
  recovery_rate <- num_correct(result$weights, W_true)
  estimate_sim <- similarity_estimated(result$weights, result$loadings)
  true_sim <- similarity_true(result$weights, P_true)
  print(abs(sum(W_true-P_true)))
  # Get sample size and number of items
  I <- nrow(X)
  J <- ncol(X)

  # Get simulation conditions from Info_simulation
  print(W_true)
  print(result$weights)
  print(P_true)
  print(result$loadings)
  # Store everything in a small data.frame
  results_list[[i]] <- data.frame(
    Dataset = i,
    Sample_Size = I,
    Num_Items = J,
    Num_latent = R,
    Recovery_Rate = recovery_rate,
    Condition1 = conditions[1],
    Condition2 = conditions[2],
    Condition3 = conditions[3],
    final_loss = result$Residual,
    est_sim = estimate_sim,
    true_sim = true_sim
  )

  cat("Dataset", i, "completed\n")
}


results_table <- do.call(rbind, results_list)
# write.csv(results_table, "CEC_PLS_SEM_simulation_results_correct_test.csv", row.names = FALSE)
results_table

# m <- matrix(1:9, nrow = 3, ncol = 3)
# m
# P_T <- matrix(1:9, nrow=3, ncol=3)
# vec_m <- as.vector(m)
# vec_m
# W_sparse_vec <- vec_m
# abs_sorted_indices <- order(abs(vec_m), decreasing = TRUE)
# print(abs_sorted_indices)
# W_sparse_vec[abs_sorted_indices[(4+1):length(vec_m)]] <- 0
# print(W_sparse_vec)
# W_sparse_matrix <- matrix(W_sparse_vec, nrow = 3, ncol = 3)
# print(W_sparse_matrix)
