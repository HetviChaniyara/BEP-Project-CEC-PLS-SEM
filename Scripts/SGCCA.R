set.seed(123)

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()


library(RGCCA)
library(qgraph)


num_correct <- function (TargetW, EstimatedW){
  total_vnumber <- dim(TargetW)[1] * dim(TargetW)[2]
  TargetW[TargetW != 0] <- 1
  EstimatedW[EstimatedW != 0] <- 1
  total_correct <- sum(TargetW == EstimatedW)
  prop_correct <- total_correct / total_vnumber
  return(prop_correct)
}

load("~/BEP-Project-CEC-PLS-SEM/Scripts/DATA-R/Info_simulation.RData")  # Contains Info_matrix
Info_matrix <- Infor_simulation$design_matrix_replication

results_list <- list()

for (i in 1:18) {
  cat("Processing Dataset", i, "\n")
  
  filename <- paste0("DATA-R/Wsparse", i, ".RData")
  load(filename)
  
  X <- out$X
  W_true <- out$W
  P_true <- out$P
  R <- ncol(W_true)
  I <- nrow(X)
  J <- ncol(X)
  conditions <- Info_matrix[i, ]

  blocks <- list(bl1 = X)
  
  n_multistarts <- 100
  best_recovery_rate <- -Inf  
  best_fit_sgcca <- NULL
  
  for (start in 1:n_multistarts) {
    cat("Multistart", start, "\n")

    fit_sgcca <- rgcca(blocks = blocks,
                       method = "sgcca",
                       sparsity = as.numeric(conditions[3])+0.00000001,  # You can vary this later
                       superblock = FALSE,
                       ncomp = R,
                       scheme = "factorial",
                       comp_orth = FALSE,
                       verbose = FALSE)

    EstimatedW <- fit_sgcca$a$bl1
    
    recovery_rate <- num_correct(W_true, EstimatedW)
    
    if (recovery_rate > best_recovery_rate) {
      best_recovery_rate <- recovery_rate
      best_fit_sgcca <- fit_sgcca
    }
  }
  
  EstimatedW_best <- best_fit_sgcca$a$bl1
  print(fit_sgcca$crit)
  results_list[[i]] <- data.frame(
    Dataset = i,
    Sample_Size = I,
    Num_Items = J,
    Num_Latent = R,
    Condition1 = conditions[1],
    Condition2 = conditions[2],
    Condition3 = conditions[3],
    Recovery_Rate = best_recovery_rate,
    Nonzero_C1 = sum(EstimatedW_best[, 1] != 0),
    Nonzero_C2 = ifelse(R >= 2, sum(EstimatedW_best[, 2] != 0), NA)
  )
}

results_table <- do.call(rbind, results_list)

write.csv(results_table, "SGCCA_simulation_results_multistart2.csv", row.names = FALSE)

print(results_table)
