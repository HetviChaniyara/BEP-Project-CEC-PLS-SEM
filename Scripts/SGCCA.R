# install.packages('qgraph')
# install.packages("RGCCA", dependencies = TRUE, INSTALL_opts = "--no-lock")

library(qgraph)
library(RGCCA)

load("DATA-R/Wsparse7.RData")
DATA <- out$X

blocks <- list(bl1 = DATA)

fit_sgcca <- rgcca(blocks = blocks, method = "sgcca", sparsity = 0.5, superblock = FALSE, ncomp = 2,
                   
                   scheme = "factorial",
                   
                   comp_orth = TRUE,
                   
                   verbose = TRUE)

nonzero1 <- sum(fit_sgcca$a$bl1[,1]!=0)

nonzero2 <- sum(fit_sgcca$a$bl1[,2]!=0)

# nonzero3 <- sum(fit_sgcca$a$bl1[,3]!=0)

# nonzero4 <- sum(fit_sgcca$a$bl1[,4]!=0)
# 
# nonzero5 <- sum(fit_sgcca$a$bl1[,5]!=0)

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

num_correct(out$W, fit_sgcca$a$bl1)
