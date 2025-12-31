rm(list=ls())
library('crosstalk')
library('rgl')
library('lba')
library(Ternary)
library(xtable)
library(R.matlab)
library(Matrix)
library(readxl)

source("code//emma//emma//c_err.R")
source("code//emma//emma//norm1.R")
source("code//emma//emma//r2.R")
source("code//emma//emma//tcalc.R")
source("code//emma//modify emma//RECAauto.R")
set.seed(135)

ranges = c('00', '05', '10', '15', '20', '25');

for (i in 1:length(ranges)) {
  X <- read.csv(paste0("D://MVCEMA//nmf//createddata//threesources//data_set_", ranges[i], ".csv"), header = FALSE)
  X <- as.matrix(X)
  X_rowname <- X[,1]
  X <- X[,2:dim(X)[2]]
  rownames(X) <- X_rowname
  X <- t(X)
  numeric_matrix <- apply(X, 2, as.numeric)
  X <- numeric_matrix
  X <- diag(1/apply(X, 1, sum)) %*% X
  
  emma <- RECAauto(X, q = 3, c1 = -6, i = 1000, c5 = 0.5)
  
  writeMat(con = paste0("D://MVCEMA//nmf//createddata//threesources/EMMAresultsfromR/EMMA_mixed_", ranges[i], "_endmem.mat"), EMMA_mixed_endmem = t(emma$B_mod))
  writeMat(con = paste0("D://MVCEMA//nmf//createddata//threesources/EMMAresultsfromR/EMMA_mixed_", ranges[i], "_abundances.mat"), EMMA_mixed_abundances = t(emma$m_mod))
  }
