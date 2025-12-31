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

# set random seed
set.seed(135)
# Before running the code, please create mixed_distributions_with_grain.csv using D:\MVCEMA\nmf\create_two_source_data.m 
X <- read.csv("D://MVCEMA//nmf//createddata//twosources//mixed_distributions_with_grain.csv", header = FALSE)

# Convert to matrix (only works if all columns are numeric)
X <- as.matrix(X)
dim(X)
X_rowname <- X[,1]
X <- X[,2:dim(X)[2]]
dim(X)
rownames(X) <- X_rowname
X <- t(X)
# Convert each column to numeric
numeric_matrix <- apply(X, 2, as.numeric)
X <- numeric_matrix

colnames(X)
rownames(X)
class(X)
dim(X)

apply(X,1,sum)

X <- diag(1/apply(X, 1, sum)) %*% X
apply(X,1,sum)

emma <- RECAauto(X, q = 2, c1 = -6, i = 1000, c5 = 0.5)

emma$B_mod[1:2,] <- emma$B_mod[2:1,]
emma$m_mod[,1:2] <- emma$m_mod[,2:1]

writeMat(con = "D:/MVCEMA/nmf/createddata/twosources/EMMAresultsfromR/EMMAtwosourceendmem.mat", EMMAtwosourceendmem = t(emma$B_mod))
writeMat(con = "D:/MVCEMA/nmf/createddata/twosources/EMMAresultsfromR/EMMAtwosourceabundances.mat", EMMAtwosourceabundances = t(emma$m_mod))

