rm(list = ls())

library(microbenchmark)
library(tictoc)
library(ggplot2)

setwd("C:\\Users\\cbjkr\\OneDrive - Aalborg Universitet\\Aalborg Universitet\\Kandidat\\P9-Projekt\\Kode")
source("FunctionsForP9.R")

#testing the two methods against each other
sim_test <- microbenchmark(
  simfBM_cholesky(m = 1, T = 1, H = 0.3, n = 100),
  simfBM_fft(m = 1, T = 1, H = 0.3, N = 100),
  times = 100
)
sim_test

autoplot(sim_test)

boxplot(sim_test)

#testing each method for various sample sizes N
chol_test <- microbenchmark(
  simfBM_cholesky(m = 1, T = 1, H = 0.3, n = 10),
  simfBM_cholesky(m = 1, T = 1, H = 0.3, n = 100),
  simfBM_cholesky(m = 1, T = 1, H = 0.3, n = 1000),
  simfBM_cholesky(m = 1, T = 1, H = 0.3, n = 10000),
  times = 100
)
