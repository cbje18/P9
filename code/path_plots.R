rm(list = ls())

library(ggplot2)
library(latex2exp)

setwd("C:\\Users\\cbjkr\\OneDrive - Aalborg Universitet\\Aalborg Universitet\\Kandidat\\P9-Projekt\\Kode")
source("FunctionsForP9.R")


t_grid <- seq(0,1,1/10000)
X <- simfBM_fft(m = 1, T = 1, H = 0.8, N = 10000)
Y <- simfBM_fft(m = 1, T = 1, H = 0.8, N = 10000)
Z <- simfBM_fft(m = 1, T = 1, H = 0.8, N = 10000)

fBMs <- data.frame(t_grid,X,Y,Z)

ggplot(data=fBMs) + geom_line(aes(x=t_grid, y=X), colour = "blue") +
  geom_line(aes(x=t_grid, y=Y), colour = "#E69F00") + 
  geom_line(aes(x=t_grid, y=Z)) 

ggplot(data = fBMs) + geom_line(aes(x=t_grid, y = X), color ="#56B4E9") +
  xlab(TeX("$t$")) + ylab(TeX("$B_{t}^{H}$")) +
  theme(axis.text = element_text(size=12))

       