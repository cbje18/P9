library("RTAQ")
library("dccmidas")
library("rumidas")
library(dplyr)
library(ggplot2)
library(latex2exp)

#data("realized_library")
#colnames(realized_library)

#data(sp500)
#plot(sp500)
#colnames(sp500)

data(rv5)
plot(rv5)

autoplot(rv5) + ylab("Realized Variance") + xlab(NULL)
autoplot(sqrt(rv5)) + ylab("Realized Volatility") + xlab(NULL) 

#Define the q's as in the article
q <- c(0.5, 1, 1.5, 2, 3)
Delta <- c(1:30)

data <- matrix(0, nrow=length(Delta), ncol = length(q))

#Compute m for every q and every Delta 
for(j in 1:length(q)){
  for(i in 1:length(Delta)){
    data(rv5)
    rv5
    
    log_vol <- log(sqrt(rv5$rv5)) %>% as.list()
    log_vol <- diff(log_vol$rv5, lag = Delta[i])
    log_vol <- log_vol[-(1:Delta[i])]
    
    m <- (1/length(log_vol))*sum(abs(log_vol)^q[j])
    m <- log(m)
    data[i,j] <- m
  }
}

data <- as.data.frame(data)
colnames(data) <- c("q = 0.5", "q = 1", "q = 1.5", "q = 2", "q = 3")


#Fit linear models to the computed data 
linear_models <- matrix(data = 0, nrow=2, ncol = length(q))
for(i in 1:length(q)){
  glm <- lm(data[,i]~log(Delta))
  linear_models[,i] <- glm$coefficients
}

linear_models

#Plot data points along with straight lines
ggplot(data=data) + geom_point(aes(x=log(Delta), y=data[,1],color="q=0.5")) + 
  geom_abline(aes(intercept=linear_models[1,1], slope = linear_models[2,1],color="q=0.5")) +
  geom_point(aes(x=log(Delta), y=data[,2],color = "q=1")) +
  geom_abline(aes(intercept=linear_models[1,2], slope = linear_models[2,2], color = "q=1")) + 
  geom_point(aes(x=log(Delta), y=data[,3],color = "q=1.5")) +
  geom_abline(aes(intercept=linear_models[1,3], slope = linear_models[2,3], color = "q=1.5")) +
  geom_point(aes(x=log(Delta), y=data[,4],color = "q=2")) +
  geom_abline(aes(intercept=linear_models[1,4], slope = linear_models[2,4],color = "q=2")) +
  geom_point(aes(x=log(Delta), y=data[,5],color = "q=3")) +
  geom_abline(aes(intercept=linear_models[1,5], slope = linear_models[2,5], color="q=3")) +
  xlab(TeX("\\log(\\Delta)")) + ylab(TeX("\\log(m(q,\\Delta))")) + 
  scale_color_manual(breaks =c("q=0.5","q=1","q=1.5","q=2","q=3"),
                     values = c("q=0.5"= "green", "q=1"="darkgreen", 
                                "q=1.5"="cornflowerblue", "q=2"="cyan", "q=3"="darkblue"))


#The estimated H values 
H <- linear_models[2,]/q; H
H_hat <- mean(H)


# Histograms of increments  -----------------------------------------------
#Histogram for Delta = 1
log_vol <- NULL
delta <- 1

data(rv5)
rv5

log_vol <- log(sqrt(rv5$rv5)) %>% as.list()
log_vol <- diff(log_vol$rv5, lag = delta)
log_vol <- log_vol[-(1:delta)]
log_vol <- log_vol %>% as.data.frame()

hist(log_vol$rv5, breaks = 40, prob = TRUE, main = TeX("\\Delta = 1"), xlab = "", ylab = "")
m <- mean(log_vol$rv5)
s <- sd(log_vol$rv5)
curve(dnorm(x, mean = m, sd = s), add = TRUE, col = "red", lwd = 1.5)

#Histogram for Delta = 5
log_vol <- NULL
delta <- 5

data(rv5)
rv5

log_vol <- log(sqrt(rv5$rv5)) %>% as.list()
log_vol <- diff(log_vol$rv5, lag = delta)
log_vol <- log_vol[-(1:delta)]
log_vol <- log_vol %>% as.data.frame()

hist(log_vol$rv5, breaks = 40, prob = TRUE, main = TeX("\\Delta = 5"), xlab = "", ylab = "")
m <- mean(log_vol$rv5)
s <- sd(log_vol$rv5)
curve(dnorm(x, mean = m, sd = s), add = TRUE, col = "red", lwd = 1.5)

#Histogram for Delta = 25
log_vol <- NULL
delta <- 25

data(rv5)
rv5

log_vol <- log(sqrt(rv5$rv5)) %>% as.list()
log_vol <- diff(log_vol$rv5, lag = delta)
log_vol <- log_vol[-(1:delta)]
log_vol <- log_vol %>% as.data.frame()

hist(log_vol$rv5, breaks = 40, prob = TRUE, main = TeX("\\Delta = 25"), xlab = "", ylab = "")
m <- mean(log_vol$rv5)
s <- sd(log_vol$rv5)
curve(dnorm(x, mean = m, sd = s), add = TRUE, col = "red", lwd = 1.5)

#Histogram for Delta = 125
log_vol <- NULL
delta <- 125

data(rv5)
rv5

log_vol <- log(sqrt(rv5$rv5)) %>% as.list()
log_vol <- diff(log_vol$rv5, lag = delta)
log_vol <- log_vol[-(1:delta)]
log_vol <- log_vol %>% as.data.frame()

hist(log_vol$rv5, breaks = 40, prob = TRUE, main = TeX("\\Delta = 125"), xlab = "", ylab = "")
m <- mean(log_vol$rv5)
s <- sd(log_vol$rv5)
curve(dnorm(x, mean = m, sd = s), add = TRUE, col = "red", lwd = 1.5)


# Parameter Estimation for Fractional Ornstein-Uhlenbeck ------------------
rm(list = ls())

data(rv5)

log_process <- log(sqrt(rv5$rv5)) %>% as.data.frame()
plot(log_process$rv5, type = "l")

#Estimating H
N <- length(rv5)
summands_num <- rep(0,N-4)
for(i in 1:(N-4)){
  summands_num[i] <- abs(log_process$rv5[i+4] - 2*log_process$rv5[i+2] + log_process$rv5[i])^2
}

summands_denom <- rep(0,N-2)
for(j in 1:(N-2)){
  summands_denom[j] <- abs(log_process$rv5[j+2] - 2*log_process$rv5[j+1] + log_process$rv5[j])^2
}

H_estimate <- 0.5*log2(sum(summands_num)/sum(summands_denom))
H_estimate

#Estimating phi
#phi_estimate <- sqrt(sum(summands_denom)/(N*(4-2^(2*H_estimate))*1^(2*H_estimate)))

tau <- 4*1^(2*H_estimate) -(2*1)^(2*H_estimate)
phi_estimate2 <- sqrt((1/5079*tau)*sum(summands_denom))


#Estimating theta 
theta_estimate <- (1/5079)*sum(log_process$rv5)

#diff_log <- diff(log_process$rv5, lag = 1)
#top <- (log_process$rv5[5079]-log_process$rv5[1])*sum(log_process$rv5^2)-sum(head(log_process$rv5,-1)*diff_log)*sum(log_process$rv5)

#bottom <- (log_process$rv5[5079]-log_process$rv5[1])*sum(log_process$rv5)-N*sum(log_process$rv5^2)-sum(head(log_process$rv5,-1)*diff_log)

#theta_estimate <- top/bottom

#Estimating lambda
lambda_num <- 5079*sum((log_process$rv5)^2) - (sum(log_process$rv5))^2
lambda_denom <- (5079^2)*(phi_estimate2^2)*H_estimate*gamma(2*H_estimate)

lambda_estimate <- (lambda_num/lambda_denom)^(-1/(2*H_estimate))

#top_lambda <- (N*sum(log_process$rv5^2)- (sum(log_process$rv5))^2)
#bottom_lambda <- (N^2)*(phi_estimate^2)*(H_estimate^2)*gamma(2*H_estimate)

#lambda_estimate <- (top_lambda/bottom_lambda)^(-1/(2*H_estimate))

#Simulating this process
Y <- sim_fOU_euler(T = 1, N = 5079, lambda = lambda_estimate, H= H_estimate, 
              phi=phi_estimate2, theta = theta_estimate, Y0=log_process$rv5[1])


Y <- head(Y, -1)
plot(Y, type = "l")
plot(exp(Y), type = "l")
sim_data <- cbind(t_grid, exp(Y)) %>% as.data.frame()

t_grid <- seq(1,5079,1)
hej <- sqrt(rv5$rv5)
data_for_plot <- cbind(t_grid, hej) %>% as.data.frame()
ggplot(data = data_for_plot) + geom_line(aes(x=t_grid, y = hej)) + xlab(NULL) + ylab(NULL) +
  ylim(0,0.061) + ggtitle("Data")


ggplot(data = sim_data) + geom_line(aes(x=t_grid, y = exp(Y)))+ xlab(NULL) + ylab(NULL) +
  ylim(0,0.061) + ggtitle("Model")


