library("RTAQ")
library("dccmidas")
library("rumidas")
library(dplyr)

data("realized_library")
colnames(realized_library)

data(sp500)
plot(sp500)
colnames(sp500)


data(rv5)
rv5
plot(rv5)

library(ggplot2)
autoplot(sqrt(rv5)) + ylab("Realized Volatility") + xlab(NULL) 

realized <- as.data.frame(rv5)
realized$rv5
plot(realized$rv5,type = "l")

#Define the q's as in the article
q <- c(0.5, 1, 1.5, 2, 3)
Delta <- c(1:30)

data <- matrix(0, nrow=length(Delta), ncol = length(q))
#Our Delta = 1 for one day, and N = 5079 observations
for(j in 1:length(q)){
  for(i in 1:length(Delta)){
    data(rv5)
    rv5
    
    log_vol <- log(rv5$rv5) %>% as.list()
    log_vol <- diff(log_vol$rv5, lag = Delta[i])
    log_vol <- log_vol[-(1:Delta[i])]
    
    m <- (1/length(log_vol))*sum(abs(log_vol[-1])^q[j])
    m <- log(m)
    data[i,j] <- m
  }
}

data <- as.data.frame(data)
colnames(data) <- c("q = 0.5", "q = 1", "q = 1.5", "q = 2", "q = 3")


linear_models <- matrix(data = 0, nrow=2, ncol = 5)
for(i in 1:5){
  glm <- lm(data[,i]~log(Delta))
  linear_models[,i] <- glm$coefficients
}

library(latex2exp)

ggplot(data=data) + geom_point(aes(x=log(Delta), y=data[,1]), color="q=0.5") + 
  geom_abline(intercept=linear_models[1,1], slope = linear_models[2,1], color="q=0.5") +
  geom_point(aes(x=log(Delta), y=data[,2]), color = "q=1") +
  geom_abline(intercept=linear_models[1,2], slope = linear_models[2,2], color ="q=1") + 
  geom_point(aes(x=log(Delta), y=data[,3]), color = "q=1.5") +
  geom_abline(intercept=linear_models[1,3], slope = linear_models[2,3], color = "q=1.5") +
  geom_point(aes(x=log(Delta), y=data[,4]), color = "q=2") +
  geom_abline(intercept=linear_models[1,4], slope = linear_models[2,4], color = "q=2") +
  geom_point(aes(x=log(Delta), y=data[,5]), color = "q=3") +
  geom_abline(intercept=linear_models[1,5], slope = linear_models[2,5], color = "q=3") +
  xlab(TeX("\\log(\\Delta)")) + ylab(TeX("\\log(m(q,\\Delta))")) + 
  scale_color_manual(breaks =c("q=0.5","q=1","q=1.5","q=2","q=3"),
                     values = c("q=0.5"= "green", "q=1"="red", 
                                "q=1.5"="blue", "q=2"="cyan", "q=3"="darkblue"))



linear_models


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
                     values = c("q=0.5"= "green", "q=1"="red", 
                                "q=1.5"="yellow", "q=2"="cyan", "q=3"="darkblue"))




