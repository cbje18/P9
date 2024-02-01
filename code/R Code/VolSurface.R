rm(list = ls())
graphics.off()

library(quantmod)
library(dplyr)
library(plotly)
library(akima)
library(reshape2)

setwd("C:\\Users\\cbjkr\\OneDrive - Aalborg Universitet\\Aalborg Universitet\\Kandidat\\P9-Projekt\\P9\\Data")
data <- read.csv(file = "option_data_10november.csv")
#First get options data from Yahoo
ticker <- "SPY"
expiration_dates <- c("2023-12-06",
                      "2023-12-07", "2023-12-08", "2023-12-11",
                      "2023-12-12", "2023-12-13", "2023-12-14",
                      "2023-12-15", "2023-12-22","2023-12-29",
                      "2024-01-19", "2024-02-16", "2024-03-15", "2024-03-28",
                      "2024-06-21", "2024-06-28", "2024-09-20", "2024-09-30",
                      "2024-12-20", "2025-01-17", "2025-03-21")

options_data <- getOptionChain(ticker, Exp = expiration_dates)

options_data_matrix <- NULL
for(i in 1:15){
  options_data_matrix <- rbind(options_data_matrix, options_data[[i]]$calls)
}

#Set current stock price
S0 <- 456.04


options_data_matrix$moneyness <- S0/options_data_matrix$Strike
options_data_matrix$logmoneyness <- log(S0/options_data_matrix$Strike)


nov13 <- data[data$Expiration == "2023-11-13",]
ggplot(data = nov13, aes(x=Strike, y=IV)) + geom_point(size = 2) +
  labs(x = "Strike",  y = "Implied Volatility") +
  ggtitle("Call Options on SPY, Expiration Date 13/11-2023")


options_data_matrix$moneyness <- S0/options_data_matrix$Strike
options_data_matrix$log_moneyness <- log(S0/options_data_matrix$Strike)

ivGridCalls <- acast(options_data_matrix, TimeToExp ~ log_moneyness, value.var = "IV")

options_data_matrix$TimeToExp <- as.numeric(as.Date(options_data_matrix$Expiration) - as.Date(Sys.Date()))/365.25

setwd("C:\\Users\\cbjkr\\OneDrive - Aalborg Universitet\\Aalborg Universitet\\Kandidat\\P9-Projekt\\P9\\Data")
hej <- read.csv("option_data_10november.csv")
View(hej)

pek <- hej[hej$Expiration == "2023-11-10",]
plot(x = pek$Strike, y = pek$IV)
pek$log_moneyness <- log(437.25/pek$Strike)
plot(x = pek$log_moneyness, y = pek$IV)
#Create plot
axx <- list(
  gridcolor='rgb(255, 255, 255)',
  zerolinecolor='rgb(255, 255, 255)',
  showbackground=TRUE,
  backgroundcolor='rgb(230, 230,230)'
)

axx <- list(
  title = "Log-moneyness"
)

axy <- list(
  title = "Time to Exp."
)

axz <- list(
  title = "Implied vol."
)

fig <- plot_ly(x = ~options_data_matrix$log_moneyness, 
               y = ~options_data_matrix$TimeToExp,
               z = ~options_data_matrix$IV,
               type = 'mesh3d') 

fig

setwd("C:\\Users\\cbjkr\\OneDrive - Aalborg Universitet\\Aalborg Universitet\\Kandidat\\P9-Projekt\\P9")
options_data_matrix <- options_data_matrix %>% as.data.frame()
write.csv(options_data_matrix, "C:\\Users\\cbjkr\\OneDrive - Aalborg Universitet\\Aalborg Universitet\\Kandidat\\P9-Projekt\\options_spy_december.csv")
#Plot volatility smile for a given expiration date
options_data_matrix$Expiration <- as.Date(options_data_matrix$Expiration)

expiration_date <- options_data_matrix[options_data_matrix$Expiration == "2024-01-31",]

ggplot(data = expiration_date, aes(x=Strike, y=IV)) + geom_point() +
  labs(x = "Strike",  y = "Implied Volatility") 



# Heston Pricing ----------------------------------------------------------
library(NMOF)

data <- options_data_matrix
callHestoncf()
charFunc <- function(u, spot, ttm, v, sigma, kappa, theta, rho){
  d <- sqrt( (rho*sigma*u*1i - kappa)^2 + (sigma^2)*(u*1i + u^2) )
  d <- -d 
  g <- (kappa - rho*sigma*u*1i + d)/(kappa - rho*sigma*u*1i - d)
  
  tempM <- 
    (kappa-rho*sigma*u*1i+d)*ttm - 2*log((1-g*exp(d*ttm))/(1-g))
  M <- (kappa*theta)/(sigma^2)*tempM
  N <- (kappa-rho*sigma*u*1i + d)/(sigma^2)*((1-exp(d*ttm))/(1-g*exp(d*ttm)))
  res <- exp(M+N*v + 1i*u*log(spot))
  
  return(res)
}


HestonCall <- function(spot, strike, ttm, v, sigma, kappa, theta, rho){
  integrand1 <- function(u){
    num1 <- charFunc(u-1i, spot, ttm, v, sigma, kappa, theta, rho)
    den1 <- charFunc(-1i, spot, ttm, v, sigma, kappa, theta, rho)
    dummy1 <- exp(-1i*u*log(strike))*num1/(1i*u*den1)
    integrand1 <- Re(dummy1)
  }
  
  integrand2 <- function(u){
    dummy2 <- exp(-1i*u*log(strike))*charFunc(u, spot, ttm, v, sigma,
                                              kappa, theta, rho)/(1i*u)
    integrand2 <- Re(dummy2)
  }
  
  Pi1 <- 0.5 + 1/pi * integrate(integrand1, 0, 1000, stop.on.error = FALSE)$value
  Pi2 <- 0.5 + 1/pi * integrate(integrand2, 0, 1000, stop.on.error = FALSE)$value
  
  res <- spot*Pi1 - strike*Pi2
  return(res)
}

lossFunction <- function(parms){
  sum <- 0
  for(i in 1:l){
    sum <- sum + ( (data$Last[i] - 
                      callHestoncf(S0, data$Strike[i], data$TimeToExp,r=0.03, q = 0.01,
                                   parms[1], parms[2], parms[3], parms[4], parms[5])))^2
  }
  return(sum)
}
callHestoncf()

BSM <- function(S, K, sigma, ttm){
  d1 <- (log(S/K) + ((sigma^2)/2)*ttm)/(sigma*sqrt(ttm))
  d2 <- d1 - sigma*sqrt(ttm)
  price <- pnorm(d1)*S - K*pnorm(d2)
  
  return(price)
}

impVol <- function(price, S, K, ttm){
  f <- function(sigma){price -BSM(S,K,sigma,ttm)}
  imp_vol <- uniroot(f, c(0,5))
  
  return(imp_vol)
}

l <- length(data$TimeToExp)

par <- optim(par = c(0.1^2, 0.8, 2.01, 0.1^2, -0.75), lossFunction,
                method = "L-BFGS-B", lower = c(0.01, 0.01, 0.1, 0.01, -0.99),
                upper = c(0.99, 2, 5, 0.99, 0.99))

BSM(100, 110, 0.05, 5)



v <- 0.08118792
sigma <- 0.78055835
kappa <- 2.12469802
theta <- 0.05161678
rho <- -0.748924


heston_prices <- rep(0, l)
heston_vol <- rep(0,l)
for(i in 123:l){
  heston_prices[i] <- callHestoncf(S_t,data$Strike[i],data$TimeToExp[i],0.02,0,v,theta,rho,kappa,sigma)
}

heston_prices <- heston_prices[-(1:123)]
for(j in 1:508){
  heston_vol[j] <- impVol(heston_prices[j], S=S_t, data$strike[j+122], data$TimeToExp[j+122])
}
