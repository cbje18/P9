#Housekeeping 
rm(list = ls())
graphics.off()


#Function for simulating a Brownian motion
#inputs are dimension m, time T, number of points n
simBM <- function(m=1, T, n){
  delta = T/n
  tmesh <- seq(0,T, delta)
  N <- length(tmesh)
  B <- matrix(0, N, 1)
  for(i in 2:N){
    B[[i,1]] <- B[[i-1,1]] + sqrt(T/n)*rnorm(m, mean = 0, sd = 1)
  }
  return(B)
}


#Simulate an Ornstein-Uhlenbeck process
#Parameters are x0, T, lambda, sigma, and number of points n
simOU <- function(x0, T, lambda, sigma, n){
  tmesh <- seq(0, T, T/n)
  N <- length(tmesh)
  OU <- matrix(x0, N, 1)
  for(i in 2:N){
    xi <- rnorm(1, mean = 0, sd = 1)
    OU[[i,1]] <- exp(-lambda*(T/n))*OU[[i-1,1]] + sigma*sqrt(1/(2*lambda)*(1-exp(-2*lambda*(T/n))))*xi
  }
  return(OU)
}


#Function which computes the Black-Scholes price of a European call option
#with parameters sigma, time to maturity ttm, rate r, strike price K, and stock price St
BlackScholes <- function(sigma, ttm, r, K, St){
  t <- ttm/365.25
  d1 <- (log(St/K) + (r+sigma^2 / 2)*t)/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  BS_price <- St*pnorm(d1) - exp(r*t)*K*pnorm(d2)
  
  return(BS_price)
}

#Objective function for calculating the implied volatility 
ImVolFun <- function(sigma, market_price, St, K, r, ttm){
  abs(market_price - BlackScholes(sigma, ttm, r, K, St))
}

#Wrapper for the optimization function
getIV <- function(x, St, r){
  result <- optimize(ImVolFun, interval = c(0,2),
                     market_price = as.numeric(x["ask"]), 
                     St = St,
                     K = as.numeric(x["strike"]),
                     r = r,
                     ttm = as.numeric(x["ttm"]))
  
  return(result$minimum)
}

#Function simulating a fractional BM using the Cholesky method
simfBM_cholesky <- function(m=1, T, H, n){
  t_grid <- seq(0, T, T/n)
  cov_function <- function(t,s){
    return(0.5*(abs(t)^(2*H) + abs(s)^(2*H) - abs(t-s)^(2*H)))
  }
  C <- matrix(0, nrow = length(t_grid), ncol = length(t_grid))
  for(i in 1:length(t_grid)){
    for(j in 1:length(t_grid)){
      C[i,j] <- cov_function(t_grid[i], t_grid[j])
    }
  }
  C <- C[-1,-1]
  L <- t(chol(C))
  W <- rnorm(n, mean = 0, sd = 1)
  simulations <- L%*%W
  simulations <- c(0, simulations)
  
  return(simulations)
}


#Simulating fBM using Circulant Embedding without FFT
simfBM_Circulant <- function(m = 1, T, H, n){
  t_grid <- seq(0, T, T/n)
  
  #Define autocovariance function for fractional Gaussian noise and the matrix C
  cov_noise <- function(h){
    return(0.5*(abs(h-1)^(2*H) - 2*abs(h)^(2*H) + abs(h+1)^(2*H)))
  }
  
  C <- matrix(0, nrow = 2*n, 2*n)
  
  #Manually compute the first row of C
  for(i in 1:n){
    C[1,i] <- cov_noise(i-1)
  }
  
  for(i in 1:(n-1)){
    C[1,i+n+1] <- cov_noise(n-i)
  }
  
  #Function for creating circulant matrix
  circulant <- function(x) {
    n <- length(x)
    C <- matrix(NA, n, n)
    
    for (i in 1:n) {
      C[i, ] <- c(x[-(1:(n + 1 - i))], x[1:(n + 1 - i)])
    }
    return(C)
    
  }
  
  #Finally creating the matrix C
  C <- circulant(C[1,])
  
  #Now onto the unitary matrix Q
  Q <- matrix(0, nrow = 2*n, ncol = 2*n)
  
  for(j in 1:(2*n)){
    for(k in 1:(2*n)){
      Q[j,k] <- (1/sqrt(2*n))*exp(complex(real = 0,
                                          imaginary = -2*pi*((j-1)*(k-1))/(2*n)))
    }
  }
  
  Q_star <- Conj(t(Q))
  
  #Computing eigenvalues and the matrix D
  lambdas <- rep(0, 2*n)
  
  for(k in 1:(2*n)){
    summands <- rep(0, 2*n)
    for(j in 1:(2*n)){
      summands[j] <-  C[1,j]*exp(complex(real = 0,
                                         imaginary = 2*pi*(j-1)*(k-1)/(2*n))) 
    }
    lambdas[k] <- sum(summands)
  }
  
  #All eigenvalues are real and positive 
  lambdas <- Re(lambdas)
  D <- diag(x = lambdas, nrow = 2*n, ncol = 2*n)
  
  #Actual simulation procedure
  V <- rnorm(n = 2*n, mean = 0, sd = 1)
  W <- Q_star%*%V
  
  D_sqrt <- sqrt(D)
  Z <- Q%*%D_sqrt%*%W 
  Z <- Re(Z)
  
  simulations <- c(0, cumsum(Z[1:n+1]))
  simulations <- n^(-H)*simulations
  simulations <- (T^H)*simulations
  return(simulations)
}





