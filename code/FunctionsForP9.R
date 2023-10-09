#Housekeeping 
rm(list = ls())
graphics.off()

#Packages
require(MASS)


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
  gamma <- function(k){
    return(0.5*((k+1)^(2*H) + abs(k-1)^(2*H) - 2*k^(2*H)))
  }
  C <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      C[i,j] <- gamma(abs(i-j))
    }
  }
  L <- t(chol(C))
  W <- rnorm(n, mean = 0, sd = 1)
  noise <- L%*%W
  noise <- c(0, noise)
  
  simulations <- cumsum(noise)
  simulations <- ((T/n)^(H))*simulations
  
  return(simulations)
}


#Simulating fBM using Circulant Embedding without FFT
simfBM_circulant <- function(m = 1, T, H, n){
  t_grid <- seq(0, T, T/n)
  
  #Define autocovariance function for fractional Gaussian noise and the matrix C
  gamma <- function(h){
    return(0.5*(abs(h-1)^(2*H) - 2*abs(h)^(2*H) + abs(h+1)^(2*H)))
  }
  
  C <- matrix(0, nrow = 2*n, 2*n)
  
  #Manually compute the first row of C
  for(i in 1:n){
    C[1,i] <- gamma(i-1)
  }
  
  for(i in 1:(n-1)){
    C[1,i+n+1] <- gamma(n-i)
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

#Function simulating fBM using circulant embedding with FFT
simfBM_fft <- function(m = 1, T, H, N){
  
  #Define autocovariance function for fractional Gaussian noise
  gamma <- function(k){
    return(0.5*((k+1)^(2*H) + abs(k-1)^(2*H) - 2*k^(2*H)))
  }
  
  # Get eigenvalues
  g = c()
  for(k in 0:(N-1)){g <- c(g, gamma(k))}
  r = c(g, 0, rev(g)[1:(N-1)])
  j = seq(0, ((2*N)-1))
  K = (2*N)-1
  i = complex(real=0, imaginary=1)
  lk = rev(fft(r*exp(2*pi*i*K*j*(1/(2*N)))))
  
  # Generate random variables
  Vj <- cbind(seq(0,0,length.out=2*N), seq(0,0,length.out=2*N))
  Vj[1,1] <- rnorm(1)
  Vj[N+1, 1] <- rnorm(1)
  Vj1 <- rnorm(N-1)
  Vj2 <- rnorm(N-1)
  Vj[2:N,1] <- Vj1
  Vj[2:N,2] <- Vj2
  Vj[(N+2):(2*N),1] <- rev(Vj1)
  Vj[(N+2):(2*N),2] <- rev(Vj2)
  
  # Compute Z (fractional Gaussian Noise)
  wk = seq(0,0,length.out=2*N)
  wk[1] <- sqrt(lk[1]/(2*N))*Vj[1,1]
  wk[2:N] <- sqrt(lk[2:N]/(4*N))*(Vj[2:N,1] + i*Vj[2:N,2])
  wk[N+1] <- sqrt(lk[N+1]/(2*N))*Vj[N+1,1]
  wk[(N+2):(2*N)] <- sqrt(lk[(N+2):(2*N)]/(4*N))*(Vj[(N+2):(2*N),1] - i*Vj[(N+2):(2*N),2])
  Z = fft(wk)
  fGn = Z[1:N]
  fBm = cumsum(fGn)*(N^(-H))
  fBM = Re((T^H)*fBm)
  return(fBM)
}

#Function for simulating price and volatility paths from Heston Model
heston_sim <- function(n_sims, T, kappa, theta, r, sigma, rho, v0, S0){
  dt <- T/n_sims
  t_grid <- seq(0,T,dt)
  
  S <- rep(0, n_sims + 1); S[1] <- S0
  V <- rep(0, n_sims + 1); V[1] <- v0
  
  cov_matrix <-matrix(data <- c(1,rho,rho,1), nrow=2,ncol=2,byrow=TRUE)
  B <- mvrnorm(n = n_sims, mu = c(0,0), Sigma = cov_matrix)
  B <- t(B)
  
  for(i in 1:n_sims){
    V_pos <- max(V[i],0)
    V[i+1] <- V[i] + kappa*(theta - V_pos)*dt + sigma*sqrt(V_pos)*sqrt(dt)*B[1,i]
    
    S[i+1] <- S[i]*exp( (r - 0.5*V[i])*dt + sqrt(V[i])*sqrt(dt)*B[2,i])
  }
  return(list(S,V))
}

