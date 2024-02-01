rm(list = ls())

#Function simulating fBM using circulant embedding with FFT
simfBM_fft <- function(m = 1, T, H, N){
  
  #Define autocovariance function for fractional Gaussian noise
  gamma <- function(k){return(0.5*((k+1)^(2*H) + abs(k-1)^(2*H) - 2*k^(2*H)))}
  
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
  fBM = c(0,fBM)
  return(fBM)
}

#Simulation Study for fOU
#Simulate Fractional Ornstein-Uhlenbeck by Euler Scheme
sim_fOU_euler <- function(T, N, alpha, H, phi, m, Y0){
  delta <- T/N
  t_grid <- seq(0, T, delta)
  
  W <- simfBM_fft(m = 1, T = T, H = H, N = N)
  dW <- diff(W, lag = 1)
  
  Y <- numeric(N+1)
  Y[1] <- Y0
  
  for(i in 2:(N+1)){
    Y[i] <- Y[i-1] + alpha*(m - Y[i-1])*delta + phi*(W[i]-W[i-1])
  }
  
  return(Y)
}

#We will now choose some parameters for our simulation study
m <- 3
alpha <- 0.79324
H <- 0.137923
phi <- 1.47829


#Simulate 1000 paths with 5000 points each
n_paths <- 2000
n_points <- 5000
sim_paths <- matrix(data = 0, nrow=n_points + 1, ncol = n_paths )

for(i in 1:n_paths){
  sim_paths[,i] <- sim_fOU_euler(T=1, N=n_points, alpha=alpha, H=H, phi=phi, m=m, Y0=3.5)
}

#Define the estimators as functions
H_hat <- function(x){
  N <- length(x)
  sum1 <- sum( (x[5:(N)] - 2*x[3:(N-2)] + x[1:(N-4)])^2)
  sum2 <- sum( (x[3:N] - 2*x[2:(N-1)] + x[1:(N-2)])^2)
  H_estimate <- abs(0.5*log2(sum1/sum2))
  
  return(H_estimate)
}

H_estimates <- rep(0, n_paths)
for(i in 1:n_paths){
  H_estimates[i] <- H_hat(sim_paths[,i])
}

mean(H_estimates)
hist(H_estimates, breaks = 25)
abline(v=H)
abline(v=mean(H_estimates), col="red")

phi_hat <- function(x, T=1){
  N <- length(x)
  delta <- 1/N
  tau <- 4*delta^(2*mean(H_estimates)) -(2*delta)^(2*mean(H_estimates))
  sum <- sum( (x[3:N] - 2*x[2:(N-1)] + x[1:(N-2)])^2)
  phi_estimate <- sqrt( (delta)/(T*tau)*sum)
  
  return(phi_estimate)
}

phi_estimates <- rep(0, n_paths)
for(i in 1:n_paths){
  phi_estimates[i] <- phi_hat(sim_paths[,i])
}

mean(phi_estimates)
hist(phi_estimates, breaks = 25)
abline(v=phi)

abline(v=mean(phi_estimates), col = "red")

m_hat <- function(x){
  N <- length(x)
  m_estimate <- (1/(N-1))*sum(x)
  
  return(m_estimate)
}

m_estimates <- rep(0,n_paths)
for(i in 1:n_paths){
  m_estimates[i] <- m_hat(sim_paths[,i])
}

mean(m_estimates)

alpha_hat <- function(x){
  N <- length(x)
  num <- N*sum(x^2) - sum(x)^2
  denom <- (N^2)*mean(phi_estimates)*mean(H_estimates)*gamma(2*mean(H_estimates))
  alpha_estimate <- (num/denom)^(-1/(2*mean(H_estimates)))
  
  return(alpha_estimate)
}

alpha_estimates <- rep(0, n_paths)
for(i in 1:n_paths){
  alpha_estimates[i] <- alpha_hat(sim_paths[,i])
}

mean(alpha_estimates) - 0.5


#Estimating parameters for SP-500 data 
data(rv5)
log_process <- log(sqrt(rv5$rv5)) %>% as.data.frame()
autoplot(log_process)

H_rv <- H_hat(log_process$rv5); H_rv

phi_hat2 <- function(x, T=5079,H, delta=1){
  N <- length(x)
  tau <- 4*delta^(2*H) - (2*delta)^(2*H)
  sum1 <- sum( (x[3:N] - 2*x[2:(N-1)] + x[1:(N-2)])^2)
  phi_estimate <- sqrt(sum1/(N*(4-2^(2*H_rv))*1^(2*H_rv) ))
  
  return(phi_estimate)
}


phi_rv <- phi_hat2(log_process$rv5, T = 5079, H =H_rv); phi_rv

m_rv <- (1/5078)*sum(log_process); m_rv

alpha_hat2 <- function(x, H, phi){
  N <- length(x)
  num <- N*sum(x^2) - sum(x)^2
  denom <- (N^2)*phi*H*gamma(2*H)
  alpha_estimate <- (num/denom)^(-1/(2*H))
  
  return(alpha_estimate)
}
alpha_rv <- alpha_hat2(log_process$rv5, H=H_rv, phi=phi_rv); alpha_rv

simulations <- sim_fOU_euler(T=1, N=5079, alpha=alpha_rv, H=H_rv, phi=0.38, m=m_rv, Y0=log_process$rv5[1])
simulations <- simulations
simulated_vol <- exp(simulations)
plot(simulated_vol, type = "l", ylim = c(0,0.061))

simulated_vol <- head(simulated_vol, -1)

t_grid <- seq(1,5079,1)
sim_data <- cbind(t_grid, simulated_vol) %>% as.data.frame()
ggplot(data = sim_data) + geom_line(aes(x=t_grid, y = simulated_vol))+ xlab(NULL) + ylab(NULL) + ylim(0,0.061) + ggtitle("Model")

alpha_rv
m_rv
H_rv
phi_rv+0.5
