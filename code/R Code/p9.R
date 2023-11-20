
ImpliedVol = function(S0,tau,K,Surface,mu,t){
  VolSurface = matrix(0,length(tau),length(K))
  for (i1 in 1:length(tau)){
    for (i2 in 1:length(K)){
      if(Surface[i1,i2] == 0){ next }
      VolSurface[i1,i2] = uniroot(FitPrice,c(-100,100), S0 = S0, K = K[i2], P = Surface[i1,i2], mu = mu, tau = tau[i1], t = t)$root
    }
  }
  return(VolSurface)
}

OUFit = function(V,H,Tau,K,S0,n,t,Skew){
  if(H[5] > 0.5){return(Inf)}
  if(H[5] < 0){return(Inf)}
  if(H[6] > 1){return(Inf)}
  if(H[6] < -1){return(Inf)}
  set.seed(2)
  S = SampleOU(S0,H[1],H[2],H[3],H[4],H[5],H[6],n,t,Tau)
  
  P = matrix(0,length(Tau),length(K))
  
  for(i in 1:length(Tau)){
    P[i,] = Price(S[i,],K)
  }
  V.est = ImpliedVol(S0,Tau,K,P,0,t)
  print(H)
  print(sum((t(V.est)[V!=0]-V[V!=0])^2))
  return(sum((t(V.est)[V!=0]-V[V!=0])^2)+10000*sum((Skew - FitPowerLaw(Tau/t,abs((V.est[,11]-V.est[,6]))/(log(K[11])-log(K[6])),Plot=FALSE)[[1]])^2))
}

FitModel = function(V,H,Tau,K,S0,n,t,Skew,model = "BS"){
  Model = optim(H,OUFit,V = V,Tau = Tau, K = K, S0 = S0, n = n, t = t, Skew = Skew,control = list(reltol=0.01))$par
  return(Model)
}

EstimateVol = function(Data,S0,Tau,K,t){
  Vol = matrix(0,length(K),length(Tau))
  for(i1 in 1:length(Tau)){
    for(i2 in 1:length(K)){
      n = 0
      for(i3 in 1:length(S0)){
        if((S0[i3] != 0)&(Data[i3+length(S0)*(i2-1),i1] != 0)){
          Vol[i2,i1] = Vol[i2,i1] + uniroot(FitPrice,c(-100,100), S0 = S0[i], K = K[i2], P = Data[i3+length(S0)*(i2-1),i1], mu = 0, tau = Tau[i1], t = t)$root
          n = n + 1
        }
      }
      Vol[i2,i1] = Vol[i2,i1] / (n + (n==0))
    }
  }
  return(Vol)
}

Skew = function(Xl,Xh,S0,Tau,t,K){
  Y = numeric(length(Tau))
  for(j in 1:length(Tau)){
    n = 0
    for(i in 1:length(S0)){
      if((S0[i] != 0)&(Xh[i,j] != 0)&(Xl[i,j] != 0)){
        Y[j] = Y[j] + uniroot(FitPrice,c(0,100), S0 = S0[i], K = K[2], P = Xh[i,j], mu = 0, tau = Tau[j], t = t)$root
        Y[j] = Y[j] - uniroot(FitPrice,c(0,100), S0 = S0[i], K = K[1], P = Xl[i,j], mu = 0, tau = Tau[j], t = t)$root
        n = n + 1
      }
    }
    Y[j] = Y[j] / (n+(n==0))
  }
  return(Y/(log(K[2])-log(K[1])))
}

FitPowerLaw = function(x,y,titel = "",Plot=TRUE){
  if( Plot ){ plot(x,y,xlab = expression(tau),ylab=expression(psi(tau)),main = titel) }
  s = sign(mean(exp(-c(1:length(y))/2)*y));  Y = y*s;  X = cbind(1,log(x[(Y>0)]))
  lines(x[(Y>0)],s*exp(X%*%solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])))
  b = solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])
  return(list(c(exp(b[1])*s,b[2]),log(Y[(Y>0)])-X%*%solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])))
}

FitPrice = function(S0,K,mu,P,tau,t,sigma0){
  if (sigma0 == 0){
    d1 = Inf
    d2 = Inf
  }
  else{
    d1 = (log(S0/K)+tau/t*(mu+sigma0^2/2))/(sigma0*sqrt(tau/t))
    d2 = (log(S0/K)+tau/t*(mu-sigma0^2/2))/(sigma0*sqrt(tau/t))
  }
  price = pnorm(d1)*S0-pnorm(d2)*K*exp(-mu*tau/t)
  return(price-P)
}

SampleOU = function(S0,Y0,lambda,theta,nu,H,rho,n,t,Tau){
  time = Sys.time()
  Sigma = matrix(0,t,n); dZ = matrix(0,t,n)
  for(i in 1:n){ D = OU(Y0,lambda,theta,nu,H,t); Sigma[,i] = D[[1]]
  dZ[,i] = rho*D[[2]]+sqrt(1-rho^2)*rnorm(t)/sqrt(t) }
  Sigma = exp(Sigma)
  S = BSV(S0,Sigma,dZ,Tau)
  if(length(Tau)>1){
    sigma = lm(log(S%*%rep(1,length(S[1,])))~Tau)$coefficients[2]*Tau%*%t(rep(1,n))
    S = S*exp(-sigma)
  }
  else{S = S/mean(S)}
  print(Sys.time()-time)
  return(S)
}


OU = function(Y0,lambda,theta,nu,H,t){
  if(H == 0.5){ dW = rnorm(t)/sqrt(t) }else{ D = fBm(H,t); dW = D[[1]]; dWW = D[[2]] }
  Y = numeric(t); Y[1] = Y0
  for(i in 1:t){ Y[i] = Y[i-(i>1)] - lambda*(Y[i-(i>1)]-theta)/t+nu*dW[i] }
  if(H<0.5){ return(list(Y,dWW)) }
  else { return(list(Y,dW)) }
}

SampleHes = function(S0,V0,kappa,theta,xi,rho,n,t,Tau){
  Sigma = matrix(0,t,n); dZ = matrix(0,t,n)
  for(i in 1:n){ D = HesVol(V0,kappa,theta,xi,t); Sigma[,i] = D[[1]]
  dZ[,i] = rho*D[[2]]+sqrt(1-rho^2)*rnorm(t)/sqrt(t) }
  S = BSV(S0,Sigma,dZ,Tau)
  sigma = lm(log(S%*%rep(1,length(S[1,])))~Tau)$coefficients[2]*Tau%*%t(rep(1,n))
  S = S*exp(-sigma)
  return(S)
}

HesVol = function(V0,kappa,theta,xi,t){
  dW = rnorm(t)/sqrt(t);  V = numeric(t);  V[1] = V0
  for(i in 1:t){ V[i] = abs(V[i-(i>1)] - kappa*(V[i-(i>1)]-theta)/t+xi*sqrt(V[i-(i>1)])*dW[i]) }
  return(list(sqrt(V),dW))
}

BSV = function(S0,sigma,dZ,Tau){
  t = length(dZ[,1])
  n = length(dZ[1,])
  Tau = c(0,Tau)
  S = matrix(0,length(Tau)-1,n)
  for (i1 in 2:length(Tau)){
    S[i1-1,] = S[i1-2+(i1==2),] + rep(1,Tau[i1]-Tau[i1-1])%*%(sigma[(Tau[i1-1]+1):Tau[i1],]*dZ[(Tau[i1-1]+1):Tau[i1],]
                                                              - 0.5*sigma[(Tau[i1-1]+1):Tau[i1],]^2/t)
  }
  S = S0*exp(S)
  return(S)
}

Price = function(S,K){
  price = numeric(length(K))
  for(i in 1:length(K)){
    price[i] = mean((abs(S-K[i])+S-K[i])/2)
  }
  return(price)
}

fBm = function(hurst,n){
  n = 2*(n - 1)
  r <- numeric(n+1)
  r[1] <- 1
  for (k in 1:n) r[k + 1] <- 0.5 * ((k + 1)^(2 * hurst) - 
                                      2 * k^(2 * hurst) + (k - 1)^(2 * hurst))
  r <- c(r, r[seq(length(r) - 1, 2)])
  lambda <- Re((fft(r)))/(2 * n)
  dW = Re(fft(sqrt(lambda) * (1+1i)*rnorm(2 * n)))[1:(n + 2)]
  W <- n^(-hurst) * (dW[2*(1:(n%/%2+1))-1]+dW[2*(1:(n%/%2+1))])
  return(list(W,dW[2*(1:(n%/%2+1))]/sqrt(n%/%2 + 1)))
}

fBmCov = function(H,t){
  Cov = matrix(0,t,t)
  for(i1 in 1:t){
    Cov[i1,i1] = (i1/t)^(2*H)
    if(i1 == t){next}
    for(i2 in (i1+1):t){
      Cov[i1,i2] = 1/2*((i1/t)^(2*H)+(i2/t)^(2*H)-abs((i1-i2)/t)^(2*H))
      Cov[i2,i1] = 1/2*((i1/t)^(2*H)+(i2/t)^(2*H)-abs((i1-i2)/t)^(2*H))
    }
  }
  return(Cov)
}


fBm2 = function (H, n,lambda=NA){
  if(is.na(lambda[1])){
    r <- numeric(n + 1)
    r[1] <- 1
    for (k in 1:n) r[k + 1] <- 0.5 * ((k + 1)^(2 * H) - 
                                        2 * k^(2 * H) + (k - 1)^(2 * H))
    r <- c(r, r[seq(length(r) - 1, 2)])
    lambda <- Re((fft(r))/(2 * n))
  }
  W <- fft(sqrt(lambda) * (rnorm(2 * n) * (1+1i)))
  W <- n^(-H) * cumsum(Re(W[1:(n + 1)]))
  return(W)
}

ComparefBm = function(H,n,K){
  t.Exact.K = numeric(length(K));  t.Circ.K = numeric(length(K))
  Error.Exact.K = numeric(length(K));  Error.Circ.K = numeric(length(K))
  for(j in 1:length(K)){
    print(K[j])
    Cov = fBmCov(H,K[j])
    mat.Exact = matrix(0,K[j],K[j])
    mat.Circ = matrix(0,K[j],K[j])
    tCov = t(chol(Cov))
    r <- numeric(K[j])
    r[1] <- 1
    for (k in 1:(K[j]-1)) r[k + 1] <- 0.5 * ((k + 1)^(2 * H) - 
                                               2 * k^(2 * H) + (k - 1)^(2 * H))
    r <- c(r, r[seq(length(r) - 1, 2)])
    lambda <- Re((fft(r))/(2 * (K[j]-1)))
    for(i in 1:n){
      time = Sys.time()
      W = tCov%*%rnorm(K[j])
      t.Exact.K[j] = t.Exact.K[j] + (Sys.time() - time)[[1]]
      mat.Exact = mat.Exact + W%*%t(W)
      time = Sys.time()
      W = fBm2(H,K[j]-1,lambda)
      t.Circ.K[j] = t.Circ.K[j] + (Sys.time() - time)[[1]]
      mat.Circ = mat.Circ + W%*%t(W)
    }
    Error.Exact.K[j] = sum((Cov-mat.Exact/n)^2)/(K[j]^2)
    Error.Circ.K[j] = sum((Cov-mat.Circ/n)^2)/(K[j]^2)
  }
  return(rbind(K,sqrt(rbind(t.Exact.K^2/n,Error.Exact.K,t.Circ.K^2/n,Error.Circ.K)/n)))
}

PriceError = function(Y0,lambda,theta,nu,H,rho,t,n=100000){
  N = c(10,20,50,100,200,500,1000,2000,5000,10000)
  M = numeric(length(N));  Mean = numeric(length(N));  Err = numeric(length(N))
  Time = numeric(length(N)); Time.std = numeric(length(N));
  for(i in 1:length(N)){
    m = 0;  M1 = 0;  M2 = 0; T1 = 0; T2 = 0
    print(N[i])
    while(m<n/N[i]){
      time = Sys.time()
      val = Price(SampleOU(1,Y0,lambda,theta,nu,H,rho,N[i],t,t),1)
      #val = uniroot(FitPrice,c(0,100), S0 = 1, K = 1, P = val, mu = 0, tau = t, t = t)$root
      time = (Sys.time() - time)[[1]]
      m = m + 1;  M1 = M1 + val;  M2 = M2 + val^2;  T1 = T1 + time;  T2 = T2 + time^2
    }
    M[i] = m;  Mean[i] = M1/m;  Err[i] = sqrt(M2/m - M1^2/m^2);  Time[i] = T1/m; Time.std[i] = sqrt(T2/m-T1^2/m^2)
  }
  return(list(N,M,Mean,Err,Time,Time.std))
}

MCDist = function(Y0,lambda,theta,nu,H,rho,t){
  val = numeric(100)
  for(i in 1:100){
    val[i] = Price(SampleOU(1,Y0,lambda,theta,nu,H,rho,10000,t,t),1)
    print(i)
  }
  return(val)
}


TestOU = function(Y0,lambda,theta,nu,dW,Type){
  t = length(dW)
  Y = numeric(t)
  time = Sys.time()
  if(Type=="Euler"){ Y[1] = Y0
  for(i in 1:t){
    Y[i] = Y[i-(i>1)] - lambda*(Y[i-(i>1)]-theta)/t+nu*dW[i]
  }
  }
  else{
    for(i in 1:t){
      Y[i] = Y[i-(i>1)] + exp(lambda*i/t)*dW[i]
    }
    Y = nu*exp(-lambda*(1:t)/t)*Y + theta + (Y0-theta)*exp(-lambda*(1:t)/t)
  }
  return(list(Y,(Sys.time()-time)[[1]]))
}


CompareOU = function(Y0,lambda,theta,nu,H,t,n,K){
  t.Exact = 0;  t.Exact.K = numeric(length(K))
  t.Euler = 0;  t.Euler.K = numeric(length(K))
  Error.Exact = 0;  Error.Exact.K = numeric(length(K))
  Error.Euler = 0;  Error.Euler.K = numeric(length(K))
  for(i in 1:n){
    dW = fBm(H,t)[[1]]
    OU.Exact = TestOU(Y0,lambda,theta,nu,dW,"Exact")
    t.Exact = t.Exact + OU.Exact[[2]]
    OU.Exact = OU.Exact[[1]]
    OU.Euler = TestOU(Y0,lambda,theta,nu,dW,"Euler")
    t.Euler = t.Euler + OU.Euler[[2]] 
    OU.Euler = OU.Euler[[1]]
    Error.Euler = Error.Euler + sum((OU.Exact-OU.Euler)^2)/t
    for(j in 1:length(K)){
      t.k = t%/%K[j]
      dW.k = diff(c(0,cumsum(dW)[K[j]*(1:t.k)]))
      OU.K = TestOU(Y0,lambda,theta,nu,dW.k,"Euler")
      t.Euler.K[j] = t.Euler.K[j] + OU.K[[2]]
      OU.K = OU.K[[1]]
      Error.Euler.K[j] = Error.Euler.K[j] + sum((OU.Exact[K[j]*(1:t.k)]-OU.K)^2)/t.k
      if(i == n){
        plot((1:t.k)/t.k,OU.Exact[K[j]*(1:t.k)],type="l",xlab="t",ylab="Y",main=paste0("Euler-Maruyama, n = ",t.k));  lines((1:t.k)/t.k,OU.K,col="red") }
      OU.K = TestOU(Y0,lambda,theta,nu,dW.k,"Exact")
      t.Exact.K[j] = t.Exact.K[j] + OU.K[[2]]
      OU.K = OU.K[[1]]
      Error.Exact.K[j] = Error.Exact.K[j] + sum((OU.Exact[K[j]*(1:t.k)]-OU.K)^2)/t.k
      if(i == n){
        plot((1:t.k)/t.k,OU.Exact[K[j]*(1:t.k)],type="l",xlab="t",ylab="Y",main=paste0("Exact, n = ",t.k));  lines((1:t.k)/t.k,OU.K,col="red") }
    }
  }
  return(sqrt(rbind(c(t,t/K)^2*n,c(t.Exact,t.Exact.K)^2/n,c(Error.Exact,Error.Exact.K),c(t.Euler,t.Euler.K)^2/n,c(Error.Euler,Error.Euler.K))/n))
}

CompareRFSV = function(S0,Y0,lambda,theta,nu,H,rho,K,n){
  t.K = numeric(length(K))
  Error.K = numeric(length(K))
  for(i in 1:length(K)){
    set.seed(2)
    time = Sys.time()
    S.k = SampleOU(S0,Y0,lambda,theta,nu,H,rho,n,K[i],1:K[i])
    P.k = matrix(0,K[i],1)
    for(j in 1:K[i]){ P.k[j,] = Price(S.k[j,],S0) }
    t.K[i] = (Sys.time() - time)[[1]]
    if(i == 1){ P = P.k; plot((1:K[i])/K[i],P,type="l",xlab="t",ylab="price",main="Prices for Different Step Sizes");  next }
    lines((1:K[i])/K[i],P.k,col=i)
    Error.K[i] = Error.K[i] + sum((P[K[i-1]/K[i]*(1:K[i])]-P.k)^2)/K[i]
    P = P.k
  }
  return(rbind(t.K,sqrt(Error.K)))
}

GetData = function(ticker,time,maturity,strikes){
  #CvQbcdE7syVkOIP5dAqccpOeK2o0iup4
  #CBQouP6k8C92g2XWhofZB5PCGgxI6lOk
  API = paste0("https://api.polygon.io/v2/aggs/ticker/",ticker, "/range/1/minute/")
  API = paste0(API,format(time,"%Y-%m-%d"),"/",format(time,"%Y-%m-%d"))
  API = paste0(API,"?adjusted=true&sort=asc&limit=5000&apiKey=CBQouP6k8C92g2XWhofZB5PCGgxI6lOk")
  raw = jsonlite::fromJSON(API)
  Data = raw$results
  minutes = Data$t/60000
  Asset = numeric(minutes[length(minutes)]-minutes[1]+1)
  Asset[minutes-minutes[1]+1] = Data$vw
  center = minutes[1]-1+330
  strikes = c(-5*(5:1)+5*floor(mean(Data$vw)/5),(5*floor(mean(Data$vw)/5)):(5*ceiling(mean(Data$vw)/5)),
              5*ceiling(mean(Data$vw)/5)+5*(1:5))
  print(strikes)
  n = (minutes[length(minutes)]-240-center)
  Data.Opt = matrix(0,length(strikes)*n,maturity)
  Maturities = rep(FALSE,maturity)
  for(i in 1:maturity){
    I = time + i
    if(weekdays(I) == "Saturday"){ next }
    if(weekdays(I) == "Sunday"){ next }
    Sys.sleep(12)
    Data = APICall(ticker,I,5*floor(mean(Asset[minutes-minutes[1]+1])/5),time)
    if(length(Data)>1){
      Maturities[i] = TRUE
      print(format(I,"%y%m%d"))
    }
    if(Maturities[i] == TRUE){
      for(j in 1:length(strikes)){
        Sys.sleep(12)
        Data = APICall(ticker,I,strikes[j],time)
        if(length(Data)>1){
          Data.Opt[j*n-n+Data$t/60000-center,i] = Data$vw
          print(strikes[j])
        }
      }
    }
  }
  return(list(Asset,(1:maturity)[Maturities],Data.Opt[,Maturities]))
}

APICall = function(ticker,I,price,time){
  API.opt = paste0("https://api.polygon.io/v2/aggs/ticker/O:",ticker)
  API.opt = paste0(API.opt,format(I,"%y%m%d"))
  API.opt = paste0(API.opt,"C00",price,"000/range/1/minute/")
  API.opt = paste0(API.opt,format(time,"%Y-%m-%d"),"/",format(time,"%Y-%m-%d"))
  API.opt = paste0(API.opt,"?adjusted=true&sort=asc&limit=5000&apiKey=CBQouP6k8C92g2XWhofZB5PCGgxI6lOk")
  raw = jsonlite::fromJSON(API.opt)
  return(raw$results)
}