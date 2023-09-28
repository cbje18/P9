#Install dependencies first
!pip install eod
!pip install nelson_siegel_svensson

#Import libraries 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.optimize import minimize
from datetime import datetime as dt

from eod import EodHistoricalData
from nelson_siegel_svensson import NelsonSiegelSvenssonCurve
from nelson_siegel_svensson.calibrate import calibrate_nss_ols

#Set all parameters
S0 = 100.0  #Initial stock price
T = 1.0  #Time horizon
r = 0.02  #Risk-free interest rate
N = 252  #Number of points per path
M = 5  # Number of paths simulated

#Heston model parameters
kappa = 3 #Rate of reversion to theta
theta = 0.2 #Long-range average volatility
v0 = 0.25  #Initial volatility
rho = 0.7  #Correlation between BM's 
sigma = 0.6 #Volatility of volatility 

dt = T/N
mu = np.array([0,0])
cov = np.array([[1,rho], [rho,1]])

S = np.full(shape = (N+1,M), fill_value = S0)
V = np.full(shape = (N+1,M), fill_value = v0)

Z = np.random.multivariate_normal(mu, cov, (N,M))

for i in range(1,N+1):
  S[i] = S[i-1]*np.exp( (r-0.5*V[i-1])*dt + np.sqrt(V[i-1]*dt)*Z[i-1,:,0])
  V[i] = np.maximum(V[i-1] + kappa*(theta-V[i-1])*dt + sigma*np.sqrt(V[i-1]*dt)*Z[i-1,:,1],0)

#Plotting the paths
fig, (ax1, ax2)  = plt.subplots(1, 2, figsize=(12,5))
time = np.linspace(0,T,N+1)
ax1.plot(time, S)
ax1.set_title("Heston Model Asset Prices")
ax1.set_xlabel("Time")
ax1.set_ylabel("Price")

ax2.plot(time, V)
ax2.set_title("Heston Model Volatility Processes")
ax2.set_xlabel("Time")
ax2.set_ylabel("Volatility")

plt.show()
