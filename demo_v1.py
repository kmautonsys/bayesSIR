#!/usr/bin/env python3
import pystan
import numpy as np

# day 0 corresponds to 2020-03-28
dat = dict()
dat["dat_ts"] =  [-24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10,  -9,  -8,  -7,  -6,  -5,  -4,  -3,  -2,  -1]
dat["dat_cases"] = [6,    22,    33,    76,   105,   142,   173,   216,   216,   421,   524,   729,   950,  1700,  2382,  4152,  7102, 10356, 15168, 20875, 25665, 30811, 37258, 44635]
dat["dat_tests"] = [54,   98,   125,  168,    197,    234,    265,    308,    308,   3200, 3303,   5272,   5493,   7206,  14597,  22284,  32427,  45437,  61401,  78289, 91270, 103479, 122104, 14575]
dat["datSD"] = 100
TotalPop  = 19440469

# outbreak priors
dat["tranMU"] = np.array([2.5]) # mean Reproduction number (number new infections/case)
dat["tranSD"] = np.array([1.4]) # standard deviation

dat["recMU"] = 14 # mean recovery time (days)
dat["recSD"] = 1 # standard deviation

dat["SevereMU"] = 0.5; # mean proportion of COVID cases that are severe enough to garner testing
dat["SevereSD"] = 0.25; # standard deviation
dat["nonCOVIDMU"] = TotalPop*0.1 # mean number of nonCOVID cases that garner testing for COVID
dat["nonCOVIDSD"] = 5e6 # standard deviation
dat["ConfirmMU"] = 0
dat["ConfirmSD"] = 2

dat["I0MU"] = max(dat["dat_cases"]) # mean current infected population
dat["I0SD"] = 1e6 # standard deviation
dat["R0MU"] = max(dat["dat_cases"]) # mean current recovered population
dat["R0SD"] = 1e6 # standard deviation

dat["M"]  = len(dat["tranMU"])
dat["Q"]  = len(dat["dat_ts"])

dat["T0"] = np.array([])
dat["TotalPop"] = TotalPop

# Forecast parameters
ts = np.arange(1,120)
N  = len(ts)#dat["ts"])
dat["N"]  = N
dat["ts"] = ts

print(dat)
sm = pystan.StanModel(file="bayesSIRv1.1.stan")
fit = sm.sampling(data=dat, iter=1000, chains=4)


# if(N>0){
#   posteriorI = summary(fit, pars = "I", probs = probs)$summary
#   if(Q>0){
#     posteriorDatI = summary(fit, pars = "dat_I", probs = probs)$summary
#     posteriorI = rbind(posteriorDatI,posteriorI)
#     posteriorI = cbind(posteriorI,c(dat_ts,ts))
#   } else{
#     posteriorI = cbind(posteriorI,ts)  
#   }
# } else {
#   posteriorI = summary(fit, pars = "dat_I", probs = probs)$summary
#   posteriorI = cbind(posteriorI,dat_ts)  
# }
# colnames(posteriorI)[ncol(posteriorI)] = "t"

#if N>0:
#    posterior = fit.summary()
#    print(posterior)
