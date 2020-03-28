#####################################################################
# Copyright 2020 AutonSystems, LLC
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#####################################################################

library(ggplot2)
library(rstan)
rstan_options(auto_write=TRUE)

# SOURCE: https://covidtracking.com/api/states/daily  NY state
# day 0 corresponds to 2020-03-28
dat_ts =  c(-24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10,  -9,  -8,  -7,  -6,  -5,  -4,  -3,  -2,  -1)
dat_cases = c(6,    22,    33,    76,   105,   142,   173,   216,   216,   421,   524,   729,   950,  1700,  2382,  4152,  7102, 10356, 15168, 20875, 25665, 30811, 37258, 44635)
dat_tests = c(54,   98,   125,  168,    197,    234,    265,    308,    308,   3200, 3303,   5272,   5493,   7206,  14597,  22284,  32427,  45437,  61401,  78289, 91270, 103479, 122104, 145753)
TotalPop = 19440469
Q = length(dat_ts)

# outbreak priors
tranMU = c(2.5) # mean Reproduction number (number new infections/case)
tranSD = c(1.4) # standard deviation
T0 = numeric(0)
M = length(tranMU)
recMU = 14 # mean recovery time (days)
recSD = 1 # standard deviation
SevereMU = 0.5; # mean proportion of COVID cases that are severe enough to garner testing
SevereSD = 0.25; # standard deviation
nonCOVIDMU = TotalPop*0.1 # mean number of nonCOVID cases that garner testing for COVID
nonCOVIDSD = 5e6 # standard deviation

I0MU = max(dat_cases) # mean current infected population
I0SD = 1e6 # standard deviation
R0MU = max(dat_cases) # mean current recovered population
R0SD = 1e6 # standard deviation

# forecast time points
ts = as.numeric(1:120)
N = length(ts)

data = list(
  M=M,
  tranMU=array(tranMU,dim=M),
  tranSD=array(tranSD,dim=M),
  T0=array(T0,dim=M-1),
  recMU=recMU,
  recSD=recSD,
  I0MU=I0MU,
  I0SD=I0SD,
  R0MU=R0MU,
  R0SD=R0SD,
  N=N,
  ts=array(ts,dim=N),
  Q=Q,
  dat_ts=array(dat_ts,dim=Q),
  dat_tests=array(dat_tests,Q),
  dat_cases=array(dat_cases,Q),
  SevereMU=SevereMU,
  SevereSD=SevereSD,
  nonCOVIDMU=nonCOVIDMU,
  nonCOVIDSD=nonCOVIDSD,
  TotalPop=TotalPop
)

fit <- stan(
  file = "bayesSIRv2.0.stan",  
  data = data,    
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4               # number of cores (could use one per chain)
)

probs0 = (5:48)/100 # posterior percentiles
probs = c(probs0,0.5,rev(1-probs0))

posteriorExpectedConfirmedCases = summary(fit, pars = "ExpectedConfirmedCases", probs = probs)$summary
posteriorExpectedConfirmedCases = cbind(posteriorExpectedConfirmedCases,dat_ts)  
colnames(posteriorExpectedConfirmedCases)[ncol(posteriorExpectedConfirmedCases)] = "t"

if(N>0){
  posteriorI = summary(fit, pars = "I", probs = probs)$summary
  if(Q>0){
    posteriorDatI = summary(fit, pars = "dat_I", probs = probs)$summary
    posteriorI = rbind(posteriorDatI,posteriorI)
    posteriorI = cbind(posteriorI,c(dat_ts,ts))
  } else{
    posteriorI = cbind(posteriorI,ts)  
  }
} else {
  posteriorI = summary(fit, pars = "dat_I", probs = probs)$summary
  posteriorI = cbind(posteriorI,dat_ts)  
}
colnames(posteriorI)[ncol(posteriorI)] = "t"
######################################

# Plot results

# plot forecast
####
df = as.data.frame(posteriorI)
# have to deal with the bad names from the summary 
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",probs0*100,"l"),"c50",rev(paste0("c",probs0*100,"u")))
# convert to % of total population
df[,colnames(df)!="t"] = df[,colnames(df)!="t"]/data$TotalPop*100

# now we're ready to plot
p = ggplot(df)
alpha = exp(seq(log(0.025),log(0.075),length.out=length(probs0)))
for( pr in probs0*100 ) p=p+geom_ribbon(aes_string(x="t",ymin=paste0("c",pr,"l"),ymax=paste0("c",pr,"u")),alpha=alpha[pr],fill="red")
p = p+geom_line(aes(x=t,y=c50),size=1)
p = p+geom_line(aes(x=t,y=c25l),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=c25u),size=0.5,lty=2)

p=p+theme_minimal(36)+xlab("Day")+ylab("% of population")
plot(p)

# plot fit
#####
df = as.data.frame(posteriorExpectedConfirmedCases)
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",probs0*100,"l"),"c50",rev(paste0("c",probs0*100,"u")))
p = ggplot(df)
alpha = exp(seq(log(0.025),log(0.075),length.out=length(probs0)))
for( pr in probs0*100 ) p=p+geom_ribbon(aes_string(x="t",ymin=paste0("c",pr,"l"),ymax=paste0("c",pr,"u")),alpha=alpha[pr],fill="red")
p = p+geom_line(aes(x=t,y=c50),size=1)
p = p+geom_line(aes(x=t,y=c25l),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=c25u),size=0.5,lty=2)
p=p+geom_point(data=data.frame(t=dat_ts,y=dat_cases),aes(x=t,y=y),size=2)
p=p+theme_minimal(36)+xlab("Day")+ylab("Count of cases")
plot(p)

####
