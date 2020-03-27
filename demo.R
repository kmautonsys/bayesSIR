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

# outbreak priors
tranMU = 2.5 # mean Reproduction number (number new infections/case)
tranSD = 1.4 # standard deviation
T0 = numeric(0)
M = length(tranMU)
recMU = 14 # mean recovery time (days)
recSD = 1 # standard deviation

I0MU = 4000 # mean current infected population
I0SD = 500 # standard deviation
R0MU = 4000 # mean current recovered population
R0SD = 1000 # standard deviation

# forecast time points
ts = as.numeric(1:120)
N = length(ts)

# observed data
dat_ts = seq(-20,-1,length.out=5) # time of observations
Q = length(dat_ts)
dat_tests = rep(0,Q) # not used
dat_cases = 5000*exp(dat_ts*((3-1)/14)) # some very fake data
datSD = 100 # confidence in data

# population data
TotalPop=19e6 # total population

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
  datSD=datSD,
  TotalPop=TotalPop
  )

fit <- stan(
  file = "bayesSIRv1.0.stan",  
  data = data,    
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4               # number of cores (could use one per chain)
)

probs0 = c(0.25) # posterior percentiles
probs = c(probs0,0.5,rev(1-probs0))

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
####
df = as.data.frame(posteriorI)
# have to deal with the bad names from the summary 
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",1:length(probs0),"l"),"c50",rev(paste0("c",1:length(probs0),"u")))
# convert to % of total population
df[,colnames(df)!="t"] = df[,colnames(df)!="t"]/data$TotalPop*100

# now we're ready to plot
p = ggplot(df)
p = p+geom_line(aes(x=t,y=c50),size=1,color='dark blue')
p = p+geom_ribbon(aes(x=t,ymin=c25l,ymax=c25u),fill='blue',alpha=0.33)

p=p+theme_minimal(36)+xlab("Day")+ylab("% of population")
p=p+geom_point(data=data.frame(t=dat_ts,y=dat_cases/data$TotalPop*100),aes(x=t,y=y))
plot(p)

plot(p+xlim(-20,20)+ylim(0,0.1))

