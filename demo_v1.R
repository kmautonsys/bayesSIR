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

# SOURCE: https://covidtracking.com/
library(jsonlite)
datdf <- fromJSON("https://covidtracking.com/api/states/daily")
datdf = datdf[datdf$state=="AS",]

dat_ts = as.Date(as.character(datdf$date),"%Y%m%d")
dat_ts = as.numeric(dat_ts-(max(dat_ts)+1))
dat_cases = as.integer(datdf$positive)
dat_hospitalizations = as.integer(datdf$hospitalized)
dat_deaths = as.integer(datdf$death)
idx = order(dat_ts); dat_ts=dat_ts[idx]; dat_cases=dat_cases[idx]; dat_hospitalizations=dat_hospitalizations[idx]; dat_deaths=dat_deaths[idx];
dat_hospNA = is.na(dat_hospitalizations)
dat_hospitalizations[dat_hospNA] = 0
dat_deathNA = is.na(dat_deaths)
dat_deaths[dat_deathNA] = 0
TotalPop = 19440469 #NY
Q = length(dat_ts)

# outbreak priors
T0 = seq(-14,-3,2)
M = length(T0)+1

# forecast time points
ts = as.numeric(1:60)
N = length(ts)

data = list(
  M=M,
  T0=array(T0,dim=M-1),
  N=N,
  ts=array(ts,dim=N),
  Q=Q,
  dat_ts=array(dat_ts,dim=Q),
  dat_cases=array(dat_cases,dim=Q),
  dat_caseNA=array(rep(0,Q),dim=Q),
  dat_hospitalizations=array(dat_hospitalizations,dim=Q),
  dat_hospNA=array(dat_hospNA*1,dim=Q), # NA deaths should be forward-filled
  dat_deaths=array(dat_deaths,dim=Q),
  dat_deathNA=array(dat_deathNA*1,dim=Q), # NA deaths should be forward-filled
  TotalPop=TotalPop,
  tranMU=array(log(2.5),dim=M),
  tranSD=array(1.0,dim=M),
  recMU=log(14),
  recSD=0.1,
  I0MU=0,
  I0SD=10,
  S0MU=log(0.5*TotalPop),
  S0SD=10,
  HospMU=-3,
  HospSD = 1,
  DeathMU = -4,
  DeathSD = 1,
  trendMU = 0,
  trendSD = 1.0,
  phi=-1
)

library(plyr)
#############
# MAP inference
vpar = function(par) c(par$S0_,par$I0_,par$rec_,par$tran_,par$HospProp_,par$DeathProp_)
stanmodel = stan_model("bayesSIRv1.1.stan")

optim_ = lapply(1:100,function(i){ 
  tryCatch({
  opt = optimizing(stanmodel, data = data, hessian=TRUE, as_vector=FALSE,iter=750)
  mxi = which.max(opt$par$I)
  mxv = opt$par$I[mxi]
  opt$mxi = mxi
  opt$mxv = mxv
  print(paste0(i,"  ",opt$value,"  ",opt$return_code))
  return(opt)
  },error=function(e){ return(NULL) },warning=function(w){})
})
optim_ = optim_[-which(unlist(lapply(optim_,function(o) is.null(o))))]
# coverage:
optim2 = list(); optim=optim_
for(i in 1:length(optim)){
  p = unlist(lapply(1:length(optim),function(k) optim[[k]]$value))
  j = which.max(p)
  optim2[[length(optim2)+1]] = optim[[j]]
  par_ = vpar(optim[[j]]$par)
  # drop similar
  drop_idx = which(unlist(lapply(optim,function(o) max(abs(vpar(o$par)-par_)/par_)<0.05 )))
  optim[drop_idx] = NULL
  if(length(optim)==0) break
}
best = optim2

v=unlist(lapply(best,function(o) o$value)); mx=max(v);
best = best[which(unlist(lapply(best,function(o) o$value>mx+log(0.01) )))] # mx+log(tau) to keep samples with probability mass at least tau-% of the most likely point
v=unlist(lapply(best,function(o) o$value)); mx=max(v);
exp(v-mx)/sum(exp(v-mx))

print(length(best))
print(unlist(lapply(best,function(o) o$value)))

