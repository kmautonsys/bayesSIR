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
library(MASS)
rstan_options(auto_write=TRUE)

# SOURCE: https://covidtracking.com/
library(jsonlite)
datdf <- fromJSON("https://covidtracking.com/api/states/daily")
state="MI"
datdf = datdf[datdf$state==state,]

dat_ts = as.Date(as.character(datdf$date),"%Y%m%d")
dat_ts = as.numeric(dat_ts-(max(dat_ts)+1))
dat_cases = as.integer(datdf$positive)
idx = order(dat_ts); dat_ts=dat_ts[idx]; dat_cases=dat_cases[idx];
#TotalPop = 19440469 #NY
#TotalPop = 10045029 #MI
#TotalPop = 39937489 #CA
#TotalPop = 12820878 #PA
location <- fromJSON("https://coronadatascraper.com/locations.json")
state_string = sprintf("%s, USA",state);
state_data = location[location$name==state_string,]$population
TotalPop = sum(state_data)

Q = length(dat_ts)

# outbreak priors
tranMU = c(2.5,2.5) # mean Reproduction number (number new infections/case)
tranSD = c(1.4,1.4) # standard deviation
T0 = c(-12) #numeric(0)
M = length(tranMU)

# forecast time points
ts = as.numeric(1:60)
N = length(ts)

data = list(
  M=M,
  tranMU=array(tranMU,dim=M),
  tranSD=array(tranSD,dim=M),
  T0=array(T0,dim=M-1),
  recMU=14,
  recSD=1,
  I0MU=max(dat_cases),
  I0SD=1e6,
  R0MU=max(dat_cases),
  R0SD=1e6,
  N=N,
  ts=array(ts,dim=N),
  Q=Q,
  dat_ts=array(dat_ts,dim=Q),
  dat_cases=array(dat_cases,Q),
  TotalPop=TotalPop,
  datSD = 100,
  ConfirmMU = 0,
  ConfirmSD = 2
)

library(plyr)
#############
# MAP inference
stanmodel = stan_model("bayesSIRv1.1.stan")

optim_ = lapply(1:500,function(i){ 
  tryCatch({
  opt = optimizing(stanmodel, data = data)
  idx = grep("^I\\[[0-9]+\\]$",names(opt$par))
  mxi = which.max(opt$par[idx])
  mxv = opt$par[idx][mxi]
  opt$mxi = mxi
  opt$mxv = mxv
  print(paste0(i,"  ",opt$value))
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
  par_ = optim[[j]]$par[1:(5+M)] # include ConfirmProportion_, since ConfirmProportion is exponentiated
  # drop similar
  drop_idx = which(unlist(lapply(optim,function(o) max(abs(o$par[1:(5+M)]-par_)/par_)<0.05 )))
  optim[drop_idx] = NULL
  if(length(optim)==0) break
}
# coverage:
best = list(); optim = optim2
for(i in 1:length(optim)){
  p = unlist(lapply(1:length(optim),function(k) optim[[k]]$value))
  j = which.max(p)
  best[[length(best)+1]] = optim[[j]]
  mxi = optim[[j]]$mxi
  mxv = optim[[j]]$mxv
  # drop similar
  drop_idx = which(unlist(lapply(optim,function(o) abs(o$mxi-mxi)<=1 & abs(o$mxv-mxv)<=TotalPop*0.01 )))
  best[[length(best)]]$value = best[[length(best)]]$value+log(sum(exp(p[drop_idx]-best[[length(best)]]$value)))
  optim[drop_idx] = NULL
  if(length(optim)==0) break
}
v=unlist(lapply(best,function(o) o$value)); mx=max(v);
best = best[which(unlist(lapply(best,function(o) o$value>mx+log(0.001) )))] # mx+log(tau) to keep samples with probability mass at least tau-% of the most likely point

v=unlist(lapply(best,function(o) o$value)); mx=max(v);
exp(v-mx)/sum(exp(v-mx))

print(length(best))
print(unlist(lapply(best,function(o) o$value)))
print(best[[1]]$par[grep("^tran\\[[0-9]+\\]$",names(best[[1]]$par))])
print(best[[1]]$par[grep("^rec$",names(best[[1]]$par))])
print(best[[1]]$par[grep("^ConfirmProportion_$",names(best[[1]]$par))])
#############
df = ldply(best,function(opt){
  idx = grep("^I\\[[0-9]+\\]$",names(opt$par))
  return(data.frame(LogProb=opt$value,y=opt$par[idx],t=ts))
})
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases")
###
df = ldply(best,function(opt){
  idx = grep("^FitCases\\[[0-9]+\\]$",names(opt$par))
  return(data.frame(LogProb=opt$value,y=opt$par[idx],t=dat_ts))
})
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=dat_ts,y=dat_cases),aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases")

#############
opt = best[[1]]
init = Reduce(c,lapply(names(opt$par),function(n){ 
  l = list(opt$par[n])
  names(l)=n
  return(l)
}))
hess_inv = solve(-optimizing(stanmodel, data = data, init=init, hessian=TRUE)$hessian)
mu = sapply(1:nrow(hess_inv),function(i) init[[i]])

samp = lapply(1:100,function(i) as.list(mvrnorm(n = 1, mu, hess_inv)))
fit <- stan(file = "bayesSIRv1.1.stan",data=data,init=samp,chains=length(samp),warmup=0,iter=1,cores=1,algorithm="Fixed_param",refresh = 0)

probs0 = (5:49)/100 # posterior percentiles
probs = c(probs0,0.5,rev(1-probs0))

FitCases = summary(fit, pars = "FitCases", probs = probs)$summary
FitCases = cbind(FitCases,dat_ts)  
colnames(FitCases)[ncol(FitCases)] = "t"

posteriorI = summary(fit, pars = "I", probs = probs)$summary
posteriorI = cbind(posteriorI,ts)  
colnames(posteriorI)[ncol(posteriorI)] = "t"
######################################

# Plot results
####

# plot forecast
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

p = p+geom_line(data=data.frame(y=opt$par[grep("^I\\[[0-9]+\\]$",names(opt$par))]/TotalPop*100,t=ts),aes(x=t,y=y),col='blue',size=1) 

p=p+theme_minimal(36)+xlab("Day")+ylab("% of population")
plot(p)

# Plot fit to data
df = as.data.frame(FitCases)
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",probs0*100,"l"),"c50",rev(paste0("c",probs0*100,"u")))
p = ggplot(df)
alpha = exp(seq(log(0.025),log(0.075),length.out=length(probs0)))
for( pr in probs0*100 ) p=p+geom_ribbon(aes_string(x="t",ymin=paste0("c",pr,"l"),ymax=paste0("c",pr,"u")),alpha=alpha[pr],fill="red")
p = p+geom_line(aes(x=t,y=c50),size=1)
p = p+geom_line(aes(x=t,y=c25l),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=c25u),size=0.5,lty=2)
p=p+geom_point(data=data.frame(t=dat_ts,y=dat_cases),aes(x=t,y=y),size=2)

p = p+geom_line(data=data.frame(y=opt$par[grep("^FitCases\\[[0-9]+\\]$",names(opt$par))],t=dat_ts),aes(x=t,y=y),col='blue',size=1) 

p=p+theme_minimal(36)+xlab("Day")+ylab("Count of cases")
plot(p)



