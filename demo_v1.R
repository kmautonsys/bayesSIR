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
#TotalPop = 10045029 #MI
#TotalPop = 39937489 #CA
#TotalPop = 12820878 #PA
#TotalPop = 11747694 #OH
Q = length(dat_ts)

# outbreak priors
tranMU = c(2.5) # mean Reproduction number (number new infections/case)
tranSD = c(1.4) # standard deviation
T0 = numeric(0)
M = length(tranMU)

tranMU = c(2.5,2.5,2.5) # mean Reproduction number (number new infections/case)
tranSD = c(1.4,1.4,1.4) # standard deviation
T0 = c(-14,-7)
M = length(tranMU)

T0 = seq(-14,-3,2)
tranMU = rep(2.5,length(T0)+1)
tranSD = rep(1.4,length(T0)+1)
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
  recSD=0.1,
  I0MU=1,
  I0SD=10,
  R0MU=1,
  R0SD=10,
  N=N,
  ts=array(ts,dim=N),
  Q=Q,
  dat_ts=array(dat_ts,dim=Q),
  dat_cases=array(dat_cases,dim=Q),
  dat_hospitalizations=array(dat_hospitalizations,dim=Q),
  dat_hospNA=array(dat_hospNA*1,dim=Q), # NA deaths should be forward-filled
  dat_deaths=array(dat_deaths,dim=Q),
  dat_deathNA=array(dat_deathNA*1,dim=Q), # NA deaths should be forward-filled
  TotalPop=TotalPop,
  ConfirmMU = -4,
  ConfirmSD = 0.5,
  DeathMU = -5,
  DeathSD = 1.5,
  HospMU=-5,
  HospSD = 1.5,
  trendMU = 0,
  trendSD = 0.05,
  fix_confirm = 0
)

library(plyr)
#############
# MAP inference
vpar = function(par) c(par$ConfirmProportion,par$HospitalProportion,par$DeathProportion,par$rec,par$I0_,par$R0_,par$tran)
lpar = function(par) list(ConfirmProportion=par[1],HospitalProportion=par[2],DeathProportion=par[3],rec=par[4],I0_=par[5],R0_=par[6],tran=par[7:(7+M-1)])

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

opt = best[[1]]; p=opt$par$S+opt$par$I+opt$par$R
c(min(p-TotalPop),max(p-TotalPop))
x=opt$par$S;c(min(x),max(x))/TotalPop
x=opt$par$I;c(min(x),max(x))/TotalPop
x=opt$par$R;c(min(x),max(x))/TotalPop

print(best[[1]]$par$tran)
print(best[[1]]$par$rec)
print(best[[1]]$par$ConfirmProportion_)
print(best[[1]]$par$HospitalProportion_)
print(best[[1]]$par$DeathProportion_)
#############
t = c(dat_ts,ts)
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$I,t=t))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Cases[1:Q],t=dat_ts))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=dat_ts,y=dat_cases),aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Hosp[1:Q],t=dat_ts))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=dat_ts,y=dat_hospitalizations)[dat_hospNA==0,],aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Deaths[1:Q],t=dat_ts))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=dat_ts,y=dat_deaths)[dat_deathNA==0,],aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
#############
opt = best[[1]]


#############
# THIS -OR- THAT
library(MASS)
eig = eigen(-opt$hessian,symmetric=TRUE)
eig$values[eig$values<0] = 0; eig$values[eig$values>0] = 1/eig$values[eig$values>0]
Sigma = eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
mu = vpar(opt$par)
samps=lapply(1:100,function(i){
    c = 10;
    s = mvrnorm(n=1,mu,Sigma*c)  
    logprob = -0.5*t(s)%*%Sigma%*%(s)*c
    s = lpar(s); s$logprob=logprob; return(s)
})
fit <- stan(file = "bayesSIRv1.1.stan",data=data,init=samps,chains=length(samps),warmup=0,iter=1,cores=2,algorithm="Fixed_param",refresh=0)
#############
# THAT
# Can be SLOW
fit <- stan(file = "bayesSIRv1.1.stan",data=data,init=list(opt$par),chains=1,warmup=1000,iter=2000,cores=1)
#############




probs0 = (5:49)/100 # posterior percentiles
probs = c(probs0,0.5,rev(1-probs0))

FitCases = summary(fit, pars = "Cases", probs = probs)$summary
FitCases = cbind(FitCases,t)  
colnames(FitCases)[ncol(FitCases)] = "t"

posteriorI = summary(fit, pars = "I", probs = probs)$summary
posteriorI = cbind(posteriorI,t)  
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

p = p+geom_line(data=data.frame(y=opt$par$I/TotalPop*100,t=t),aes(x=t,y=y),col='blue',size=1) 

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

p = p+geom_line(data=data.frame(y=opt$par$Cases,t=t),aes(x=t,y=y),col='blue',size=1) 
#p=p+coord_cartesian(ylim=c(0,max(dat_cases)*1.1),xlim=c(min(dat_ts),0))
p=p+theme_minimal(36)+xlab("Day")+ylab("Count of cases")
plot(p)

# Plot fit to data
FitCases = summary(fit, pars = "Cases", probs = probs)$summary
FitCases = cbind(FitCases,t)  
colnames(FitCases)[ncol(FitCases)] = "t"

df = as.data.frame(FitCases)
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",probs0*100,"l"),"c50",rev(paste0("c",probs0*100,"u")))
p = ggplot(df)
alpha = exp(seq(log(0.025),log(0.075),length.out=length(probs0)))
for( pr in probs0*100 ) p=p+geom_ribbon(aes_string(x="t",ymin=paste0("c",pr,"l"),ymax=paste0("c",pr,"u")),alpha=alpha[pr],fill="red")
p = p+geom_line(aes(x=t,y=c50),size=1)
p = p+geom_line(aes(x=t,y=c25l),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=c25u),size=0.5,lty=2)
p=p+geom_point(data=data.frame(t=dat_ts,y=dat_cases),aes(x=t,y=y),size=2)

p = p+geom_line(data=data.frame(y=opt$par$Cases,t=t),aes(x=t,y=y),col='blue',size=1) 
#p=p+coord_cartesian(ylim=c(0,max(dat_cases)*1.1),xlim=c(min(dat_ts),0))
p=p+theme_minimal(36)+xlab("Day")+ylab("Count of cases")
plot(p)


