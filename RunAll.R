library(rstan)
library(plyr)
library(zoo)
rstan_options(auto_write=TRUE)

vpar = function(par) c(par$S0_,par$I0_,par$rec_,par$tran_,par$HospProp_,par$DeathProp_)
lpar = function(par,M) list(S0_=par[1],I0_=par[2],rec_=par[3],tran_=par[4:(4+M-1)],HospProp_=par[4+M],DeathProp_=par[5+M])
stanmodel = stan_model("bayesSIRv1.1.stan")

StatePop = list(AL=4908621,AK=734002,AZ=7378494,AR=3038999,CA=39937489,CO=5845526,CT=3563077,DE=982895,DC=720687,FL=21992985,GA=10736059,HI=1412687,ID=1826156,IL=12659682,IN=6745354,IA=3179849,KS=2910357,KY=4499692,LA=4645184,ME=1345790,MD=6083116,MA=6976597,MI=10045029,MN=5700671,MS=2989260,MO=6169270,MT=1086759,NE=1952570,NV=3139658,NH=1371246,NJ=8936574,NM=2096640,NY=19440469,NC=10611862,ND=761723,OH=11747694,OK=3954821,OR=4301089,PA=12820878,PR=3032165,RI=1056161,SC=5210095,SD=903027,TN=6897576,TX=29472295,UT=3282115,VT=628061,VA=8626207,WA=7797095,WV=1778070,WI=5851754,WY=567025)
# SOURCE: https://covidtracking.com/
library(jsonlite)
datdf_all <- fromJSON("https://covidtracking.com/api/states/daily")
states = unique(datdf_all$state)


ret_all = lapply(states,function(state){ 
ret_list = list()
for(LAG in c(0,3,7,14,21)){
  tryCatch({
datdf = datdf_all[datdf_all$state==state,]
dat_ts = as.Date(as.character(datdf$date),"%Y%m%d")
dat_ts = as.numeric(dat_ts-(max(dat_ts)+1))
dat_cases = as.integer(datdf$positive)
dat_hospitalizations = as.integer(datdf$hospitalized)
dat_deaths = as.integer(datdf$death)
idx = order(dat_ts); dat_ts=dat_ts[idx]; dat_cases=dat_cases[idx]; dat_hospitalizations=dat_hospitalizations[idx]; dat_deaths=dat_deaths[idx];
dat_caseNA = is.na(dat_cases)
dat_hospNA = is.na(dat_hospitalizations)
dat_deathNA = is.na(dat_deaths)
dat_cases = na.locf(dat_cases,na.rm=FALSE) # forward fill NA
dat_cases[is.na(dat_cases)]=0 # zero out leading NAs
dat_hospitalizations = na.locf(dat_hospitalizations,na.rm=FALSE) # forward fill NA
dat_hospitalizations[is.na(dat_hospitalizations)]=0 # zero out leading NAs
dat_deaths = na.locf(dat_deaths,na.rm=FALSE) # forward fill NA
dat_deaths[is.na(dat_deaths)]=0 # zero out leading NAs
Q = length(dat_ts)

if(LAG>0){
dat_caseNA[(Q-LAG+1):Q] = 1
dat_hospNA[(Q-LAG+1):Q] = 1
dat_deathNA[(Q-LAG+1):Q] = 1
}

# forecast time points
if(LAG==0){
  ts = as.numeric(1:60)
} else {
  ts = as.numeric(1)
}
N = length(ts)

# knots
T0 = seq(-14,-3,2)
T0 = T0[T0< -LAG-3]
M = length(T0)+1

data = list(
  M=M,
  T0=array(T0,dim=M-1),
  N=N,
  ts=array(ts,dim=N),
  Q=Q,
  dat_ts=array(dat_ts,dim=Q),
  dat_cases=array(dat_cases,dim=Q),
  dat_caseNA=array(dat_caseNA,dim=Q),
  dat_hospitalizations=array(dat_hospitalizations,dim=Q),
  dat_hospNA=array(dat_hospNA*1,dim=Q), # NA deaths should be forward-filled
  dat_deaths=array(dat_deaths,dim=Q),
  dat_deathNA=array(dat_deathNA*1,dim=Q), # NA deaths should be forward-filled
  TotalPop=StatePop[[state]],
  tranMU=array(log(2.5),dim=M),
  tranSD=array(1.0,dim=M),
  recMU=log(14),
  recSD=0.1,
  I0MU=0,
  I0SD=10,
  S0MU=log(0.5*StatePop[[state]]),
  S0SD=10,
  HospMU=-3,
  HospSD = 1,
  DeathMU = -4,
  DeathSD = 1,
  trendMU = 0,
  trendSD = 1.0,
  phi=-1
)
#############
# MAP inference
optim_ = lapply(1:50,function(i){ 
  tryCatch({
    opt = optimizing(stanmodel, data = data, hessian=TRUE, as_vector=FALSE,iter=100)
    Sys.sleep(0.01)
    mxi = which.max(opt$par$I)
    mxv = opt$par$I[mxi]
    opt$mxi = mxi
    opt$mxv = mxv
    print(paste0(state,"_",LAG,"  ",i,"  ",opt$value,"  ",opt$return_code))
    return(opt)
  },error=function(e){ print(e);return(NULL) },warning=function(w){print(w);return(NULL)})#,
})
idx = which(unlist(lapply(optim_,function(o) is.null(o))))
if(length(idx)>0) optim_=optim_[-idx]
# coverage:
optim2 = list(); optim=optim_
for(i in 1:length(optim)){
  p = unlist(lapply(1:length(optim),function(k) optim[[k]]$value))
  j = which.max(p)
  optim2[[length(optim2)+1]] = optim[[j]]
  opt = optim[[j]]
  par_ = vpar(optim[[j]]$par)
  # drop similar
  #drop_idx = which(unlist(lapply(optim,function(o) max(abs(vpar(o$par)-par_)/par_)<0.05 )))
  drop_idx = which(unlist(lapply(optim,function(o) abs(o$mxi-opt$mxi)<=2 &  abs(o$mxv-opt$mxv)<=0.01*opt$mxv )))
  optim[drop_idx] = NULL
  if(length(optim)==0) break
}
best = optim2

v=unlist(lapply(best,function(o) o$value)); mx=max(v);
best = best[which(unlist(lapply(best,function(o) o$value>mx+log(0.01) )))] # mx+log(tau) to keep samples with probability mass at least tau-% of the most likely point
v=unlist(lapply(best,function(o) o$value)); mx=max(v);
prob = exp(v-mx)/sum(exp(v-mx))
for(i in 1:length(best)) best[[i]]$Prob = prob[i]

ret_list[[length(ret_list)+1]] = list(LAG=LAG,modes=best,data=data,state=state)
},error=function(e){print(e);return(e)})
}
return(ret_list)
})

idx = which(5==unlist(lapply(ret_all,function(r) length(r))))
ret_all_ = ret_all[idx]
m =Reduce(rbind,lapply(ret_all_,function(r){
  dat_cases = r[[1]]$data$dat_cases
  Q = r[[1]]$data$Q
  lag = unlist(lapply(r,function(o) o$LAG))
  residual = unlist(lapply(r,function(o) {
    sum(unlist(lapply(o$modes,function(opt) opt$Prob*(dat_cases[Q]-opt$par$Cases[Q]))))
    }))
  error = matrix(abs(residual)/dat_cases[Q],nrow=1)
  colnames(error) = as.character(lag)
  return(error)
}))
colMeans(m)
qnorm(1-0.05/2)*sqrt(colMeans(m^2)-colMeans(m)^2)/sqrt(nrow(m))
unlist(lapply(1:ncol(m),function(i) median(m[,i])))

save(list="ret_all",file="RetAllStates_dev.RData")



st=0

st=st+1
print(ret_all[[st]][[1]]$state)
data = ret_all[[st]][[1]]$data
best = ret_all[[st]][[1]]$modes
#############
t = c(data$dat_ts,data$ts)
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$I,t=t))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
ggplot(df[df$t<=0,])+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Cases[1:data$Q],t=data$dat_ts))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=data$dat_ts,y=data$dat_cases),aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Hosp[1:data$Q],t=data$dat_ts))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=data$dat_ts,y=data$dat_hospitalizations)[data$dat_hospNA==0,],aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Deaths[1:data$Q],t=data$dat_ts))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=data$dat_ts,y=data$dat_deaths)[data$dat_deathNA==0,],aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
#############

df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$HospLoad,t=t))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(LogProb=opt$value,y=opt$par$Deaths,t=t))
df$LogProb = as.character( match(df$LogProb, sort(unique(df$LogProb),decreasing=TRUE)) )
ggplot(df)+geom_point(data=data.frame(t=data$dat_ts,y=data$dat_deaths)[data$dat_deathNA==0,],aes(x=t,y=y),size=2)+geom_line(aes(x=t,y=y,col=LogProb))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
#############
#############
#############
#############
probs0 = (5:49)/100 # posterior percentiles
probs = c(probs0,0.5,rev(1-probs0))
t = c(data$dat_ts,data$ts)
#############
posteriorI = NULL; FitCases=NULL
library(MASS)
for(k in 1:length(best)){
  opt=best[[k]]
  eig = eigen(-opt$hessian,symmetric=TRUE)
  eig$values[eig$values<0] = 0; eig$values[eig$values>0] = 1/eig$values[eig$values>0]
  Sigma = eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
  mu = vpar(opt$par)
  samps=lapply(1:1,function(i){
    s = mvrnorm(n=1,mu,Sigma)  
    s_ = lpar(s,data$M);
    if(s_$S0_<0) s_$S0_=0
    if(s_$S0_>log(data$TotalPop)) s_$S0_=log(data$TotalPop)
    if(s_$I0_<0) s_$I0_=0
    if(s_$I0_>s_$S0_) s_$I0_=s_$S0_
    return(s_)
  })
  fit <- stan(file = "bayesSIRv1.1.stan",data=data,init=samps,chains=length(samps),warmup=0,iter=1,cores=1,algorithm="Fixed_param",refresh=0)
  posteriorI_ = summary(fit, pars = "I", probs = probs)$summary
  posteriorI_ = cbind(posteriorI_,t); posteriorI_=cbind(posteriorI_,rep(k,nrow(posteriorI_)))
  posteriorI_ = cbind(posteriorI_,opt$par$I);
  colnames(posteriorI_)[(ncol(posteriorI_)-2):ncol(posteriorI_)] = c("t","idx","mode")
  if(is.null(posteriorI)){ posteriorI=posteriorI_
  }else{ posteriorI = rbind(posteriorI,posteriorI_)}
  FitCases_ = summary(fit, pars = "Cases", probs = probs)$summary
  FitCases_ = cbind(FitCases_,t); FitCases_=cbind(FitCases_,rep(k,nrow(FitCases_)))
  FitCases_ = cbind(FitCases_,opt$par$Cases);
  colnames(FitCases_)[(ncol(FitCases_)-2):ncol(FitCases_)] = c("t","idx","mode")
  if(is.null(FitCases)){ FitCases=FitCases_
  }else{ FitCases = rbind(FitCases,FitCases_)}
}
k#############
df = as.data.frame(posteriorI); row.names(df)=NULL; df$idx=as.factor(df$idx)
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",probs0*100,"l"),"c50",rev(paste0("c",probs0*100,"u")))
p = ggplot(df)
alpha = exp(seq(log(0.025),log(0.075),length.out=length(probs0)))
for( pr in probs0*100 ) p=p+geom_ribbon(aes_string(x="t",ymin=paste0("c",pr,"l"),ymax=paste0("c",pr,"u"),group="idx"),alpha=alpha[pr],fill="red")
p = p+geom_line(aes(x=t,y=c50,group=idx),size=1)
p = p+geom_line(aes(x=t,y=c25l,group=idx),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=c25u,group=idx),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=mode,group=idx),size=1,col='blue')
p=p+theme_minimal(36)+xlab("Day")+ylab("Infections")
plot(p)
###
df = as.data.frame(FitCases); row.names(df)=NULL; df$idx=as.factor(df$idx)
names(df)[names(df)%in%paste0(probs*100,"%")]=c(paste0("c",probs0*100,"l"),"c50",rev(paste0("c",probs0*100,"u")))
p = ggplot(df)
alpha = exp(seq(log(0.025),log(0.075),length.out=length(probs0)))
for( pr in probs0*100 ) p=p+geom_ribbon(aes_string(x="t",ymin=paste0("c",pr,"l"),ymax=paste0("c",pr,"u")),alpha=alpha[pr],fill="red")
p = p+geom_line(aes(x=t,y=c50),size=1)
p = p+geom_line(aes(x=t,y=c25l),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=c25u),size=0.5,lty=2)
p = p+geom_line(aes(x=t,y=mode,group=idx),size=1,col='blue')
p=p+geom_point(data=data.frame(t=data$dat_ts,y=data$dat_cases),aes(x=t,y=y),size=2)
p=p+theme_minimal(36)+xlab("Day")+ylab("Count of cases")
plot(p)


