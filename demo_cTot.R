library(ggplot2)
library(rstan)
library(zoo)
library(plyr)
rstan_options(auto_write=TRUE)

logit = function(x) -log(1/x-1)
inv_logit = function(x) 1/(1+exp(-x))

stanmodel = stan_model("bayesSIRcTot.stan")


# SOURCE: https://covidtracking.com/
library(jsonlite)
datdf <- fromJSON("https://covidtracking.com/api/states/daily")
state="NY"
datdf = datdf[datdf$state==state,]

ts = as.Date(as.character(datdf$date),"%Y%m%d")
ts = as.numeric(ts-(max(ts)+1))
dat_cases = as.integer(datdf$positive)
dat_hospitalizations = as.integer(datdf$hospitalized)
dat_deaths = as.integer(datdf$death)
idx = order(ts); ts=ts[idx]; dat_cases=dat_cases[idx]; dat_hospitalizations=dat_hospitalizations[idx]; dat_deaths=dat_deaths[idx];
dat_caseNA = is.na(dat_cases)
dat_hospNA = is.na(dat_hospitalizations)
dat_deathNA = is.na(dat_deaths)
dat_cases = na.locf(dat_cases,na.rm=FALSE) # forward fill NA
dat_cases[is.na(dat_cases)]=0 # zero out leading NAs
dat_hospitalizations = na.locf(dat_hospitalizations,na.rm=FALSE) # forward fill NA
dat_hospitalizations[is.na(dat_hospitalizations)]=0 # zero out leading NAs
dat_deaths = na.locf(dat_deaths,na.rm=FALSE) # forward fill NA
dat_deaths[is.na(dat_deaths)]=0 # zero out leading NAs
dat_cases = c(dat_cases[1],diff(dat_cases))
dat_hospitalizations = c(dat_hospitalizations[1],diff(dat_hospitalizations))
dat_deaths = c(dat_deaths[1],diff(dat_deaths))

# forecast time points
ts_ = as.numeric(1:60)

KnotsR = numeric(0)
K = 4
col_indices = lapply(1:K,function(k){ x=k+c(-1,0,1);x[x>=1 & x<=K]})
row_starts = cumsum(c(1,unlist(lapply(col_indices,length))))
col_indices = unlist(col_indices)

StatePop = list(AL=4908621,AK=734002,AZ=7378494,AR=3038999,CA=39937489,CO=5845526,CT=3563077,DE=982895,DC=720687,FL=21992985,GA=10736059,HI=1412687,ID=1826156,IL=12659682,IN=6745354,IA=3179849,KS=2910357,KY=4499692,LA=4645184,ME=1345790,MD=6083116,MA=6976597,MI=10045029,MN=5700671,MS=2989260,MO=6169270,MT=1086759,NE=1952570,NV=3139658,NH=1371246,NJ=8936574,NM=2096640,NY=19440469,NC=10611862,ND=761723,OH=11747694,OK=3954821,OR=4301089,PA=12820878,PR=3032165,RI=1056161,SC=5210095,SD=903027,TN=6897576,TX=29472295,UT=3282115,VT=628061,VA=8626207,WA=7797095,WV=1778070,WI=5851754,WY=567025)
rec_min = 1; rec_max = 21;
data = list(
  Knr=length(KnotsR),KnR=array(KnotsR,dim=length(KnotsR)),
  KnR_SD=0.05,
  K=K,Q=length(col_indices),P=length(row_starts),
  TotalPop=StatePop[[state]],
  col_indices=array(col_indices,dim=length(col_indices)),row_starts=row_starts,
  T=length(ts),ts=ts,ts_=ts_,S=length(ts_),
  dat_Cases=dat_cases,dat_caseNA=dat_caseNA,
  dat_Deaths=dat_deaths,dat_deathNA=dat_deathNA,
  rec_min = rec_min, rec_max = rec_max,
  gammaMU =           logit((1/7-1/rec_max)/(1/rec_min-1/rec_max)), #=> 14days
  gammaSD =           10.5,#0.025,
  within_MU =  log(2.5),#-log(10/2.5-1),
  within_SD =  20.5,
  between_MU = log(0.1),#-log(10/0.01-1),
  between_SD = 20.5,
  death_prop_MU =     -log(1/0.05-1),
  death_prop_SD =     0.1,
  I0MU =              array(rep(log(0.01),K),dim=K),
  I0SD =              array(rep(3,K),dim=K),
  phi = -1
)


#############
# MAP inference
init = list( gamma_=data$gammaMU, DeathProportion_=data$death_prop_MU,
             #I0_ = data$I0MU,
             S0_ = array(1/K,dim=K) )

#samp = stan("bayesSIRcTot.stan",data=data,iter=200,chains=1)

#opt = optimizing(stanmodel, data = data, hessian=TRUE, as_vector=FALSE,verbose=TRUE,algorithm="BFGS")

best = lapply(0:4,function(J){
  data0 = data
  if(J>0){
    tmpidx = (data$T-7*J):data$T
    data0$dat_caseNA[tmpidx]=TRUE; data0$dat_deathNA[tmpidx]=TRUE;
  }
  opt_=lapply(1:50,function(i){
    opt0=optimizing(stanmodel, data = data0, #init=init, hessian=FALSE, 
               as_vector=FALSE,verbose=TRUE,algorithm="BFGS")
    Sys.sleep(0.01)
    return(opt0)
  })
  opt = opt_[[which.max(unlist(lapply(opt_,function(o)o$value)))]]
  opt$jdx = -(J+1)
  return(opt)
})

#############
t = c(data$ts,data$ts_)
df = ldply(best,function(opt) data.frame(grp=opt$jdx,y=rowSums(opt$par$CIR[,,2]),t=t))
ggplot(df)+geom_line(aes(x=t,y=y,col=as.factor(grp)))+theme_minimal(12)+xlab("Day")+ylab("Count of cases") + theme(legend.position="none")
###
df = ldply(best,function(opt) data.frame(grp=opt$jdx,y=rowSums(opt$par$CIR[,,1]),t=t))
p=ggplot(df)
p=p+geom_point(data=data.frame(t=data$ts,y=cumsum(data$dat_Cases),grp=pmax(floor(data$ts/7),-5)),
                 aes(x=t,y=y,col=as.factor(grp)),shape = 21,  size = 2, stroke = 1)+
   # coord_cartesian(xlim=c(min(t),30),ylim=c(0,max(cumsum(data$dat_Cases))*1.25))
   scale_y_log10()
p=p+geom_line(aes(x=t,y=y,col=as.factor(grp)))+
  theme_minimal(18)+xlab("Day")+ylab("Count of cases")+
  theme(legend.position="none")
plot(p)
ggsave(paste0("Weekly_Forecasts_Cases_",state,"_K",K,"c.png"),p,width=6,height=6*3/5)
###
df = ldply(best,function(opt) data.frame(grp=opt$jdx,y=rowSums(opt$par$Deaths),t=t))
p=ggplot(df)
p=p+geom_point(data=data.frame(t=data$ts,y=cumsum(data$dat_Deaths),grp=pmax(floor(data$ts/7),-5)),
               aes(x=t,y=y,col=as.factor(grp)),shape = 21,  size = 2, stroke = 1)+
  coord_cartesian(xlim=c(min(t),30),ylim=c(0,max(cumsum(data$dat_Deaths))*1.25))
p=p+geom_line(aes(x=t,y=y,col=as.factor(grp)))+
  theme_minimal(18)+xlab("Day")+ylab("Count of deaths")+
  theme(legend.position="none")
plot(p)
ggsave(paste0("Weekly_Forecasts_Deaths_",state,"_K",K,".png"),p,width=6,height=6*3/5)



