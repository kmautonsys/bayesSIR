///////////////////////////////////////////////////////////////////////////////
// Copyright 2020 AutonSystems, LLC
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

functions {
  
  real bnd_(real x, real l, real u){
    return x>l ? (x<u?x:u) : l;
  }
  
  real SIR0(real y,real[] theta){
    return bnd_(-theta[1]*y+theta[2]*expm1(y)-theta[3], -50,0);
  }

  real[,] Forecast(real I0, real S0, data real[] ts, data int[] endpts, data int[] ts_is_output, real gamma, real[] beta, data int N){
    real out[N,2];
    real I1 = I0;
    real S1 = S0;
    real theta[3];
    int k = 1;
    int j = 1;
    real t0 = ts[1];
    real t = t0;
    real y = 0;
    real t1;
    real k1;
    real k2;
    real k3;
    real k4;
    if( ts_is_output[1] ){
      out[j,1] = I1;
      out[j,2] = S1;
      j = j+1;
    }
    theta[1] = gamma;
    theta[2] = beta[k];
    theta[3] = beta[k]*I1/S1;
    for(i in 2:size(ts)){
      t1 = ts[i];
      while(t<t1){
        real DT = 0.1;
        real dt = t1-t>DT?DT:t1-t;
        if(dt<=0) break; // just being paranoid about round-off error
        k1 = SIR0(y,theta);
        k2 = SIR0(y+0.5*dt*k1,theta);
        k3 = SIR0(y+0.5*dt*k2,theta);
        k4 = SIR0(y+dt*k3,theta);
        y=y+dt*(k1+2*k2+3*k3+k4)/6.0;
        t = t+dt;
      }
      if( ts_is_output[i] ){
        out[j,1] = -S1/beta[k]*SIR0(y,theta);
        out[j,2] = S1*exp(y);
        j=j+1;
        if(j>N) break;
      }
      if(t>=ts[endpts[k]]){ 
        I1 = -S1/beta[k]*SIR0(y,theta);
        S1 = S1*exp(y);
        k=k+1;
        theta[2] = beta[k];
        theta[3] = beta[k]*I1/S1;
        y = 0;
        t0 = t;
      }
    }
    return out;
  }

}

data {
  real<lower=1> TotalPop;
  
  // Outbreak data
  real I0MU; // mean initial infected population
  real<lower=0> I0SD; // standard deviation
  real S0MU; // mean initial recovered population
  real<lower=0> S0SD; // standard deviation
  
  int<lower=1> M; // the number of pieces in a piece-wise constant rate model
  real T0[M-1]; // knots, time points between pieces
  real<lower=0> tranMU[M]; // mean of each R transmision value (new infections/case)
  real<lower=0> tranSD[M]; // standard deviation
  real<lower=0> recMU; // mean recovery time
  real<lower=0> recSD; // standard deviation
  real HospMU; // expected proportion of cases requiring hospitalization
  real<lower=0> HospSD; // standard deviation
  real DeathMU; // expected proportion of cases resulting in death
  real<lower=0> DeathSD; // standard deviation
  real trendMU;
  real<lower=0> trendSD;
  
  // Forecasting data
  int<lower=0> N; // number of time steps to report
  real ts[N]; // forecast time-stamps
  
  // Historical data
  int<lower=0> Q; // number of observations
  real dat_ts[Q]; // time-stamps for historical data
  int<lower=0> dat_cases[Q]; // number of infections at times dat_ts
  int<lower=0,upper=1> dat_caseNA[Q]; // 1 = dat_cases[i] should be ignored.
  int<lower=0> dat_hospitalizations[Q];
  int<lower=0,upper=1> dat_hospNA[Q]; // 1 = dat_hospitalizations[i] should be ignored.
  int<lower=0> dat_deaths[Q];
  int<lower=0,upper=1> dat_deathNA[Q]; // 1 = dat_deaths[i] should be ignored.
  
  real phi; //negative binomial overdispersion
}

transformed data {

  real x_r[0]; // unused input to standard ODE signature
  int x_i[0]; // unused input to standard ODE signature
  
  real ts_[N+Q+M-1]; // time-stamps to get forward solution for SIR ODE
  int ts_is_output[N+Q+M-1]; // flag identifying if an endpoint is also a forecast time-stamp
  int endpts[M]; // endpoints of piece-wise periods
  
  int delta_cases[Q];
  int delta_hospitalizations[Q];
  int delta_deaths[Q];
  int dat_caseNA_[Q];
  int dat_hospNA_[Q];
  int dat_deathNA_[Q];

  // compute delta cases
  {// local variables
  int a=0;
  int b=0;
  int c=0;
  for(l in 1:Q){
    dat_caseNA_[l] = dat_caseNA[l];
    dat_hospNA_[l] = dat_hospNA[l];
    dat_deathNA_[l] = dat_deathNA[l];
    
    delta_cases[l] = dat_cases[l]-a;
    delta_hospitalizations[l] = dat_hospitalizations[l] - b;
    delta_deaths[l] = dat_deaths[l] - c;
    
    if(delta_cases[l]<0) dat_caseNA_[l] = 1;
    if(delta_hospitalizations[l]<0) dat_hospNA_[l] = 1;
    if(delta_deaths[l]<0) dat_deathNA_[l] = 1;
    
    a = dat_cases[l];
    b = dat_hospitalizations[l];
    c = dat_deaths[l];
  }
  }

  // Fill forward ODE time stamps
  if( N>0 ){// local variables
  int i;
  int j;
  int k;
  for(l in 1:(N+M-1)) ts_[l] = 0;
  for(l in 1:(N+M-1)) ts_is_output[l] = 0;
  for(l in 1:M) endpts[l] = 0;

  i=1;
  j=1;
  k=1;
  while(!(i>(N+Q) && j>M-1)){
    real T_ = (j<=M-1) ? (T0[j]<ts[N]?T0[j]:ts[N]) : ts[N];
    real t_ = (i<=Q) ? dat_ts[i] : (i<=(N+Q) ? ts[i-Q] : ts[N]);
    if(t_<T_){
      ts_[k] = t_;
      ts_is_output[k] = 1;
      i=i+1;
    } else {
      ts_[k] = T_;
      endpts[j] = k;
      if(T_==t_) i=i+1;
      if(T_==t_) ts_is_output[k] = 1;
      j=j+1;
    } 
    k=k+1;
  }
  }

}

parameters {

  real<lower=0,upper=log(TotalPop)> S0_; // S0_~Normal; S0=exp(S0_) => log_normal prior on S0
  real<lower=0,upper=log(TotalPop)> I0_; // I0_~Normal; I0=exp(I0_) => log_normal prior on I0
  real<lower=1,upper=4> rec_; // rec_~Normal; rec=exp(rec_) => log_normal prior on rec
  real<lower=-5,upper=5> tran_[M]; // tran_~Normal; tran=exp(tran_) => log_normal prior on tran

  real<lower=-10,upper=10> HospProp_; //logit parameter for proportion of cases resulting in death
  real<lower=-10,upper=10> DeathProp_; //logit parameter for proportion of cases resulting in death
}

transformed parameters {
  
  real S0 = exp(S0_);
  real I0 = exp(I0_);
  real gamma = 1.0/exp(rec_);
  real tran[M];
  real beta[M];
  
  real HospitalProportion = inv_logit(HospProp_);
  real DeathProportion = inv_logit(DeathProp_);
  
  for(l in 1:M) tran[l] = exp(tran_[l]);
  for(l in 1:M) beta[l] = tran[l]*gamma;

}

model {

  // simulate given initial conditions, and compare with observed data
  if( Q>0 ){
    real a = S0;
    real b = I0;
    
    real I_[Q,2] = Forecast(I0,S0,ts_,endpts,ts_is_output,gamma,beta,Q);
 
    for(i in 1:Q){
      real new_infections = a-I_[i,2]; // new infections
      real new_hosp = HospitalProportion * new_infections; // new infections
      real new_death = DeathProportion * (new_infections - (I_[i,1]-b)); // new resolved cases = new_infection - delta_infections
      
      if(new_infections<1e-50) new_infections=1e-50;
      if(new_hosp<1e-50) new_hosp=1e-50;
      if(new_death<1e-50) new_death=1e-50;
      
      if(phi<=0){
        
        if(dat_caseNA_[i]==0) target+=poisson_lpmf(delta_cases[i]|new_infections);
        if(dat_hospNA_[i]==0) target+=poisson_lpmf(delta_hospitalizations[i]|new_hosp);
        if(dat_deathNA_[i]==0) target+=poisson_lpmf(delta_deaths[i]|new_death);
        
      } else {
        
        if(dat_caseNA_[i]==0) target+=neg_binomial_2_lpmf(delta_cases[i]|new_infections,phi);
        if(dat_hospNA_[i]==0) target+=neg_binomial_2_lpmf(delta_hospitalizations[i]|new_hosp,phi);
        if(dat_deathNA_[i]==0) target+=neg_binomial_2_lpmf(delta_deaths[i]|new_death,phi);
        
      }
      a = I_[i,2];
      b = I_[i,1];
    }
  }
  
  if( M>1 ){
    for(i in 2:M){
      real d = beta[i]-beta[i-1];
      d ~ normal( trendMU, trendSD );
    }
  }
  
  // priors
  S0_ ~ normal(S0MU, S0SD);
  I0_ ~ normal(I0MU, I0SD);
  rec_ ~ normal(recMU, recSD);
  tran_ ~ normal(tranMU, tranSD);
  
  HospProp_ ~ normal(HospMU,HospSD);
  DeathProp_ ~ normal(DeathMU,DeathSD);
  
}

generated quantities {

  // forecast and back-cast 
  real Cases[Q+N];
  real Hosp[Q+N];
  real Deaths[Q+N];
  real HospLoad[Q+N];
  real I[Q+N];
  real S[Q+N];
  {
    real I_[N+Q,2] = Forecast(I0,S0,ts_,endpts,ts_is_output,gamma,beta,N+Q);

    for(i in 1:(N+Q)){
      Cases[i] = S0-I_[i,2];
      Hosp[i] = HospitalProportion*Cases[i];
      HospLoad[i] = HospitalProportion* I_[i,1];
      Deaths[i] = DeathProportion * (Cases[i] - (I_[i,1]-I0));
      I[i] = I_[i,1];
      S[i] = I_[i,2];
    }
  }

}
