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
  
  real[] SIR(real t,real[] y,real[] theta,real[] x_r,int[] x_i) {
    real dydt[1];
    real gamma = theta[1];
    real beta = theta[2];
    real C1 = theta[3];
    //-g*y+bP*exp(y)- [ b*C0+g*ln(P) ]
    dydt[1] = -gamma*y[1]+beta*exp(y[1])-C1;
    return dydt;
  }

  real[,] Forecast(real I0, real R0, data real t0, data real[] ts, data int[] endpts, data int[] ts_is_output, real gamma, real[] beta, data real TotalPop, data real[] x_r, data int[] x_i, data int M, data int N){
    real I[N,3];
    real TotalPop_;
    real R0_;
    real S0;
    real I0_;
    real C0;
    real theta[3];
    real y0[1];
    int j;
    
    theta[1] = gamma;
    I0_ = I0;
    R0_ = R0;
    S0 = TotalPop - I0 - R0;
    j = 1;
    for(i in 1:M){
      real t1;
      real beta_ = beta[i]*TotalPop/S0;
      int P = i>1 ? endpts[i]-endpts[i-1] : endpts[i];
      real S[P,1];
      TotalPop_ = TotalPop-R0_;
      C0 = I0_+S0-gamma/beta_*TotalPop_*log(S0/TotalPop_);
      t1 = ts[endpts[i]];
      theta[2] = beta_;
      theta[3] = beta_/TotalPop_*C0;
      y0[1] = log(S0/TotalPop_);
      if(I0_>0){
      if(i==1){
        S = integrate_ode_rk45(SIR, y0, t0, ts[1:endpts[i]], theta, x_r, x_i);
      } else {
        S = integrate_ode_rk45(SIR, y0, ts[endpts[i-1]], ts[(endpts[i-1]+1):endpts[i]], theta, x_r, x_i);
      }
      } else {
        for(k in 1:P) S[k,1] = y0[1];
      }
      for( k in 1:P ){
        int kk = i>1 ? endpts[i-1] : 0;
        if( ts_is_output[kk+k]==0) continue;
        I[j,2] = TotalPop_*exp(S[k,1]); // Susceptible population
        //g*ln(P*exp(y))/b-P*exp(y)+C0
        //g/b*y-P*exp(y)+C0+g/b*ln(P)
        I[j,1] = gamma/beta_*TotalPop_*S[k,1]-I[j,2]+C0; // Infected population
        if(I[j,1]<0) I[j,1]=0;
        I[j,3] = R0_+TotalPop_ - (I[j,2] + I[j,1]); // Recovered population
        if(j>=N) break;
        j = j + 1;
      }
      if(j>=N) break;
      S0 = TotalPop_*exp(S[P,1]);
      I0_ = gamma/beta_*TotalPop_*S[P,1]-S0+C0;
      R0_ = R0_+TotalPop_-S0-I0_;
      if(I0_<0) I0_=0;
    }
    return I;
  }

}

data {
  
  // Outbreak data
  real I0MU; // mean initial infected population
  real<lower=0> I0SD; // standard deviation
  real R0MU; // mean initial recovered population
  real<lower=0> R0SD; // standard deviation
  
  int<lower=1> M; // the number of pieces in a piece-wise constant rate model
  real T0[M-1]; // knots, time points between pieces
  real<lower=0> tranMU[M]; // mean of each R transmision value (new infections/case)
  real<lower=0> tranSD[M]; // standard deviation
  real<lower=0> recMU; // mean recovery time
  real<lower=0> recSD; // standard deviation
  real ConfirmMU; // expected proportion of cases confirmed
  real<lower=0> ConfirmSD; // standard deviation
  real DeathMU; // expected proportion of cases confirmed
  real<lower=0> DeathSD; // standard deviation
  real trendMU;
  real<lower=0> trendSD;
  
  // Population data
  real<lower=0> TotalPop; // total population size
  
  // Forecasting data
  int<lower=0> N; // number of time steps to report
  real ts[N]; // forecast time-stamps
  
  // Historical data
  int<lower=0> Q; // number of observations
  real dat_ts[Q]; // time-stamps for historical data
  int<lower=0> dat_cases[Q]; // number of infections at times dat_ts
  int<lower=0> dat_deaths[Q];
  int<lower=0,upper=1> dat_deathNA[Q]; // 1 = dat_deaths[i] should be ignored.
  
}

transformed data {

  real x_r[0]; // unused input to standard ODE signature
  int x_i[0]; // unused input to standard ODE signature
  
  real t0;
  real ts_[N+Q+M-1]; // time-stamps to get forward solution for SIR ODE
  int ts_is_output[N+Q+M-1]; // flag identifying if an endpoint is also a forecast time-stamp
  int endpts[M]; // endpoints of piece-wise periods
  
  int delta_cases[Q];
  int delta_deaths[Q];

  // compute delta cases
  {// local variables
  int a=0;
  int b=0;
  for(l in 1:Q){
    delta_cases[l] = dat_cases[l]-a;
    delta_deaths[l] = dat_deaths[l] - b;
    a = dat_cases[l];
    b = dat_deaths[l];
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
  t0 = ts_[1]-1;
  }

}

parameters {

  real ConfirmProportion; //logit parameter for proportion of cases confirmed
  real DeathProportion; //logit parameter for proportion of cases resulting in death
  real<lower=1> rec; // recovery time
  //real<lower=0,upper=TotalPop-2> I0; // initial infected population
  //real<lower=1,upper=TotalPop-I0-1> R0; // initial recovered population
  real<upper=log(TotalPop-2)> I0_;
  real<upper=log(TotalPop-1-exp(I0_))> R0_;
  
  // tran must be large enough for the initial conditions I0,R0,rec to represent a physically possible outbreak
  //real<lower=TotalPop/(I0+R0)*(log(TotalPop)-log(TotalPop-I0-R0))> tran[M]; // Reproduction value (number new infections/case)
  real<lower=0> tran[M];
}

transformed parameters {
  
  real I0 = exp(I0_);
  real R0 = exp(R0_);
  real DeathProportion_ = inv_logit(DeathProportion);
  real ConfirmProportion_ = inv_logit(ConfirmProportion);
  real gamma = 1/rec;
  real beta[M];
  for( l in 1:M ) beta[l] = tran[l]/rec;

}

model {

  // simulate given initial conditions, and compare with observed data
  if( Q>0 ){
    real a = TotalPop;
    real b = 0;
    real I_[Q,3] = Forecast(I0,R0,t0,ts_,endpts,ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,Q);
    for(i in 1:Q){
      real lambda = ConfirmProportion_ * (a-I_[i,2]);
      real lambda_death = DeathProportion_ * (I_[i,3]-b);
      if(lambda<1e-50) lambda=1e-50;
      if(lambda_death<1e-50) lambda_death=1e-50;
      if(dat_deathNA[i]==0) target+=poisson_lpmf(delta_deaths[i]|lambda_death);
      target+=poisson_lpmf(delta_cases[i]|lambda);
      a = I_[i,2];
      b = I_[i,3];
    }
  }
  
  if( M>1 ){
    for(i in 2:M){
      real d = tran[i]-tran[i-1];
      d ~ normal( trendMU, trendSD );
    }
  }
  
  // priors
  rec ~ normal(recMU, recSD);
  tran ~ normal(tranMU,tranSD);
  //for(i in 1:M) tran[i] ~ normal(tranMU[i], tranSD[i]) T[TotalPop/(I0+R0)*(log(TotalPop)-log(TotalPop-I0-R0)),];
  I0_ ~ normal(I0MU,I0SD);
  R0_ ~ normal(R0MU,R0SD) T[,log(TotalPop-1-I0)];
  ConfirmProportion ~ normal(ConfirmMU,ConfirmSD);
  DeathProportion ~ normal(DeathMU,DeathSD);
  
}

generated quantities {

  // forecast and back-cast 
  real FitCases[Q];
  real FitDeaths[Q];
  real FitI[Q];
  real FitS[Q];
  real FitR[Q];
  real Deaths[N];
  real I[N];
  real S[N];
  real R[N];
  {
    real I_[N+Q,3] = Forecast(I0,R0,t0,ts_,endpts,ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,N+Q);
    
    if( Q>0 ){
      real a = TotalPop;
      real b = 0;
      for(i in 1:Q){
        real lambda = ConfirmProportion_*(a-I_[i,2]);
        FitCases[i] = b + lambda;
        FitDeaths[i] = DeathProportion_ * I_[i,3];
        FitI[i] = I_[i,1];
        FitS[i] = I_[i,2];
        FitR[i] = I_[i,3];
        a = I_[i,2];
        b = FitCases[i];
      }
      
    }
    if( N>0 ){
      
      for(i in 1:N){
        Deaths[i] = DeathProportion_ * I_[Q+i,3];
        I[i] = I_[Q+i,1];
        S[i] = I_[Q+i,2];
        R[i] = I_[Q+i,3];
      }
      
    }
  }

}
