///////////////////////////////////////////////////////////////////////////////
// Copyright 2020 Auton Systems LLC 
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
    real hack = 0; 
    real gamma = theta[1];
    real beta = theta[2];
    real C0 = theta[3];
    if(y[1]>=0) hack = lmultiply(y[1],y[1]);
    dydt[1] = -gamma*hack+beta*y[1]*y[1]-beta*C0*y[1]; // lmultiply(y,y) = y*log(y)
    return dydt;
  }

  real[,] Forecast(real I0, real R0, data real[] ts, data int[] endpts, data int[] ts_is_output, real gamma, real[] beta, data real TotalPop, data real[] x_r, data int[] x_i, data int M, data int N, real time_reverse){
    real I[N,1];
    real I0_;
    real C0;
    real S0;
    real theta[3];
    real y0[1];
    int j;
    
    theta[1] = time_reverse*gamma;
    I0_ = I0;
    S0 = TotalPop - I0 - R0;
    j = 1;
    for(i in 1:M){
      real t1;
      if(endpts[i]<=0) continue;
      t1 = ts[endpts[i]];
      if(t1<=0) continue;
      {
      real beta_ = beta[time_reverse>0 ? i : M-i+1];
      int P = i>1 ? endpts[i]-endpts[i-1] : endpts[i];
      real S[P,1];
      C0 = I0_+S0-gamma/beta_*TotalPop*log(S0);
      theta[2] = time_reverse*beta_/TotalPop;
      theta[3] = C0;
      y0[1] = S0;
      if(i==1){
        S = integrate_ode_rk45(SIR, y0, 0, ts[1:endpts[i]], theta, x_r, x_i);
      } else if(ts[endpts[i-1]]<0) {
        S = integrate_ode_rk45(SIR, y0, 0, ts[(endpts[i-1]+1):endpts[i]], theta, x_r, x_i);
      } else {
        S = integrate_ode_rk45(SIR, y0, ts[endpts[i-1]], ts[(endpts[i-1]+1):endpts[i]], theta, x_r, x_i);
      }
      for( k in 1:P ){
        int kk = i>1 ? endpts[i-1] : 0;
        if( ts_is_output[kk+k]==0) continue;
        I[j,1] = gamma/beta_*TotalPop*log(S[k,1])-S[k,1]+C0;
        j = j + 1;
      }
      S0 = S[P,1];
      I0_ = gamma/beta_*TotalPop*log(S0)-S0+C0;
      }
    }
    if( time_reverse<0 ){ // flip time back
      int N2 = N/2;
      for( i in 1:N2 ){
        real v = I[i,1];
        I[i,1] = I[N-i+1,1];
        I[N-i+1,1] = v;
      }
    }
    return I;
  }

}

data {
  
  // Outbreak data
  real<lower=0> I0MU; // mean initial infected population
  real<lower=0> I0SD; // standard deviation
  real<lower=0> R0MU; // mean initial recovered population
  real<lower=0> R0SD; // standard deviation
  
  int<lower=1> M; // the number of pieces in a piece-wise constant rate model
  real T0[M-1]; // knots, time points between pieces
  real<lower=0> tranMU[M]; // mean of each R transmision value (new infections/case)
  real<lower=0> tranSD[M]; // standard deviation
  real<lower=0> recMU; // mean recovery time
  real<lower=0> recSD; // standard deviation
  
  // Population data
  real<lower=0> TotalPop; // total population size
  
  // Forecasting data
  int<lower=0> N; // number of time steps to report
  real<lower=0> ts[N]; // forecast time-stamps
  
  // Historical data
  int<lower=0> Q; // number of observations
  real<upper=0> dat_ts[Q]; // time-stamps for historical data
  real dat_tests[Q]; // unused
  real dat_cases[Q]; // number of infections at times dat_ts
  real<lower=0> datSD; // confidence in data
  
}

transformed data {

  real x_r[0]; // unused input to standard ODE signature
  int x_i[0]; // unused input to standard ODE signature
  
  real ts_[N+M-1]; // time-stamps to get forward solution for SIR ODE
  int ts_is_output[N+M-1]; // flag identifying if an endpoint is also a forecast time-stamp
  int endpts[M]; // endpoints of piece-wise periods
  
  real dat_ts_[Q+M-1]; // time-stamps to get backward solution for SIR ODE
  int dat_ts_is_output[Q+M-1]; // flag identifying if an endpoint is also a data time-stamp
  int dat_endpts[M]; // endpoints of piece-wise periods
  
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
  while(!(i>N && j>M-1)){
    real t_ = (i<=N) ? ts[i] : ts[N];
    real T_ = (j<=M-1) ? (T0[j]<ts[N]?T0[j]:ts[N]) : ts[N];
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
  
  // Fill backward ODE time stamps
  if( Q>0 ){// local variables
  int i;
  int j;
  int k;
  real dat_ts_rev[Q];
  real T0_rev[M-1];
  for(l in 1:Q) dat_ts_rev[l] = -dat_ts[Q-l+1];
  for(l in 1:(M-1)) T0_rev[l] = -T0[M-l];
  for(l in 1:(Q+M-1)) dat_ts_[l] = 0;
  for(l in 1:(Q+M-1)) dat_ts_is_output[l] = 0;
  for(l in 1:M) dat_endpts[l] = 0;

  i=1;
  j=1;
  k=1;
  while(!(i>Q && j>M-1)){
    real t_ = (i<=Q) ? dat_ts_rev[i] : dat_ts_rev[Q];
    real T_ = (j<=M-1) ? (T0_rev[j]<dat_ts_rev[Q]?T0_rev[j]:dat_ts_rev[Q]) : dat_ts_rev[Q];
    if(t_<T_){
      dat_ts_[k] = t_;
      dat_ts_is_output[k] = 1;
      i=i+1;
    } else {
      dat_ts_[k] = T_;
      dat_endpts[j] = k;
      if(T_==t_) i=i+1;
      if(T_==t_) dat_ts_is_output[k] = 1;
      j=j+1;
    } 
    k=k+1;
  }
  }
  
}

parameters {

  real<lower=1> rec; // recovery time
  real<lower=0> tran[M]; // Reproduction value (number new infections/case)
  real<lower=0,upper=TotalPop/2> I0; // initial infected population
  real<lower=0,upper=TotalPop/2> R0; // initial recovered population

}

transformed parameters {
  
  real gamma = 1/rec;
  real beta[M];
  for( l in 1:M ) beta[l] = tran[l]/rec;

}

model {
  
  // back-cast initial conditions, and compare with observed data
  real dat_I[Q,1];
  if( Q>0 ){
    
    dat_I = Forecast(I0,R0,dat_ts_,dat_endpts,dat_ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,Q,-1);
    
    for(i in 1:Q) dat_cases[i] ~ normal(dat_I[i,1],datSD);
  }
  
  // priors
  rec ~ normal(recMU, recSD);
  tran ~ normal(tranMU, tranSD);
  I0 ~ normal(I0MU,I0SD);
  R0 ~ normal(R0MU,R0SD);
  
}

generated quantities {

  // forecast and back-cast 
  real dat_I[Q,1];
  real I[N,1];
  {
    if( Q>0 ){
      
      dat_I = Forecast(I0,R0,dat_ts_,dat_endpts,dat_ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,Q,-1);
      
    }
    if( N>0 ){
      
      I = Forecast(I0,R0,ts_,endpts,ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,N,1);
      
    }
  }

}
