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

  real ConfirmationProb( real TestingPop, real[] I, real[] lastI ) {
    //return I[1]/(TestingPop + I[1]);
    real new_cases = -(I[2]-lastI[2]);
    return new_cases/(TestingPop+new_cases);
  }

  real[] SIR(real t,real[] y,real[] theta,real[] x_r,int[] x_i) {
    real dydt[1];
    real gamma = theta[1];
    real beta = theta[2];
    dydt[1] = -gamma*y[1]+beta*exp(y[1])-beta;
    return dydt;
  }

  real[,] Forecast(real I0, real R0, data real[] ts, data int[] endpts, data int[] ts_is_output, real gamma, real[] beta, data real TotalPop, data real[] x_r, data int[] x_i, data int M, data int N, real time_reverse){
    real I[N,3];
    real S0;
    real theta[2];
    real y0[1];
    int j;
    
    theta[1] = time_reverse*gamma;
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
      theta[2] = time_reverse*beta_;
      y0[1] = log(S0/TotalPop);
      if(i==1){
        S = integrate_ode_bdf(SIR, y0, 0, ts[1:endpts[i]], theta, x_r, x_i);
      } else if(ts[endpts[i-1]]<0) {
        S = integrate_ode_bdf(SIR, y0, 0, ts[(endpts[i-1]+1):endpts[i]], theta, x_r, x_i);
      } else {
        S = integrate_ode_bdf(SIR, y0, ts[endpts[i-1]], ts[(endpts[i-1]+1):endpts[i]], theta, x_r, x_i);
      }
      for( k in 1:P ){
        int kk = i>1 ? endpts[i-1] : 0;
        if( ts_is_output[kk+k]==0) continue;
        I[j,2] = TotalPop*exp(S[k,1]); // Susceptible population
        I[j,1] = gamma/beta_*TotalPop*S[k,1]-I[j,2]+TotalPop; // Infected population
        I[j,3] = TotalPop - (I[j,2] + I[j,1]); // Recovered population
        j = j + 1;
      }
      S0 = TotalPop*exp(S[P,1]);
      }
    }
    if( time_reverse<0 ){ // flip time back
      int N2 = N/2;
      for( i in 1:N2 ){
        real v = I[i,1];
        real w = I[i,2];
        real z = I[i,3];
        I[i,1] = I[N-i+1,1];
        I[i,2] = I[N-i+1,2];
        I[i,3] = I[N-i+1,3];
        I[N-i+1,1] = v;
        I[N-i+1,2] = w;
        I[N-i+1,3] = z;
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
  int<lower=0> future_tests[N]; // anticipated number of new tests in future (not cumulative)
  
  // Historical data
  int<lower=0> Q; // number of observations
  real<upper=0> dat_ts[Q]; // time-stamps for historical data
  int dat_tests[Q]; // total number of tests performed, to-date
  int dat_cases[Q]; // total number of confirmed infections, to-date

  real ConfirmMU; // mean proportion of cases which are severe enough to garner testing
  real ConfirmSD; // standard deviation
  real nonCOVIDMU; // mean proportion of cases which are severe enough to garner testing
  real nonCOVIDSD; // standard deviation

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
  
  int delta_cases[Q];
  int delta_tests[Q];

  // compute delta cases & tests
  {// local variables
  int a=0;
  int b=0;
  for(l in 1:Q){
    delta_cases[l] = dat_cases[l]-a;
    delta_tests[l] = dat_tests[l]-b;
    a = dat_cases[l];
    b = dat_tests[l];
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
  real<lower=0,upper=TotalPop-2> I0; // initial infected population
  real<lower=1,upper=TotalPop-I0-1> R0; // initial recovered population
  
  // tran must be large enough for the initial conditions I0,R0,rec to represent a physically possible outbreak
  real<lower=TotalPop/(I0+R0)*(log(TotalPop)-log(TotalPop-I0-R0))> tran[M]; // Reproduction value (number new infections/case)
  
  real ConfirmProportion;
  real TestingPop;
}

transformed parameters {
  
  real TestingPop_ = exp(TestingPop);
  real ConfirmProportion_ = inv_logit(ConfirmProportion);
  real gamma = 1/rec;
  real beta[M];
  for( l in 1:M ) beta[l] = tran[l]/rec;

}

model {
  
  // back-cast initial conditions, and compare with observed data
  if( Q>0 ){
    real dat_I[Q,3] = Forecast(I0,R0,dat_ts_,dat_endpts,dat_ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,Q,-1);
    real lastI[3]; lastI[1]=0; lastI[2]=TotalPop; lastI[3]=0;
    
    for(i in 1:Q){
      real new_cases = -(dat_I[i,2]-lastI[2]);
      real theta = ConfirmationProb(TestingPop_,dat_I[i,],lastI);
      delta_cases[i] ~ binomial( delta_tests[i], theta );
      delta_tests[i] ~ poisson( ConfirmProportion_*new_cases );
      lastI = dat_I[i,];
      
    }
  }
  
  // priors
  rec ~ normal(recMU, recSD);
  for(i in 1:M) tran[i] ~ normal(tranMU[i], tranSD[i]) T[TotalPop/(I0+R0)*(log(TotalPop)-log(TotalPop-I0-R0)),];
  I0 ~ normal(I0MU,I0SD);
  R0 ~ normal(R0MU,R0SD) T[1,TotalPop-I0-1];
  TestingPop ~ normal(nonCOVIDMU,nonCOVIDSD);
  ConfirmProportion ~ normal(ConfirmMU,ConfirmSD);

}

generated quantities {

  // forecast and back-cast 
  real FitI[Q];
  real FitR[Q];
  real FitS[Q];
  real FitCases[Q];
  real I[N];
  real R[N];
  real S[N];
  real Cases[N];
  {
    if( Q>0 ){
      real c=0;
      real dat_I_[Q,3] = Forecast(I0,R0,dat_ts_,dat_endpts,dat_ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,Q,-1);
      real lastI[3]; lastI[1]=0; lastI[2]=TotalPop; lastI[3]=0;
      
      for(i in 1:Q){
        real theta = ConfirmationProb(TestingPop_,dat_I_[i,],lastI);
        FitI[i] = dat_I_[i,1];
        FitR[i] = dat_I_[i,3];
        FitS[i] = dat_I_[i,2];
        FitCases[i] = c + delta_tests[i]*theta;
        c = FitCases[i];
        lastI=dat_I_[i,];
      }
    }
    if( N>0 ){
      real c=0;
      real I_[N,3] = Forecast(I0,R0,ts_,endpts,ts_is_output,gamma,beta,TotalPop,x_r,x_i,M,N,1);
      real lastI[3]; lastI[1]=0; lastI[2]=TotalPop; lastI[3]=0;
      
      for(i in 1:N){
        real theta = ConfirmationProb(TestingPop_,I_[i,],lastI);
        I[i] = I_[i,1];
        R[i] = I_[i,3];
        S[i] = I_[i,2];
        Cases[i] = c + future_tests[i]*theta;
        c = Cases[i];
        lastI = I_[i,];
      }
    }
  }

}
