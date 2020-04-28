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
    
    vector bnd_(vector x, real l, vector u){
      vector[rows(x)] out;
      for(i in 1:rows(x)) out[i] = x[i] > l ? (x[i]<u[i]?x[i]:u[i]) : l;
      return out;
    }
    
    matrix SIRc(real gamma, vector beta_values, data int[] col_indices, data int[] row_starts, matrix X){
      
      // X is Kx3, [S,I,R]
      int K = rows(X);
      matrix[K,4] out;
      vector[K] deltaS = X[,1].*csr_matrix_times_vector(K,K,beta_values,
                                                        col_indices,row_starts, X[,2]); // M*b, M is mxn, with values w
      deltaS = bnd_( deltaS, 0.0, 100.0*X[,1] ); // prevent over-shoot with extreme parameters
      out[,1] = -deltaS;
      out[,2] = deltaS - gamma*X[,2];
      out[,3] = gamma*X[,2];
      out[,4] = deltaS;
      return out;
      
    }
    
    real[,,] integrate_rk4(data real[] ts,matrix X0,data int[] col_indices,data int[] row_starts,matrix beta_values, real gamma, vector TotalPop,data real[] KnR){
      
      int Knr = size(KnR);
      int T = size(ts);
      int K = rows(X0);
      int r_idx = 1;
      real out[T,K,3];
      real t=ts[1];
      matrix[K,4] X;
      matrix[K,4] k1;
      matrix[K,4] k2;
      matrix[K,4] k3;
      matrix[K,4] k4;
      
      X[,1:3] = X0;    
      for(k in 1:K){
        out[1,k,1] = 0;
        out[1,k,2] = X[k,2];
        out[1,k,3] = X[k,3];
      }
      for(i in 2:size(ts)){
        real t1 = ts[i];
        for(k in 1:K) X[k,4] = 0;
        while(t<t1){
          real dt = t1-t>0.1?0.1:t1-t;
          if(dt<=0) break; // just being paranoid about round-off error
          
          if(r_idx<=Knr){
            if(t>=KnR[r_idx]) r_idx=r_idx+1;
          }
          
          k1 = SIRc(gamma,beta_values[,r_idx],col_indices,row_starts,X);
          k2 = SIRc(gamma,beta_values[,r_idx],col_indices,row_starts,X+0.5*dt*k1);
          k3 = SIRc(gamma,beta_values[,r_idx],col_indices,row_starts,X+0.5*dt*k2);
          k4 = SIRc(gamma,beta_values[,r_idx],col_indices,row_starts,X+dt*k3);
          X = X+dt*(k1+2*k2+3*k3+k4)/6.0;
          
          // check for overshoot
          X[,1] = bnd_(X[,1],0,TotalPop);
          X[,2] = bnd_(X[,2],0,TotalPop);
          //X[,3] = bnd_(X[,3],0,TotalPop);
          
          t = t+dt;
        }
        for(k in 1:K){
          out[i,k,1] = X[k,4];
          out[i,k,2] = X[k,2];
          out[i,k,3] = X[k,3];
        }
      }
      return out;
    }
    
  }

data {
  
  // Forecast time points
  int<lower=0> S; // number of additional time-stamps to predict
  real ts_[S]; // additional time points to predict
  
  // County data
  int<lower=1> K; // number of compartments
  real<lower=1> TotalPop; // total population
  
  // CSR data structure for county-to-county interactions
  int<lower=0> Q; // number of non-zero entries
  int<lower=0> P; // number of row-starts
  int<lower=1> col_indices[Q]; // CSR column indices
  int<lower=1> row_starts[P]; // CSR row starts
  
  // Outbreak data
  int<lower=0> Knr; // number of knots for transmission rate
  real KnR[Knr]; // knots for transmission rate (note, these will be rounded to the nearest 0.1)
  real<lower=0> KnR_SD; // trend filtering standard deviation
  real<lower=0> rec_min;
  real<lower=0> rec_max;
  
  int<lower=1> T; // number of time-stamps of measured data
  real ts[T]; // time-stamps of measured data
  
  int<lower=0> dat_Cases[T]; // new (incremental) cases at time t
  int<lower=0,upper=1> dat_caseNA[T]; // if 1, dat_Cases[t] should be ignored
  
  int<lower=0> dat_Deaths[T]; // new (incremental) deaths at time t
  int<lower=0,upper=1> dat_deathNA[T]; // if 1, dat_Deaths[t] should be ignored
  
  // priors
  real gammaMU;
  real<lower=0> gammaSD;
  real within_MU;
  real<lower=0> within_SD;
  real between_MU;
  real<lower=0> between_SD;
  real death_prop_MU;
  real<lower=0> death_prop_SD;
  vector[K] I0MU;
  vector<lower=0>[K] I0SD;
  
  real phi; //Negative binomial overdispersion. if phi<=0, then Poisson is used
}

transformed data {
  
  real all_ts[T+S];
  
  for(t in 1:T) all_ts[t] = ts[t];
  for(t in 1:S) all_ts[T+t] = ts_[t];
  all_ts = sort_asc(all_ts);
  
}

parameters {
  
  real gamma_;
  //matrix[Q,Knr+1] R_;
  matrix[K+1,Knr+1] R_;
  simplex[K] S0_;
  vector[K] I0_; // proportion of susceptible-population infected, initially
  
  real DeathProportion_;
}

transformed parameters {
  
  real gamma = inv_logit(gamma_)*(1.0/rec_min-1.0/rec_max)+1.0/rec_max;
  real DeathProportion = inv_logit(DeathProportion_);
  vector[K] S0 = S0_*TotalPop;
  vector[K] I0 = exp(I0_);
  
  matrix[Q,Knr+1] beta_values;
  {
    int k=1;
    for(i in 1:K){
      int J = row_starts[i+1]-row_starts[i]; // number of entries in this row
      for(j_ in 1:J){
        int j = col_indices[k];
        for(p_ in 1:(Knr+1)){
          beta_values[k,p_] = exp(R_[i==j?1:(i<j?i+1:j+1),p_]) ./ S0[i];
        }
        k = k+1;
      }
    }
  }
}

model {
  
  // observed counts
  real b = 0;
  matrix[K,3] X0;
  real CIR[T,K,3];
  
  X0[,1] = S0;
  X0[,2] = I0;
  for(k in 1:K) X0[k,3] = 0;
  CIR = integrate_rk4(ts,X0,col_indices,row_starts,beta_values,gamma,S0,KnR);
  
  for(t in 1:T){
      real new_cases = 1e-10 + sum(CIR[t,,1]);
      real new_deaths = 1e-10 + DeathProportion * (sum(CIR[t,,3])-b);
      b = sum(CIR[t,,3]);
      if( phi<=0 ){
        if(dat_caseNA[t]==0)  target+=poisson_lpmf(dat_Cases[t]|new_cases);
        if(dat_deathNA[t]==0) target+=poisson_lpmf(dat_Deaths[t]|new_deaths);
      } else {
        if(dat_caseNA[t]==0)  target+=neg_binomial_2_lpmf(dat_Cases[t]|new_cases, phi);
        if(dat_deathNA[t]==0) target+=neg_binomial_2_lpmf(dat_Deaths[t]|new_deaths, phi);
      }
  }
  
  // priors
  gamma_ ~ normal(gammaMU,gammaSD);
  DeathProportion_ ~ normal(death_prop_MU, death_prop_SD);
  
  if(Knr>0){
    for(i in 1:Knr){
      vector[Q] delta = R_[,i+1] - R_[,i];
      target+=normal_lpdf(delta|0,KnR_SD);
    }
  }
  
  I0_ ~ normal(I0MU,I0SD);
  {
    R_[1,] ~ normal(within_MU,within_SD);
    for(i in 1:(Knr+1)) R_[2:(K+1),i] ~ normal(between_MU,between_SD);
//    int k=1;
//    for(i in 1:K){
//      int J = row_starts[i+1]-row_starts[i]; // number of entries in this row
//      for(j_ in 1:J){
//        int j = col_indices[k];
//        if(i==j){ 
//          R_[k,] ~ normal(within_MU,within_SD);
//        } else {
//          R_[k,] ~ normal(between_MU,between_SD);
//        }
//        k = k+1;
//      }
//    }
  }
  
}

generated quantities {
  
  real CIR[T+S,K,3];
  real Deaths[T+S,K]; // cumulative deaths
  
  {
    vector[K] b;
    matrix[K,3] X0;
    
    for(k in 1:K) b[k] = 0;
    X0[,1] = S0;
    X0[,2] = I0;
    X0[,3] = b;
    CIR = integrate_rk4(all_ts,X0,col_indices,row_starts,beta_values,gamma,S0,KnR);
    
    for(t in 1:(T+S)){
      for(k in 1:K){
        real new_deaths = DeathProportion * (CIR[t,k,3]-b[k]); // new resolved cases
        Deaths[t,k] = t>1 ? Deaths[t-1,k]+new_deaths : new_deaths;
        if(t>1) CIR[t,k,1] = CIR[t-1,k,1]+CIR[t,k,1]; // return cumulative cases
        b[k] = CIR[t,k,3];
      }
    }
    
  }
}
