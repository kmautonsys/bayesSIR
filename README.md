Copyright 2020 AutonSystems, LLC

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# bayesSIR

Bayesian data fusion and forecasting for contagious disease.

# Prerequisites

- Rstudio 3.4.4 or higher
  - Including packages: rstan, jsonlite, ggplot2, and openxlsx
- rstan version 2.19.3 or higher

# Summary

This is a data driven Bayesian forecasting model for contagious disease. It is written in stan, with a front end in R. It can use available data to estimate outbreak parameters, forecast future infections, and give uncertainty estimates of current and future infections.

The basic parameters are:
- Infected Population
- Recovered Population
- Number of new infections/case

Each of these basic parameters is specified with a mean and standard deviation and can be defined in a piece wise linear manner to model the impacts of different scenarios and interventions.

# Data
The expected input data for the stan script is:

  // Outbreak data
  real I0MU; // mean initial infected population
  real I0SD; // standard deviation
  real R0MU; // mean initial recovered population
  real R0SD; // standard deviation
  
  int<lower=1> M; // the number of pieces in a piece-wise constant rate model
  real T0[M-1]; // knots, time points between pieces
  real tranMU[M]; // mean of each R transmision value (new infections/case)
  real tranSD[M]; // standard deviation
  real recMU; // mean recovery time
  real recSD; // standard deviation
  
  // Population data
  real TotalPop; // total population size
  
  // Forecasting data
  int N; // number of time steps to report
  real ts[N]; // forecast time-stamps
  
  // Historical data
  int Q; // number of observations
  rea dat_ts[Q]; // time-stamps for historical data
  real dat_tests[Q]; // unused
  real dat_cases[Q]; // number of infections at times dat_ts
  real datSD; // confidence in data
  
# Running
The bayesSIR model is driven by a front end that feeds it data and processes the results. We have provided a demonstration front end that uses synthetic data and produces plots as output. To run our demo script call "source(demo.R)" from a R session or "Rscript demo.R" from the command line. The scripts can also be invoked from RStudio with line-by-line execution or using the source buttion.




