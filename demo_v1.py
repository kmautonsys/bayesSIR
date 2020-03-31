#!/usr/bin/env python3
import pystan
import numpy as np
import json
import requests
import pandas as pd

r=requests.get("https://covidtracking.com/api/states/daily")
state_data = pd.DataFrame(r.json())
config = []
# outbreak priors

def load_config(file,data):
    config = dict()
    with open(file) as f:
        config = json.load(f)

    #Adjust the data types
    for key in config.keys():
        if isinstance(config[key],list):
            config[key] = np.array(config[key])

    # look up time series for the state
    state = config['state']
    state_data = data[data.state==state]
    
    config['dat_cases'] = np.flip(state_data.positive.values.astype('int'))
    
    config['dat_tests'] = np.flip(state_data.totalTestResults.values.astype('int'))
    #config['dat_hospitalized'] = state_data.hospitalized.values.astype('int')

    dat_ts = np.flip(np.array(state_data.date.values))
    
    config['dat_ts'] = dat_ts - np.max(dat_ts) -1

    config["TotalPop"] = 19440469 # this is new york specific
    config["nonCOVIDMU"] = 0.1*config['TotalPop']

    #fill in the derived parameters
    config["I0MU"] = max(config["dat_cases"]).astype('int') # mean current infected population
    config["R0MU"] = max(config["dat_cases"]).astype('int') # mean current recovered population

    config["M"]  = len(config["tranMU"])
    config["Q"]  = len(config["dat_ts"])
    config["T0"] = np.array([])

    # Forecast parameters
    ts = 1+np.arange(config['N'])
    config["ts"] = ts

    config.pop('state',None)
    return config


dat = load_config('ny_conf.json',state_data)
sm = pystan.StanModel(file="bayesSIRv1.1.stan")
fit = sm.sampling(data=dat, iter=1000, chains=4)

