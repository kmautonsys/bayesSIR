#!/usr/bin/env python3
import pystan
import numpy as np
import json
import requests
import pandas as pd
import matplotlib.pyplot as plt

config = []

r=requests.get("https://covidtracking.com/api/states/daily")
state_data = pd.DataFrame(r.json())


def load_config(file,data):
    """
    This loads a config file for a stan run. Config file is assumed to be in JSON format.
    The data file is assumed to have positive, deaths, and hospitalized features. Postive refers to positive tests, a proxy for confirmed cases.
    This need some error checking, particularlly for the case that a required parameter is missing.
    """
    config = dict()
    
    with open(file) as f:
        config = json.load(f)
    #Adjust the data types
    for key in config.keys():
        if isinstance(config[key],list):
            config[key] = np.array(config[key])

    with open('./state_populations.json') as f:
        state_pops = json.load(f)
    for key in state_pops.keys():
        state_pops[key] = int(state_pops[key])

    # look up time series for the state
    state = config['state']
    state_data = data[data.state==state]
    state_pop = state_pops[state]

    # These parameters can be loaded from the data
    config['dat_cases'] = np.flip(state_data.positive.values.astype('int'))
    config['dat_caseNA'] = np.isnan(state_data.positive.values).astype('int')

    config['dat_deaths'] = np.flip(np.nan_to_num(state_data.death.values).astype('int'))
    config['dat_deathNA'] = np.isnan(state_data.death.values).astype('int') # NA deaths should be forward-filled

    config['dat_hospitalized'] = np.flip(state_data.hospitalized.values.astype('int'))
    config['dat_hospNA'] = np.isnan(state_data.hospitalized.values).astype('int') #NA hospitalizations should be forward filled

    dat_ts = np.flip(np.array(state_data.date.values))

    # fill in the derived parameters and adjust data types
    config['dat_ts'] = dat_ts - np.max(dat_ts) -1
    config['TotalPop'] = float(state_pop)
#    config["nonCOVIDMU"] = 0.1*config['TotalPop']


    config['M']  = len(config["tranMU"])
    config['Q']  = len(config["dat_ts"])
    #config['recMU']=np.log(config['recMU'])

    config['S0MU'] = np.log(0.5*state_pop)

    # Forecast parameters
    ts = 1+np.arange(config['N'])
    config["ts"] = ts
    config.pop('state',None)
    #print(config)
    return config


def describe_fit(fit):
    """
        Simple analysis of the fitted model and display of the population evaluation
    """
    summary = fit.summary()
    summary = pd.DataFrame(summary['summary'], index = summary['summary_rownames'], columns = summary['summary_colnames'])
    print("Estimated transmited rate: {:.2f} ({:.2f})".format(summary['mean']['tran[1]'], summary['sd']['tran[1]']))
    print("Recovery time: {:.2f} ({:.2f}) days".format(summary['mean'].rec, summary['sd'].rec))
    print("Proportion population affected: {:.2f}% ({:.2f})".format(summary['mean'].ConfirmProportion_* 100, summary['sd'].ConfirmProportion_ * 100))
    # Evolution Cases
    for (group, color, id) in [('Infected', 'orange', 'I'), ('Dead', 'red', 'Deaths')]:
        training, prediction = summary.loc[[i for i in summary.index if 'Fit{}['.format(id) in i]], summary.loc[[i for i in summary.index if ('{}['.format(id) in i) and ('Fit' not in i)]]
        plt.plot(np.arange(-len(training),0) + 1, training['mean'], alpha = 0.5, color = color)
        plt.fill_between(np.arange(-len(training),0) + 1, training['2.5%'], training['97.5%'], alpha = 0.25, color = color)
        plt.plot(np.arange(len(prediction) + 1), [training['mean'].iloc[-1]] + prediction['mean'].tolist(), label = group, ls = ':', color = color)
        plt.fill_between(np.arange(len(prediction) + 1), [training['2.5%'].iloc[-1]] + prediction['2.5%'].tolist(), [training['97.5%'][-1]] + prediction['97.5%'].tolist(), alpha = 0.5, color = color)
    plt.xlabel('Days')
    plt.ylabel('Number')
    plt.legend()
    plt.show()
dat = load_config('ny_conf.json',state_data)
#print(dat)
sm = pystan.StanModel(file="bayesSIRv1.1.stan")
#fit = sm.sampling(data=dat, iter=1000, chains=4)
fit = sm.optimizing(data=dat)
#describe_fit(fit)


# #!/usr/bin/env python3
# import pystan
# import numpy as np
# import json
# import requests
# import pandas as pd
# import matplotlib.pyplot as plt

# r=requests.get("https://covidtracking.com/api/states/daily")
# state_data = pd.DataFrame(r.json())
# config = []

# # outbreak priors
# def load_config(file,data):
#     config = dict()
#     with open(file) as f:
#         config = json.load(f)
#     #Adjust the data types
#     for key in config.keys():
#         if isinstance(config[key],list):
#             config[key] = np.array(config[key])
#     # look up time series for the state
#     state = config['state']
#     state_data = data[data.state==state]
#     config['dat_cases'] = np.flip(state_data.positive.values.astype('int'))
#     config['dat_tests'] = np.flip(state_data.totalTestResults.values.astype('int'))
#     config['dat_deaths'] = np.flip(np.nan_to_num(state_data.death.values).astype('int'))
#     config['dat_deathNA'] = np.isnan(state_data.death.values).astype('int')
#     #config['dat_hospitalized'] = state_data.hospitalized.values.astype('int')
#     dat_ts = np.flip(np.array(state_data.date.values))
#     config['dat_ts'] = dat_ts - np.max(dat_ts) -1
#     config["TotalPop"] = 19440469 # this is new york specific
#     config["nonCOVIDMU"] = 0.1*config['TotalPop']
#     #fill in the derived parameters
#     config["I0MU"] = max(config["dat_cases"]).astype('int') # mean current infected population
#     config["R0MU"] = max(config["dat_cases"]).astype('int') # mean current recovered population
#     config["M"]  = len(config["tranMU"])
#     config["Q"]  = len(config["dat_ts"])
#     config["DeathMU"] = 0
#     config["DeathSD"] = 2
#     # Forecast parameters
#     ts = 1+np.arange(config['N'])
#     config["ts"] = ts
#     config.pop('state',None)
#     return config

# def describe_fit(fit):
#     """
#         Simple analysis of the fitted model and display of the population evaluation
#     """
#     summary = fit.summary()
#     summary = pd.DataFrame(summary['summary'], index = summary['summary_rownames'], columns = summary['summary_colnames'])
#     print("Estimated transmited rate: {:.2f} ({:.2f})".format(summary['mean']['tran[1]'], summary['sd']['tran[1]']))
#     print("Recovery time: {:.2f} ({:.2f}) days".format(summary['mean'].rec, summary['sd'].rec))
#     print("Proportion population affected: {:.2f}% ({:.2f})".format(summary['mean'].ConfirmProportion_* 100, summary['sd'].ConfirmProportion_ * 100))
#     # Evolution Cases
    
#     for (group, color, id) in [('Infected', 'orange', 'I'), ('Dead', 'red', 'Deaths')]:
#         training, prediction = summary.loc[[i for i in summary.index if 'Fit{}['.format(id) in i]], summary.loc[[i for i in summary.index if ('{}['.format(id) in i) and ('Fit' not in i)]]
#         plt.plot(np.arange(-1*len(training),0) + 1, training['mean'], alpha = 0.5, color = color)
#         plt.fill_between(np.arange(-len(training),0) + 1, training['2.5%'], training['97.5%'], alpha = 0.25, color = color)
#         plt.plot(np.arange(len(prediction) + 1), [training['mean'][-1]] + prediction['mean'].tolist(), label = group, ls = ':', color = color)
#         plt.fill_between(np.arange(len(prediction) + 1), [training['2.5%'][-1]] + prediction['2.5%'].tolist(), [training['97.5%'][-1]] + prediction['97.5%'].tolist(), alpha = 0.5, color = color)
#     plt.xlabel('Days')
#     plt.ylabel('Number')
#     plt.legend()
#     plt.show()


# dat = load_config('ny_conf.json',state_data)
# sm = pystan.StanModel(file="bayesSIRv1.1.stan")
# fit = sm.sampling(data=dat, iter=4000, chains=4)
# describe_fit(fit)

