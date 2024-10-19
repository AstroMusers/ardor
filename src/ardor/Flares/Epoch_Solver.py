# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 22:44:59 2024

@author: Nate Whitsett
"""

import pandas as pd
import numpy as np

T3_data = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Grand_Lists/All_Exoplanet_Flares_T3.csv')
MCMC_data = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Grand_Lists/All_Exoplanets_MCMC_Flares.csv')

new_epochs = []
MCMC_Hosts = MCMC_data['Host_Name']
MCMC_Flares = MCMC_data['Flare_#']
MCMC_Epochs = MCMC_data['Flare_Epoch']
T3_Hosts = T3_data['HostID']
T3_Flares = T3_data['#']
T3_Epochs = T3_data['Flare_Epoch']
for indexs in range(len(MCMC_data)):
    MCMC_Host = MCMC_Hosts[indexs]
    MCMC_Flare = MCMC_Flares[indexs]
    for indexes in range(len(T3_Hosts)):
        T3_Host = T3_Hosts[indexes]
        T3_Flare = T3_Flares[indexes]
        T3_Epoch = T3_Epochs[indexes]
        
        if T3_Host == MCMC_Host and T3_Flare == MCMC_Flare:
            T3_Hosts.drop(indexes, inplace = True)
            T3_Hosts.reset_index(drop = True,inplace = True)
            T3_Flares.drop(indexes, inplace = True)     
            T3_Flares.reset_index(drop = True,inplace = True)
            T3_Epochs.drop(indexes, inplace = True)
            T3_Epochs.reset_index(drop = True,inplace = True)
            MCMC_Epochs[indexs] = T3_Epoch
            break
    print(len(T3_Flares))
    print(indexs/len(MCMC_data)*100)
MCMC_data.to_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Grand_Lists/All_Exoplanets_MCMC_Flares_Epochs.csv', index=False)