# -*- coding: utf-8 -*-
"""
Created on Sun May 26 00:37:47 2024

@author: Nate Whitsett
"""

import pandas as pd
import numpy as np

periastron = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_Parameter_Reference.csv")
data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_MCMC_Flares.csv")


hosts = periastron['Host_Name_Time'].tolist()
data_hosts = data['Host_ID'].tolist()
# # time_hosts = data['Host_Name_Time']
# obs_time = []
# alfven_list = []
# alfven_list_l = []
# alfven_list_u = []
# lower_phase_list = []
# upper_phase_list = []
observation_time_list = []
# # count = 0
# found = False
for index, data_stars in enumerate(data_hosts):
    # try:
    data_stars = data_stars.replace(' ', '')
    observation_time = np.array(periastron.loc[(periastron['Host_Name_Time'] == str(data_stars)), 'Obs_Time (min)'])[0]
    observation_time_list.append(observation_time)
        # period = np.array(periastron.loc[(data['Host_ID'] == data_stars), 'pl_orbper'])[0]
        # err1 = np.abs(np.array(periastron.loc[(data['Host_ID'] == data_stars) , 'pl_orbtpererr1'])[0]/period)
        # err2 =  np.abs(np.array(periastron.loc[(data['Host_ID'] == data_stars) , 'pl_orbtpererr2'])[0]/period)
        # lower_err.append(err2)
        # upper_err.append(err1)
        # rad = np.array(periastron.loc[periastron['hostname'] == data_stars, 'Alfven_Rad'])[0]
        # rad_l = np.array(periastron.loc[periastron['hostname'] == data_stars, 'Alfven_Rad_Lower'])[0]
        # rad_u = np.array(periastron.loc[periastron['hostname'] == data_stars, 'Alfven_Rad_Upper'])[0]
        # phase_l = np.array(periastron.loc[periastron['hostname'] == data_stars, 'Sub_Alfv_lphase'])[0]
        # phase_u = np.array(periastron.loc[periastron['hostname'] == data_stars, 'Sub_Alfv_uphase'])[0]
        # alfven_list.append(rad)
        # alfven_list_l.append(rad_l)
        # alfven_list_u.append(rad_u)
        # lower_phase_list.append(phase_l)
        # upper_phase_list.append(phase_u)
    # except:
    #     lower_err.append(np.nan)
    #     upper_err.append(np.nan)
    # else:
    #     alfven_list.append(np.nan)
    #     alfven_list_l.append(np.nan)
    #     alfven_list_u.append(np.nan)
    #     lower_phase_list.append(np.nan)
    #     upper_phase_list.append(np.nan)
# data['Alfven_Rad'] = alfven_list
# data['Alfven_Rad_Lower'] = alfven_list_l
# data['Alfven_Rad_Upper'] = alfven_list_u
# data['Lower_Phase'] = lower_phase_list
# periastron['Flare_Epoch_BJD'] = epochs
data['Obs_Time'] = observation_time_list
data.to_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_MCMC_Flares2.csv", index = False)