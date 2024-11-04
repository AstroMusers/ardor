# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 15:43:08 2024

@author: Nate Whitsett
"""

import ardor.Flares.Flare as Flares
import numpy as np
import os
import pandas as pd
import warnings
import time
warnings.filterwarnings("ignore")
parameters = pd.read_csv('C:/Users/whitsett.n/Desktop/TOI_List.csv')
directory = 'C:/Users/whitsett.n/Desktop/TOI_Data/'
host_list = os.listdir(directory)
flare_dir = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TOI_Flare_csvs/'

for hosts in host_list:
    hosts = float(hosts)
    error_list = []
    const = 0
    try:
        input_dir = directory + str(hosts)
        output_dir = flare_dir
        try:
            os.mkdir(flare_dir + str(hosts), mode=0o777)
        except FileExistsError:
            print('Folder Exists! Skipping.')
            continue
        file_list = os.listdir(input_dir)
        if len(file_list) == 0:
            print('No Data! Skipping.')
        for index, files in enumerate(file_list):
            tier_0_timer = []
            tier_1_timer = []
            tier_2_timer = []
            print('LC #: ' + str(index) )
            Teff = float(parameters.loc[parameters['hostname'] == hosts]['st_teff'].item())
            period = float(parameters.loc[parameters['hostname'] == hosts]['pl_orbper'].item())
            radius = float(parameters.loc[parameters['hostname'] == hosts]['st_rad'].item())
            epoch = float(parameters.loc[parameters['hostname'] == hosts]['pl_orbtper'].item())
            if np.isnan(epoch) == True:
                epoch = float(parameters.loc[parameters['hostname'] == hosts]['pl_tranmid'].item())
                if np.isnan(epoch) == True:
                    epoch = np.nan
            tau1 = time.time()
            lc = Flares.tier0(input_dir + '/' + str(files))
            tau2 = time.time()
            tier_0_timer.append(tau2-tau1)                  
            flares = Flares.tier1(lc.detrended_flux, 3)
            tau3 = time.time()
            tier_1_timer.append(tau3 - tau2)
            ZZ, flare_num = Flares.tier2(lc.time, lc.detrended_flux, lc.error, flares.index, flares.length,
                          chi_square_cutoff=15, output_dir=output_dir, host_name=hosts,
                          T = Teff, host_radius=radius, csv=True, planet_period=period, 
                          planet_epoch=epoch, const = const, obs_time = lc.obs_time)
            const += flare_num
            tau4 = time.time()
            try:
                tier_2_timer.append((tau4 - tau3)/len(flares.index))
            except:
                tier_2_timer.append(np.nan)
            with open('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/TOI_T2_Flares.csv', "a") as f:
                np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                f.close()
            XX = np.column_stack((tier_0_timer, tier_1_timer, tier_2_timer))
            with open('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/ardor_time_stats.csv', "a") as f:
                np.savetxt(f, XX, delimiter=",", fmt='%s')
                f.close()
    except:
        error_list.append(hosts)
        print('Error in pipeline! Consider rerunning: ' + str(hosts))
        continue