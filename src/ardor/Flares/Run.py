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

### HOSTS TIER 2 RUN ###
# parameters = pd.read_csv('/ugrad/whitsett.n/Flare_Data/Reference_Files/Alfven_Catalog.csv')
# directory = '/ugrad/whitsett.n/data/Hosts/'
# host_list = os.listdir(directory)
# flare_dir = '/ugrad/whitsett.n/Flare_Data/Tier_2/TOI/'

# lists = os.listdir(flare_dir)
# for hosts in lists:
#     if len(os.listdir(flare_dir + hosts)) == 0:
#         os.rmdir(flare_dir + hosts)
# for hosts in host_list:
#     error_list = []
#     const = 0
#     try:
#         input_dir = directory + str(hosts)
#         output_dir = flare_dir
#         try:
#             os.mkdir(flare_dir + str(hosts), mode=0o777)
#         except FileExistsError:
#             print('Folder Exists! Skipping.')
#             continue
#         file_list = os.listdir(input_dir)
#         if len(file_list) == 0:
#             print('No Data! Skipping.')
#         for index, files in enumerate(file_list):
#             tier_0_timer = []
#             tier_1_timer = []
#             tier_2_timer = []
#             print('LC #: ' + str(index) )
#             try:
#                 Teff = float(parameters.loc[parameters['Host_ID'] == hosts]['st_teff'].item())
#             except:
#                 Teff = np.nan
#             try:
#                 period = float(parameters.loc[parameters['Host_ID'] == hosts]['pl_orbper'].item())
#             except:
#                 period = np.nan
#             try:
#                 radius = float(parameters.loc[parameters['Host_ID'] == hosts]['st_rad'].item())
#             except:
#                 radius = np.nan
#             try:
#                 peri_epoch = float(parameters.loc[parameters['Host_ID'] == hosts]['pl_orbtper'].item())
#             except:
#                 peri_epoch = np.nan
#             try:
#                 trans_epoch = float(parameters.loc[parameters['Host_ID'] == hosts]['pl_tranmid'].item())
#             except:
#                 trans_epoch = np.nan
#             tau1 = time.time()
#             lc = Flares.tier0(input_dir + '/' + str(files))
#             tau2 = time.time()
#             tier_0_timer.append(tau2-tau1)                  
#             flares = Flares.tier1(lc.detrended_flux, 3)
#             tau3 = time.time()
#             tier_1_timer.append(tau3 - tau2)
#             ZZ, flare_num = Flares.tier2(lc.time, lc.detrended_flux, lc.error, flares.index, flares.length,
#                             chi_square_cutoff=20, output_dir=output_dir, host_name=hosts,
#                             T = Teff, host_radius=radius, csv=True, planet_period=period, 
#                             transit_epoch=trans_epoch, peri_epoch = peri_epoch, const = const)
#             const += flare_num
#             tau4 = time.time()
#             try:
#                 tier_2_timer.append((tau4 - tau3)/len(flares.index))
#             except:
#                 tier_2_timer.append(np.nan)
#             with open('/ugrad/whitsett.n/Flare_Data/Flare_Lists/Hosts_T2_Flares.csv', "a") as f:
#                 np.savetxt(f, ZZ, delimiter=",", fmt='%s')
#                 f.close()
#             XX = np.column_stack((tier_0_timer, tier_1_timer, tier_2_timer))
#             with open('/ugrad/whitsett.n/Flare_Data/Flare_Lists/ardor_time_stats.csv', "a") as f:
#                 np.savetxt(f, XX, delimiter=",", fmt='%s')
#                 f.close()
#     except:
#         error_list.append(hosts)
#         print('Error in pipeline! Consider rerunning: ' + str(hosts))
#         continue


parameters = pd.read_csv('/ugrad/whitsett.n/Flare_Data/Reference_Files/TOI_List.csv')
directory = '/ugrad/whitsett.n/data/TOIs/'
host_list = os.listdir(directory)
flare_dir = '/ugrad/whitsett.n/Flare_Data/Tier_2/TOI/'

lists = os.listdir(flare_dir)
for hosts in lists:
    if len(os.listdir(flare_dir + hosts)) == 0:
        os.rmdir(flare_dir + hosts)
for hosts in host_list:
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
            try:
                Teff = float(parameters.loc[parameters['Host_ID'] == hosts]['st_teff'].item())
            except:
                Teff = np.nan
            try:
                period = float(parameters.loc[parameters['Host_ID'] == hosts]['pl_orbper'].item())
            except:
                period = np.nan
            try:
                radius = float(parameters.loc[parameters['Host_ID'] == hosts]['st_rad'].item())
            except:
                radius = np.nan
            try:
                trans_epoch = float(parameters.loc[parameters['Host_ID'] == hosts]['pl_tranmid'].item())
            except:
                trans_epoch = np.nan
            tau1 = time.time()
            lc = Flares.tier0(input_dir + '/' + str(files))
            tau2 = time.time()
            tier_0_timer.append(tau2-tau1)                  
            flares = Flares.tier1(lc.detrended_flux, 3)
            tau3 = time.time()
            tier_1_timer.append(tau3 - tau2)
            ZZ, flare_num = Flares.tier2(lc.time, lc.detrended_flux, lc.error, flares.index, flares.length,
                            chi_square_cutoff=20, output_dir=output_dir, host_name=hosts,
                            T = Teff, host_radius=radius, csv=True, planet_period=period, 
                            transit_epoch=trans_epoch, const = const)
            const += flare_num
            tau4 = time.time()
            try:
                tier_2_timer.append((tau4 - tau3)/len(flares.index))
            except:
                tier_2_timer.append(np.nan)
            with open('/ugrad/whitsett.n/Flare_Data/Flare_Lists/TOI_T2_Flares.csv', "a") as f:
                np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                f.close()
            XX = np.column_stack((tier_0_timer, tier_1_timer, tier_2_timer))
            with open('/ugrad/whitsett.n/Flare_Data/Flare_Lists/ardor_time_stats.csv', "a") as f:
                np.savetxt(f, XX, delimiter=",", fmt='%s')
                f.close()
    except:
        error_list.append(hosts)
        print('Error in pipeline! Consider rerunning: ' + str(hosts))
        continue