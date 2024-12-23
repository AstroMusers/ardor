# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:42:31 2024

@author: whitsett.n
"""

import numpy as np
import pandas as pd
import os
from scipy.stats import uniform
from scipy.stats import ks_1samp
from skgof import ad_test
import csv
from itertools import zip_longest
from astropy.stats import kuiper
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

# sigma = 0.05
num = str(0.5)
# target_dir = 'M_Flaring_SPI_' + num + '.csv'


# for files in os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/G_Type_Inverse_Cubic'):
TOI_flares = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv')
KS_Sample_list = []
AD_Sample_list = []
N = 200

interest_hosts = [['Host_ID', 'p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N', 'dlog(Z)', 'Periastron?', 'Obs_Time', 'Period', 'Close_Approach']]
# interest_hosts = [['p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N']]
All_CDF = []
All_CDFx = []
All_Energy = []
set_hosts = set(TOI_flares['Host_ID'])

obs_time = []
host_list = []

df = pd.DataFrame()
for index, hosts in enumerate(set_hosts):
    name = str(hosts).replace(' ', '')
    planet_period = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Period'])[0]
    phases = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Transit_Phase'], dtype=float)
    phases1 = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Transit_Phase'], dtype=float)
    # obs_time = 700300
    peri_phases = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Phase'])
    KU_peri_phases = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Phase'])
    peri_phases_lower = np.abs(np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Lower'])[0])
    peri_phases_upper = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Upper'])[0]
    # observation_time = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Observation_Time'])[0]
    dlogz = np.median(np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'log(Z)']))
    energy = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Energy'])
    All_Energy.append(energy)
    if len(phases) < 5 or len(peri_phases) < 5:
        continue
    if np.isnan(np.mean(peri_phases)) == False:
        peri = True
    if np.isnan(np.mean(peri_phases)) == True:
        peri = False
    if np.isnan(np.mean(phases)) == True and peri == False:
        continue
    if np.isnan(peri_phases_lower) == True or np.isnan(peri_phases_upper) == True:
        error = False
    if np.isnan(peri_phases_lower) == False and np.isnan(peri_phases_upper) == False:
        error = True
    if peri == True:
        peri_phases = np.sort(peri_phases)
        D, p_KS = ks_1samp(peri_phases, uniform.cdf, args=(0, 1))
        A, p_AD = ad_test(peri_phases, uniform(0,1), assume_sorted=True)
        a, p_KU = kuiper(peri_phases, uniform.cdf, args=(0, 1))
        for phase_sample in range(N):
            if error == False:
                peri_phases = peri_phases + np.random.normal(0, scale=0.05)
            elif error == True:
                check = np.random.random()
                if check > 0.5:
                    peri_phases = peri_phases + np.abs(np.random.normal(0, scale=np.abs(peri_phases_upper)))
                elif check < 0.5:
                    peri_phases = peri_phases - np.abs(np.random.normal(0, scale=np.abs(peri_phases_lower)))
                for index in range(len(peri_phases)):
                    if peri_phases[index] > 1:
                        peri_phases[index] = peri_phases[index] - 1 
                    if peri_phases[index] < 0:
                        peri_phases[index] = peri_phases[index] + 1
                peri_phases = np.sort(peri_phases)
                D, p_ks_samp = ks_1samp(peri_phases, uniform.cdf, args=(0, 1))
                A, p_ad_samp = ad_test(peri_phases, uniform(0,1), assume_sorted=True)
                KS_Sample_list.append(p_ks_samp)
                AD_Sample_list.append(p_ad_samp)
    if peri == False:
        phases = np.sort(phases)
        D, p_KS = ks_1samp(phases, uniform.cdf, args=(0, 1))
        A, p_AD = ad_test(phases, uniform(0,1))
        a, p_KU = kuiper(phases, uniform.cdf, args=(0, 1))
        for phase_sample in range(N):
            phases += 1/N
            for index in range(len(phases)):
                if phases[index] > 1:
                    phases[index] = phases[index] - 1 
                if phases[index] < 0:
                    phases[index] = phases[index] + 1
            phases = np.sort(phases)
            D, p_ks_samp = ks_1samp(phases, uniform.cdf, args=(0, 1))
            A, p_ad_samp = ad_test(phases, uniform(0,1), assume_sorted=True)
            KS_Sample_list.append(p_ks_samp)
            AD_Sample_list.append(p_ad_samp)
    if p_KU < 1:
        interest_hosts.append([hosts, p_KU, p_KS, p_AD, np.median(KS_Sample_list), np.median(AD_Sample_list), len(peri_phases), np.mean(dlogz), peri])
        # interest_hosts.append([p_KU, p_KS, p_AD, np.median(KS_Sample_list), np.median(AD_Sample_list), len(peri_phases)])
        x = (np.sort(peri_phases)).tolist()
        y = (np.arange(len(x))/float(len(x))).tolist()
        x.insert(0, str(hosts) + '_x')
        y.insert(0, str(hosts) + '_y')
        x.append(1)
        y.append(1)
        All_CDF.append(x)
    KS_Sample_list = []
    AD_Sample_list = []

# with open("C:/Users/whitsett.n/Desktop/Host_CDFs.csv","w+", newline='') as f:
#     writer = csv.writer(f)
#     for values in zip_longest(*All_CDF):
#         writer.writerow(values)
with open("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/KS_KU_Tests/K_Tests/Exo/K_All_Exo.csv","w+", newline='') as f:
    writer = csv.writer(f)
    for values in interest_hosts:
        writer.writerow(values)
