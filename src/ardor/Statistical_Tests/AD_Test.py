# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:42:31 2024

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
import os
from skgof import ks_test, ad_test
from matplotlib import pyplot as plt
from scipy.stats import uniform
from scipy.stats import ks_1samp
from scipy.stats import ksone
import csv
from itertools import zip_longest
import ardor.Flares.Flare as Flare

def ks_critical_value(n_trials, alpha):
    return ksone.ppf(1-alpha/2, n_trials)
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

# sigma = 0.05
num = str(0.5)
# target_dir = 'M_Flaring_SPI_' + num + '.csv'
TOI_flares = pd.read_csv('C:/Users/Nate Whitsett/Desktop/All_Exoplanet_MCMC_Flares.csv')
# alfven = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Alfven_Catalog_New.csv")
# Exoplanet_flares = pd.read_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_Exoplanet_MCMC_Flares.csv")
# lower_phases = alfven["Sub_Alfv_lphase"]
# hosts_list = alfven["hostname"].to_list()
const_list = []
# for index, phases in enumerate(lower_phases):
#     if np.isnan(phases) == False:
#         const_list.append(hosts_list[index])
# targets = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Final_Targets.csv")
count = 0
KS_p_values = []
AD_p_values = []
test = []
count1 = 0
# fig, ax = plt.subplots(2,2)
N = 50

interest_hosts = []
host_num = 0
All_CDF = []
All_CDFx = []
All_Energy = []
set_hosts = set(TOI_flares['Host_ID'])

obs_time = []
host_list = []
# directory = os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/All_Exoplanet_Hosts')


# hosts2 = targets["Host_ID"].tolist() 
for hosts in set_hosts:
    # if str(hosts) in hosts2:
        print(hosts)
        name = str(hosts).replace(' ', '')
        # planet_period = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Orb_Per'])[0]
        phases = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Transit_Phase'], dtype=float)
        phases1 = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Transit_Phase'], dtype=float)
        # obs_time = 700300
        peri_phases = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Phase'])
        peri_phases_lower = np.abs(np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Lower'])[0])
        peri_phases_upper = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Periastron_Upper'])[0]
        # observation_time = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Observation_Time'])[0]
        dlogz = np.median(np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'log(Z)']))
        energy = np.array(TOI_flares.loc[TOI_flares['Host_ID'] == hosts, 'Energy'])
        All_Energy.append(energy)
        if len(phases) < 3 or len(peri_phases) < 3:
            continue
        if np.isnan(np.mean(peri_phases)) == False:
            peri = True
            phases = peri_phases    
        if np.isnan(np.mean(peri_phases)) == True:
            peri = False
        if np.isnan(peri_phases_lower) == True or np.isnan(peri_phases_upper) == True:
            error = False
        if np.isnan(peri_phases_lower) == False and np.isnan(peri_phases_upper) == False:
            error = True
        if np.isnan(np.mean(phases)) == True and peri == False:
            continue
        for phase_sample in range(N):
            if peri == False:
                phases = phases + (1/N)
            elif peri == True and error == False:
                peri_phases = peri_phases + np.random.normal(0, scale=0.05)
            elif peri == True and error == True:
                check = np.random.random()
                if check > 0.5:
                    peri_phases = peri_phases + np.abs(np.random.normal(0, scale=peri_phases_upper))
                elif check < 0.5:
                    peri_phases = peri_phases - np.abs(np.random.normal(0, scale=peri_phases_lower))
            for index in range(len(phases)):
                if phases[index] > 1:
                    phases[index] = phases[index] - 1 
                if phases[index] < 0:
                    phases[index] = phases[index] + 1
            phases = np.sort(phases)
            a = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
            b = ad_test(phases, uniform(0,1), assume_sorted=True)
            KS_p_values.append(a.pvalue)
            AD_p_values.append(b.pvalue)
        x = np.sort(phases)
        y = np.arange(len(x))/float(len(x))
        All_CDF.append(x)
        All_CDF.append(y)
        test = []
        host_num += 1
        phases = phases1
        interest_hosts.append([hosts,np.median(KS_p_values), np.median(AD_p_values)])
        KS_p_values = []
        AD_p_values = []
# elif count == len(peri_phases) and hosts != 'AUMic b':
#     count = 0
# with open("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_type_SPI_Sim/eCDFS/M_Type_SPI_Sim_CDF_" + num + ".csv","w+", newline='') as f:
#     writer = csv.writer(f)
#     for values in zip_longest(*All_CDF):
#         writer.writerow(values)
with open("C:/Users/Nate Whitsett/Desktop/KS_AD.csv","w+", newline='') as f:
    writer = csv.writer(f)
    # writer.writerow(['KS p Values', 'KS Test Statistic', 'AD p Values', 'AD Test Statistic'])
    for values in interest_hosts:
        writer.writerow(values)

        




























# M005 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_0.05.csv', header=None)
# M01 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_0.1.csv', header=None)
# M05 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_0.5.csv', header=None)
# M1 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_1.csv', header=None)
# M2 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_2.csv', header=None)
# M5 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_5.csv', header=None)
# M10 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_10.csv', header=None)
# MCDFS = [M005, M01, M05, M5, M10]
# label = [0.05, 0.1, 0.5, 5, 10]
# index = 0
# G005 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_0.05.csv', header=None)
# G01 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_0.1.csv', header=None)
# G05 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_0.5.csv', header=None)
# G1 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_1.csv', header=None)
# G2 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_2.csv', header=None)
# G5 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_5.csv', header=None)
# G10 = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_10.csv', header=None)
# GCDFS = [G005, G01, G05, G5, G10]

# fig = plt.figure(figsize=(7,7))
# gs = fig.add_gridspec(2, 2)
# ax1 = plt.subplot(gs[0, :])
# ax1.plot(np.linspace(0,1), np.linspace(0,1), linestyle='--', c='red', label='Uniform')
# ax2 = plt.subplot(gs[1, 0], sharex=ax1)
# ax2.plot(np.linspace(0,1), np.linspace(0,1), linestyle='--', c='red', label='Uniform')
# ax3 = plt.subplot(gs[1, 1])
# ax3.plot(np.linspace(0,1), np.linspace(0,1), linestyle='--', c='red', label='Uniform')
# ax1.set_xlabel('Orbital Phase',fontdict=font)
# ax2.set_xlabel('Orbital Phase',fontdict=font)
# ax3.set_xlabel('Orbital Phase', fontdict=font)
# ax1.set_ylabel('Cumulative Frequency',fontdict=font)
# ax2.set_ylabel('Cumulative Frequency',fontdict=font)
# plt.subplots_adjust(wspace=.0)
# plt.setp(ax3.get_yticklabels(), visible=False)
# ax1.plot(M005[6], M005[7], label = 'eCDF')
# ax1.axvspan(0.45, 0.55,alpha = 0.5, color = 'green')
# ax1.text(0.2, 0.55, 'Sub-Alfvenic', c = 'green', fontdict=font)
# ax1.text(0.58, 0.67, 'D', c = 'black', fontdict = font, size = 16)
# ax1.text(-.03,0.93, '(a)', fontdict=font, size=12)
# ax1.axvline(0.555, 0.555, 0.735, c = 'black', linewidth=2)
# ax1.scatter(0.555, 0.555, s=40, zorder=10, c='black')
# ax1.scatter(0.555, 0.765, s=40, zorder=10, c='black')
# ax2.text(-.02,0.93, '(b)', fontdict=font, size=12)
# ax3.text(-.02,0.93, '(c)', fontdict=font, size=12)
# index = 0
# for CDFS in MCDFS:
#     ax2.plot(CDFS[0], CDFS[1], label = '$\sigma = $' + str(label[index]))
#     index +=1 
# index = 0
# for CDFS in GCDFS:
#     ax3.plot(CDFS[4], CDFS[5],label = '$\sigma = $' + str(label[index]))
#     index +=1 
# L = ax1.legend(ncol=1, borderpad=0.4, labelspacing = 0.2, loc='lower right')
# plt.setp(L.texts, family='Serif', size = 14)
# L1 = ax2.legend(ncol=1, borderpad=0.2, labelspacing = 0.1, loc='lower right')
# plt.setp(L1.texts, family='Serif', size = 10)
# plt.savefig('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/eCDF_Sim.png', dpi=400, bbox_inches='tight')
# plt.show()