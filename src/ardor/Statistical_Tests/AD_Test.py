# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:42:31 2024

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
import os
from skgof import ad_test
from matplotlib import pyplot as plt
from scipy.stats import uniform
from scipy.stats import ks_1samp
import csv
from itertools import zip_longest
import ardor.Flares.Flare as Flare
from astropy.stats import kuiper
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

# sigma = 0.05
num = str(0.5)
# target_dir = 'M_Flaring_SPI_' + num + '.csv'
TOI_flares = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_New.csv')
count = 0
kuiper_list = []
KS_list = []
KS_Sample_list = []
test = []
count1 = 0
N = 1000

interest_hosts = []
host_num = 0
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
    if np.isnan(np.mean(peri_phases)) == True:
        peri = False
    if np.isnan(np.mean(phases)) == True and peri == False:
        continue
    if np.isnan(peri_phases_lower) == True or np.isnan(peri_phases_upper) == True:
        error = False
    if np.isnan(peri_phases_lower) == False and np.isnan(peri_phases_upper) == False:
        error = True
    if peri == True:
        for phase_sample in range(N):
            if peri == True and error == False:
                peri_phases = peri_phases + np.random.normal(0, scale=0.03)
            elif peri == True and error == True:
                check = np.random.random()
                if check > 0.5:
                    peri_phases = peri_phases + np.abs(np.random.normal(0, scale=np.abs(peri_phases_upper)))
                elif check < 0.5:
                    peri_phases = peri_phases - np.abs(np.random.normal(0, scale=np.abs(peri_phases_lower)))
                D, p_ks_samp = ks_1samp(peri_phases, uniform.cdf, args=(0, 1))
                
    # if peri == False:
    #     # for index in range(len(phases)):
    #     #     if phases[index] > 1:
    #     #         phases[index] = phases[index] - 1 
    #     #     if phases[index] < 0:
    #     #         phases[index] = phases[index] + 1
    #     phases = np.sort(phases)
    #     a, p = kuiper(phases, uniform.cdf, args=(0, 1))
    if peri == True:
        # for index in range(len(phases)):
        #     if peri_phases[index] > 1:
        #         peri_phases[index] = peri_phases[index] - 1 
        #     if peri_phases[index] < 0:
        #         peri_phases[index] = peri_phases[index] + 1
        peri_phases = np.sort(peri_phases)
        a, p_ku = kuiper(peri_phases, uniform.cdf, args=(0, 1))
        a, p_ks = ks_1samp(peri_phases, uniform.cdf, args=(0, 1))
    if p_ku < 1:
        x = (np.sort(phases)).tolist()
        y = (np.arange(len(x))/float(len(x))).tolist()
        x.insert(0, str(hosts) + '_x')
        y.insert(0, str(hosts) + '_y')
        x.append(1)
        y.append(1)
        All_CDF.append(x)
        All_CDF.append(y)
        test = []
        host_num += 1
        interest_hosts.append([hosts, p_ku, p_ks, len(x), np.mean(dlogz)])
    kuiper_list = []
    KS_list = []

# with open("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Host_CDFs.csv","w+", newline='') as f:
#     writer = csv.writer(f)
#     for values in zip_longest(*All_CDF):
#         writer.writerow(values)
with open("C:/Users/Nate Whitsett/Desktop/KS_AD_Hosts_Peri.csv","w+", newline='') as f:
    writer = csv.writer(f)
    for values in interest_hosts:
        writer.writerow(values)

        




























# M005 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_0.05.csv', header=None)
# M01 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_0.1.csv', header=None)
# M05 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_0.5.csv', header=None)
# M1 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_1.csv', header=None)
# M2 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_2.csv', header=None)
# M5 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_5.csv', header=None)
# M10 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/M_type_SPI_Sim/M_type_SPI_Sim_CDF_10.csv', header=None)
# MCDFS = [M005, M01, M05, M5, M10]
# label = [0.05, 0.1, 0.5, 5, 10]
# index = 0
# G005 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_0.05.csv', header=None)
# G01 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_0.1.csv', header=None)
# G05 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_0.5.csv', header=None)
# G1 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_1.csv', header=None)
# G2 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_2.csv', header=None)
# G5 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_5.csv', header=None)
# G10 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/G_type_SPI_Sim/G_type_SPI_Sim_CDF_10.csv', header=None)
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
# plt.savefig('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/eCDF_Sim.png', dpi=400, bbox_inches='tight')
# plt.show()