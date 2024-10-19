# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 20:41:03 2023

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import uniform
from scipy.stats import ks_1samp
TOI_flares = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_TOI_Flares.csv")
Exoplanet_flares = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_Exoplanet_MCMC_Flares.csv")
count = 0
p_values = []
center_statistic = []
test = []

fig, ax = plt.subplots(3,3, figsize = (12,12))
for hosts in np.array(TOI_flares['Host_Name']):
    name = hosts
    planet_period = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Period'])[0]
    phases = np.array(TOI_flares.loc[TOI_flares['Host_Name'] == hosts, 'Phase'])/planet_period
    count += 1
    r = uniform.rvs(size=50)
    if count == len(phases):
        a = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
        p_values.append(a.pvalue)
        center_statistic.append(a.statistic_location)
        if count > 10 and a.pvalue < 0.05 and planet_period < 100:
            x = np.sort(phases)
            y = np.arange(len(phases))/float(len(phases))
            ax[2, 0].plot(x,y,label=hosts)
            
        count = 0
ax[2, 1].hist(center_statistic, bins=np.linspace(0,1, num=10), color='r', edgecolor='black', linewidth=1.2)
ax[2, 2].hist(p_values, bins=np.linspace(0,1, num=10), color='r',edgecolor='black', linewidth=1.2)

count = 0
p_values = []
center_statistic = []
test = []
for hosts in np.array(Exoplanet_flares['Host_Name']):
    name = hosts
    planet_period = np.array(Exoplanet_flares.loc[Exoplanet_flares['Host_Name'] == hosts, 'Period'])[0]
    phases = np.array(Exoplanet_flares.loc[Exoplanet_flares['Host_Name'] == hosts, 'Phase'])/planet_period
    count += 1
    r = uniform.rvs(size=50)
    if count == len(phases):
        a = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
        p_values.append(a.pvalue)
        center_statistic.append(a.statistic_location)
        if count > 10 and a.pvalue < 0.05 and planet_period < 100:
            x = np.sort(phases)
            y = np.arange(len(phases))/float(len(phases))
            ax[0, 0].plot(x,y,label=hosts)
            
        count = 0
        
        
ax[0,0].axvline(x=0.5, color = 'c', alpha = 0.5, linestyle = ':', linewidth = 1.75)
ax[0, 1].hist(center_statistic, bins=np.linspace(0,1, num=10), color='g',edgecolor='black', linewidth=1.2)
ax[0, 2].hist(p_values, bins=np.linspace(0,1, num=10), color='g',edgecolor='black', linewidth=1.2)

count = 0
p_values = []
center_statistic = []
test = []
for hosts in np.array(Exoplanet_flares['Host_Name']):
    name = hosts
    planet_period = np.array(Exoplanet_flares.loc[Exoplanet_flares['Host_Name'] == hosts, 'Period'])[0]
    phases = np.array(Exoplanet_flares.loc[Exoplanet_flares['Host_Name'] == hosts, 'Period_Phase'])/planet_period
    count += 1
    r = uniform.rvs(size=50)
    if count == len(phases):
        a = ks_1samp(phases, uniform.cdf, args=(0, 1), alternative='two-sided')
        p_values.append(a.pvalue)
        center_statistic.append(a.statistic_location)
        if count > 10 and a.pvalue < 0.05 and planet_period < 100:
            x = np.sort(phases)
            y = np.arange(len(phases))/float(len(phases))
            ax[1, 0].plot(x,y,label=hosts)
        count = 0
ax[1, 1].hist(center_statistic, bins=np.linspace(0,1, num=10), color='c',edgecolor='black', linewidth=1.2)
ax[1, 2].hist(p_values, bins=np.linspace(0,1, num=10), color='c',edgecolor='black', linewidth=1.2)


r = uniform.cdf(np.linspace(0, 1))
x = np.sort(r)
y = np.arange(len(x))/float(len(x))
ax[0, 0].plot(x, y, label='Uniform Distribution', linestyle = '--', linewidth=1.5, color = 'black')
ax[1, 0].plot(x, y, label='Uniform Distribution', linestyle = '--', linewidth = 1.5, color='black')
ax[2, 0].plot(x, y, label='Uniform Distribution', linestyle = '--', linewidth = 1.5, color='black')

ax[0,0].set_title('Flare Phase CDF ($p \leq 0.05$)', fontsize = 'x-large')
ax[0,1].set_title('Statistic Location Distribution', fontsize = 'x-large')
ax[0,2].set_title('p-value Distribution', fontsize = 'x-large')


ax[0,0].text(0.25, 0.75, 'Periastron Phase', horizontalalignment='center', transform=ax[0,0].transAxes, fontsize = 'medium', color = 'c')


fig.text(0.08, 0.68, 'Constrained Periastron', rotation='vertical', fontsize = 'x-large')
fig.text(0.08, 0.41, 'Other Exoplanet Hosts', rotation='vertical', fontsize = 'x-large')
fig.text(0.08, 0.2, 'All TOI Hosts', rotation='vertical', fontsize = 'x-large')

# ax[0,0].set_xlabel('Planet Phase')
# ax[0,0].set_ylabel('Cumulative Probability', fontsize = 'x-large')
ax[2,0].set_xlabel('Planet Phase', fontsize = 'x-large')
# ax[1,0].set_ylabel('Cumulative Probability', fontsize = 'large')

# ax[0,1].set_xlabel('Statistic Location')
# ax[0,1].set_ylabel('Count')
ax[2,1].set_xlabel('Statistic Location', fontsize = 'x-large')
# ax[1,1].set_ylabel('Count')

# ax[0,2].set_xlabel('p-value')
# ax[0,2].set_ylabel('Count')
ax[2,2].set_xlabel('p-value', fontsize = 'x-large')
# ax[1,2].set_ylabel('Count')

plt.subplots_adjust(wspace = 0.15, hspace = 0.1)

fig.savefig('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Publication Documents/CDF.png', dpi=fig.dpi)