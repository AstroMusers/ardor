# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 13:26:00 2024

@author: Nate Whitsett
"""

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from matplotlib import rcParams
def linear(x, alpha, beta):
    return alpha*x + beta
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 13,
        }
font_small ={'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
rcParams["mathtext.default"] = 'rm'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
data = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All TOIs/All_TOI_MCMC_Flares.csv')
catalog = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All TOIs/All_TOI_Reference.csv')
# alfven = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Alfven_Catalog_New.csv')
hosts = np.array(data['Host_ID'])
color = ['#d73027',
'#225ea8',]

color2 = ['#b2182b',
'#01665e']
color3 = ['#252525', '#bf812d']
custom_lines = [Line2D([], [], color='#d73027', lw=4)]
custom_lines2 = [Line2D([], [], color='#225ea8', lw=4)]
custom_lines3 =[Line2D([], [], color='#b2182b', lw=4),
                Line2D([], [], color = '#01665e', lw=4)]
st_class = ['M', 'KGF']
temp = [2000, 3850, 20000]
age = [0, 1, 20]
met = [-1, 0, 1]
alpha_list = []
beta_list = []
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,8), sharey = True, sharex = True)
fig.subplots_adjust(wspace=0, hspace = 0)
host_ = 0
for host in hosts:
    # periastron = np.array(pd.Series(data.loc[data['Host_ID'] == host, 'Periastron_Epoch']))[0]
    transit = np.array(pd.Series(data.loc[data['Host_ID'] == host, 'Transit_Epoch']))[0]
    if host != host_:
        energy_count_list = []
        energy_count_list2 = []
        energy_count_list3 = []
        err_list = []
        err_list2 = []
        err_list3 = []
        Teff = np.array(data.loc[data['Host_ID'] == host, 'Teff'])[0]
        obs_time = np.array(data.loc[data['Host_ID'] == host, 'Observation_Time'])[0]/(24*60)
        if obs_time == 0 or np.isnan(obs_time) == True:
            continue
        # energies = np.array(data.loc[(data['Host_Name'] == host) & (d), 'Flare_Energy'])
        energies = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'Energy']).dropna())
        energies2 = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'Energy']).dropna())
        dlogz = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'dlogZ']).dropna())
        metal = np.array(catalog.loc[(catalog['Host_Name'] == host), 'Metallicity'])[0]
        # try:
        #     ages = np.array(alfven.loc[(alfven['Host_ID'] == host), 'st_age'])[0]
        # except:
        #     continue
        energies[::-1].sort()
        energies = energies[~np.isnan(energies)]
        energies = np.log10(energies)
        energies2[::-1].sort()
        energies2 = energies2[~np.isnan(energies2)]
        energies2 = np.log10(energies2)
        for index, labels in enumerate(st_class):
            lower = temp[index]
            upper = temp[index+1]
            if len(energies) > 15 and Teff > lower and Teff < upper and np.median(dlogz) > 7:
                a = np.linspace(energies[0], energies[-1], num = 10)
                c = np.linspace(30, 37.5, num = 100)
                for bins in a:
                    b = sum(x >= bins for x in energies)
                    energy_count_list.append(np.log10(b/(obs_time)))
                    err_list.append(np.log10(np.sqrt(b/(obs_time))))
                coeff, pcov = curve_fit(linear, a, energy_count_list, sigma = err_list, absolute_sigma=True)
                print(coeff)
                if index == 0:
                    ax[0][0].scatter(a, energy_count_list, marker = '.', c = color[index], alpha = 0.75)
                    ax[0][0].plot(c, linear(c, coeff[0], coeff[1]), ':', color = color[index], alpha = 0.75)
                elif index == 1:
                    ax[0][1].scatter(a, energy_count_list, marker = '.', c = color[index], alpha = 0.75)
                    ax[0][1].plot(c, linear(c, coeff[0], coeff[1]), ':', color = color[index], alpha = 0.75)
                
        for index in range(len(met) - 1):
            lower = met[index]
            upper = met[index+1]
            if len(energies) > 5 and metal > lower and metal < upper and np.median(dlogz) > 7:
                a = np.linspace(energies[0], energies[-1], num = 5)
                c = np.linspace(30, 37.5, num = 100)
                for bins in a:
                    b = sum(x >= bins for x in energies)
                    energy_count_list2.append(np.log10(b/(obs_time)))
                    err_list2.append(np.log10(np.sqrt(b/(obs_time))))
                coeff, pcov = curve_fit(linear, a, energy_count_list2,sigma = err_list2, absolute_sigma=True)
                print(coeff)
                if index == 0:
                    ax[1][0].scatter(a, energy_count_list2, marker = '.', c = color2[index], alpha = 0.75)
                    ax[1][0].plot(c, linear(c, coeff[0], coeff[1]), ':', color = color2[index], alpha = 0.75)
                elif index == 1:
                    ax[1][1].scatter(a, energy_count_list2, marker = '.', c = color2[index], alpha = 0.75)
                    ax[1][1].plot(c, linear(c, coeff[0], coeff[1]), ':', color = color2[index], alpha = 0.75)
        # for index in range(len(age) - 1):
        #     lower = age[index]
        #     upper = age[index+1]
        #     if len(energies) > 15 and ages > lower and ages < upper and np.median(dlogz) > 7:
        #         a = np.linspace(energies[0], energies[-1], num = 10)
        #         c = np.linspace(30, 37.5, num = 100)
        #         for bins in a:
        #             b = sum(x >= bins for x in energies)
        #             energy_count_list3.append(np.log10(b/(obs_time)))
        #             err_list3.append(np.log10(np.sqrt(b/(obs_time))))
        #         coeff, pcov = curve_fit(linear, a, energy_count_list3,sigma = err_list3, absolute_sigma=True)
        #         print(coeff)
        #         if index == 0:
        #             ax[2][0].scatter(a, energy_count_list3, marker = '.', c = color3[index], alpha = 0.75)
        #             ax[2][0].plot(c, linear(c, coeff[0], coeff[1]), ':', color = color3[index], alpha = 0.75)
        #         elif index == 1:
        #             ax[2][1].scatter(a, energy_count_list3, marker = '.', c = color3[index], alpha = 0.75)
        #             ax[2][1].plot(c, linear(c, coeff[0], coeff[1]), ':', color = color3[index], alpha = 0.75)
            
    host_ = host

ax[0][0].set_ylabel('log Flare Rate $\geq E_{\mathrm{bol}} (\mathrm{d}^{-1}) $', fontdict = font)
ax[1][0].set_ylabel('log Flare Rate $\geq E_{\mathrm{bol}} (\mathrm{d}^{-1}) $', fontdict = font)
# ax[2][0].set_ylabel('log Flare Rate $\geq E_{\mathrm{bol}} (\mathrm{d}^{-1}) $', fontdict = font)
ax[0][0].xaxis.set_label_position('top')
ax[0][0].xaxis.tick_top()
ax[0][1].xaxis.set_label_position('top')
ax[0][0].tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)
ax[0][1].tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)
ax[1][1].tick_params(top=True, labeltop=False, bottom=True, labelbottom=True)
ax[1][0].tick_params(top=True, labeltop=False, bottom=True, labelbottom=True)
ax[0][1].xaxis.tick_top()
# ax[1][0].xaxis.tick_top()
# ax[1][1].xaxis.tick_top()

ax[0][1].set_xlabel('log $(E_{\mathrm{bol}})$', fontdict = font)
ax[0][0].set_xlabel('log $(E_{\mathrm{bol}})$', fontdict = font)
ax[1][0].set_xlabel('log $(E_{\mathrm{bol}})$', fontdict = font)
# ax[2][1].set_xlabel('log $(E_{\mathrm{bol}})$', fontdict = font)
# ax[2][0].set_xlabel('log $(E_{\mathrm{bol}})$', fontdict = font)
ax[0][0].set_ylim(-3.9, 1)
ax[0][1].set_ylim(-3.9, 1)
ax[0][0].set_xlim(30.75, 37.5)
ax[0][1].set_xlim(30.75, 37.5)
ax[1][0].set_xlim(30.75, 37.5)
ax[1][1].set_xlim(30.75, 37.5)
# ax[2][0].set_xlim(30.75, 37.5)
# ax[2][1].set_xlim(30.75, 37.5)
# ax[0][0].legend(custom_lines, ['M8-M5', 'M4-M0'], prop={"family":"serif"})
ax[1][1].set_xlabel('log $(E_{\mathrm{bol}})$', fontdict = font)
# ax[0][1].legend(custom_lines2, ['FGK'], prop={"family":"serif"})
ax[0][0].text(34.55, .4, '$T_{eff} \leq 3850 \; \mathrm{K}$', size = 14)
ax[0][1].text(34.5, .4, '$T_{eff} > 3850 \; \mathrm{K}$', size = 14)
ax[1][0].text(35.5, .4, '$\mathrm{\\frac{[Fe]}{H}} \leq 0$',size = 14)
ax[1][1].text(35.5, .4, '$\\frac{[Fe]}{H} \geq 0$',size = 14)
# ax[2][0].text(34.75, 1.4, '$\mathrm{t_{Age}} \leq 1 \; \mathrm{Gyr}$', size = 14)
# ax[2][1].text(34.75, 1.4, '$\mathrm{t_{Age}} \geq 1 \; \mathrm{Gyr}$', size = 14)
ax[0][0].text(31, -3.5, '(a)', fontdict= font_small)
ax[0][1].text(31, -3.5, '(b)', fontdict= font_small)
ax[1][0].text(31, -3.5, '(c)',fontdict= font_small)
ax[1][1].text(31, -3.5, '(d)',fontdict= font_small)
# ax[2][0].text(31, -2.5, '(e)',fontdict= font_small)
# ax[2][1].text(31, -2.5, '(f)',fontdict= font_small)
plt.show()
fig.savefig('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Publication Documents/Publication_2024/FFD_TOI.png', dpi=400, bbox_inches = 'tight')