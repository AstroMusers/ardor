# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:08:22 2024

@author: natha
"""

import Flare
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import gca
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 15,
        }
font_label = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 18,
        }
font_small ={'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

targets = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Final_Targets.csv")
flares = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_MCMC_Flares.csv")
flares_TOI = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All TOIs/All_TOI_MCMC_Flares.csv")
params = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/Final_Targets.csv")


hosts = targets["Host_ID"].tolist()
directory = os.listdir("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/Final_Data")
print(directory)
directory = ['2222','6708', 'TOI-833', 'HAT-P-2', 'HD 215497', 'HD 76700']
shift = 0
shift2 = 0
color = ('#d7191c',
'#fdae61',
'#0571b0',
'#7b3294',
'#1a9641','#d01c8b', '#018571', '#a6611a', '#a50026', '#006837')
fig, ax = plt.subplots()
fig.set_size_inches(10, 6)
ax.set_xlabel('BJD - 2450000', fontdict = font_label)

ax.get_yaxis().set_ticks([])
ax.set_ylabel('Rel. Flux', fontdict=font_label)
min_time = 2000
max_time = 0
c_index = 0
counter2 = 3
counter3 = 3
num = 6
count = 0
for stars in directory:
    star = stars.replace(' ', '')
    counter = 1
    if star in hosts:
        print(1)
        if count < 7:
            
            lower_phase = np.array(targets.loc[targets['Host_ID'] == star, 'Sub_Alfv_lphase'])[0]
            period = np.array(targets.loc[targets['Host_ID'] == star, 'Period'])[0]
            periastron_epoch = np.array(targets.loc[targets['Host_ID'] == star, 'Periastron_Epochs'])[0]
            flare_epochs = np.array(flares.loc[flares['Host_ID'] == star, 'Flare_Epoch_TESS'])
            
            if np.isnan(periastron_epoch) == True:
                peri = False
                try:
                    flare_epochs2 = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == int(star), 'Flare_Epoch'])
                except:
                    flare_epochs2 = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == star, 'Flare_Epoch'])
            elif np.isnan(periastron_epoch) == False:
                peri = True
                flare_epochs = np.array(flares.loc[flares['Host_ID'] == star, 'Flare_Epoch_TESS'])
            upper_phase = np.array(params.loc[params['Host_ID'] == star, 'Sub_Alfv_uphase'])[0]
            KS = np.array(targets.loc[targets['Host_ID'] == star, 'KS'])[0]
            AD = np.array(targets.loc[targets['Host_ID'] == star, 'AD'])[0]
            # if np.isnan(lower_phase) == True or lower_phase == 0 or period > 300:
            #     continue
            files = os.listdir("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/Final_Data/" + str(stars))
            c = color[c_index]
            time = np.linspace(1250, 4000, num=100000)
            phases = np.mod((np.linspace(2457000 + 1250, 2457000 + 4000, num=100000) - (periastron_epoch + period/2)), period)/period
            check = 0
            # if peri == True:
            #     for index, phase in enumerate(phases):
            #         if phase > lower_phase and phase < upper_phase and check == 0:
            #             lower = index
            #             check = 1
            #         if check == 1 and phase > upper_phase:
            #             upper = index
            #             ax.axvspan(time[lower], time[upper], 0 + shift2, 1/float(num) + shift2, color = c, alpha = 0.1)
            #             check = 0
            
            for data in files:
                try:
                    data_dir = "C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/Final_Data/" + str(stars) + '/' + str(data)
                    
                    time, flux, error = Flare.TESS_data_extract(data_dir, PDCSAP_ERR=True)
                    time, flux, error = Flare.delete_nans(time, flux, error)
                    if np.min(time) < min_time:
                        min_time = np.min(time)
                    if np.max(time) > max_time:
                        max_time = np.max(time)
                    ax.scatter(time, (((0.04)*(flux-np.min(flux))/(np.max(flux)-np.min(flux))) + (shift+.955)), s=0.5, c = 'black')
                except:
                    continue
            ax.hlines(1 + shift, xmin =1250, xmax=4000, color = 'black')
            if peri == True:
                for epochs in flare_epochs:
                    
                    phase = np.mod((epochs + 2457000) - (periastron_epoch + period/2), period)/period
                    lower_phase = 0.5 - 1.25*np.abs(upper_phase - lower_phase)/2 
                    upper_phase = 0.5 + 1.25*np.abs(upper_phase - lower_phase)/2 
                    if phase < upper_phase and phase > lower_phase:
                        ax.vlines(epochs, ymin=(0.952 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--', color='cyan')
                    else:
                        ax.vlines(epochs, ymin=(0.952 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--', color='red')
                ax.text(3280, .96 + shift,str('Periastron Candidate ') + str(counter3) +'\n' + str(len(flare_epochs)) +' Flares' + '\n' + 'KS: ' + str(round(KS, 3)) + ' AD: ' + str(round(AD, 3)), fontdict=font_small)
                counter3 -= 1
                shift += 0.05
                c_index += 1
                shift2 += 1/float(num)
            elif peri == False:
                for epochs in flare_epochs2:
                    ax.vlines(epochs, ymin=(0.952 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--', color='red')
                for epochs in flare_epochs:
                    ax.vlines(epochs, ymin=(0.952 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--', color='red')
                if counter < 2 and star != 'TOI-833':
                    ax.text(3280, .96 + shift, str('Transit Candidate ') + str(counter2) +'\n' + str(len(flare_epochs2)+len(flare_epochs)) +' Flares' + '\n' + 'KS: ' + str(round(KS, 4)) + ' AD: ' + str(round(AD, 3)), fontdict=font_small)
                    counter2 -= 1
                elif counter >= 2 or star == 'TOI-833':
                    ax.text(3280, .96 + shift, str('Transit Candidate ') + str(counter2) +'\n' + str(len(flare_epochs2)+len(flare_epochs)) +' Flares' + '\n' + 'KS: ' + str(round(KS, 4)) + ' AD: ' + str(round(AD, 3)), fontdict=font_small)
                    counter2 -= 1
                shift += 0.05
                c_index += 1
                shift2 += 1/float(num)
        count += 1
    counter += 1
ax.hlines(1, 1500, 1500, linestyle = '--', color = 'red', label = 'Flare Epoch')
ax.hlines(1, 1500, 1500, linestyle = '--', color = 'cyan', label = 'Sub-Alfvenic Flare')
ax.set_xticklabels(ax.get_xticklabels(), fontdict = font_small)
ax.set_ylim(0.95, 0.95+ shift)
ax.set_xlim(1250, 4000)
ax.legend(loc = 'lower left', prop={"family":"serif"})
plt.savefig('CandidateLCs.png', dpi=1000, bbox_inches='tight' )    
plt.show()