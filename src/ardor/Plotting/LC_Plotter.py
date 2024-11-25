# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:08:22 2024

@author: natha
"""

from ardor.Flares.Flare import tier0
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import gca
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties

for fontpath in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    if 'lmroman10-regular'.lower() in fontpath.lower():
        path = fontpath
font = FontProperties(fname=path)
font_label = FontProperties(fname=path)
font_small = FontProperties(fname=path)
targets = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/KS_KU_Tests/K_Hosts.csv")
flares = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")
flares_TOI = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares.csv")



hosts = targets["Host_ID"].tolist()
directory = os.listdir("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Hosts")
print(directory)
directory = ['TOI-1062','AUMic', 'HD153557', 'HD191939', 'TOI-833', 'TOI-540']
proper = ['TOI-1062','AU Mic', 'HD 153557', 'HD 191939', 'TOI-833', 'TOI-540']
shift = 0
shift2 = 0
color = ('#d7191c',
'#fdae61',
'#0571b0',
'#7b3294',
'#1a9641','#d01c8b', '#018571', '#a6611a', '#a50026', '#006837')
fig, ax = plt.subplots()
fig.set_size_inches(7, 8)
ax.set_xlabel('BJD - 2450000', font=font, fontsize = 12)

ax.get_yaxis().set_ticks([])
ax.set_ylabel('Rel. Flux', font=font, fontsize = 12)
min_time = 2000
max_time = 0
c_index = 0
counter2 = 3
counter3 = 3
num = 6
count = 0
for index, stars in enumerate(directory):
    star = stars.replace(' ', '')
    counter = 1
    if star in hosts:
        print(1)
        if count < 7:
            
            # lower_phase = np.array(targets.loc[targets['Host_ID'] == star, 'Sub_Alfv_lphase'])[0]
            period = np.array(flares.loc[flares['Host_ID'] == star, 'Period'])[0]
            periastron_epoch = np.array(flares.loc[flares['Host_ID'] == star, 'Periastron_Epoch'])[0]
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
            KS = np.array(targets.loc[targets['Host_ID'] == star, 'Peri_KS'])[0]
            KU = np.array(targets.loc[targets['Host_ID'] == star, 'Peri_KU'])[0]
            # if np.isnan(lower_phase) == True or lower_phase == 0 or period > 300:
            #     continue
            files = os.listdir("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/Final_Data/" + str(stars))
            c = color[c_index]
            time = np.linspace(1250, 4000, num=100000)
            phases = np.mod((np.linspace(2457000 + 1250, 2457000 + 4000, num=100000) - (periastron_epoch + period/2)), period)/period
            check = 0        
            for data in files:
                try:
                    data_dir = "C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/Final_Data/" + str(stars) + '/' + str(data)
                    lc = tier0(data_dir)
                    if np.min(time) < min_time:
                        min_time = np.min(time)
                    if np.max(time) > max_time:
                        max_time = np.max(time)
                    ax.scatter(lc.time, (((0.04)*(lc.detrended_flux-np.min(lc.detrended_flux))/(np.max(lc.detrended_flux)-np.min(lc.detrended_flux))) + shift+.955), s=0.5, c = 'black')
                except:
                    continue
            ax.hlines(1 + shift, xmin =1250, xmax=4500, color = 'black')
            for epochs in flare_epochs:
                ax.vlines(epochs, ymin=(0.952 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--', color='red')
            if np.isnan(KS) == False:
                ax.text(3500, .96 + shift,str(proper[index]) +'\n' + str(len(flare_epochs)) +' Flares' + '\n' + 'KU: ' + str("{:.2E}".format(round(KU, 3))) + ' KS: ' + str("{:.2E}".format(round(KS, 4))), font=font,fontsize = 14)
            elif np.isnan(KS) == True:
                ax.text(3500, .96 + shift,str(proper[index]) +'\n' + str(len(flare_epochs)) +' Flares' + '\n' + 'KU: ' + str("{:.2E}".format(round(KU, 5))), font=font,fontsize = 14)
            
            counter3 -= 1
            shift += 0.05
            c_index += 1
            shift2 += 1/float(num)
        count += 1
    counter += 1
ax.hlines(1, 1500, 1500, linestyle = '--', color = 'red', label = 'Flare Epoch')
ax.set_xticklabels(ax.get_xticklabels(), font=font,fontsize = 12)
ax.set_ylim(0.95, 0.95+ shift)
ax.set_xlim(1250, 4500)
ax.legend(loc = 'lower left', prop={"family":"serif", "size":12})
plt.savefig('CandidateLCs.png', dpi=1000, bbox_inches='tight' )    
plt.show()