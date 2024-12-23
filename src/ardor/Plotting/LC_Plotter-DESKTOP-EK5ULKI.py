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
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib.colors import Normalize
for fontpath in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    if 'lmroman10-regular'.lower() in fontpath.lower():
        path = fontpath
font = FontProperties(fname=path)
font_label = FontProperties(fname=path)
font_small = FontProperties(fname=path)
rcParams["mathtext.fontset"] = "cm"
targets = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/KS_KU_Tests/K_Tests/Exo/K_All_Exo.csv")
targets_TOI = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/KS_KU_Tests/K_Tests/TOI/K_All_TOI.csv")
flares = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")
flares_TOI = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares.csv")
directory = os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs')

targets.sort_values(by='p_KU', inplace=True, ascending = False)
hosts = targets["Host_ID"].tolist()
hosts_TOI = targets_TOI["Host_ID"].tolist()
shift = 0
shift2 = 0
fig, ax = plt.subplots()
fig.set_size_inches(10.5, 7.5)
ax.set_xlabel('BJD - 2450000', font=font, fontsize = 12)

ax.get_yaxis().set_ticks([])
ax.set_ylabel('Relative Flux', font=font, fontsize = 12)
min_time = 2000
max_time = 0
c_index = 0
counter2 = 3
counter3 = 3
num = 6
count = 0
norm = Normalize(vmin=0.5, vmax=1)
cmap = 'plasma'
m = cm.ScalarMappable(norm=norm, cmap=cmap)













##################### HOST LOOP ###############################################
for index, star in enumerate(hosts):
    counter = 1
    stars = star.replace(' ', '')
    
    period = np.array(flares.loc[flares['Host_ID'] == stars, 'Period'])[0]
    periastron_epoch = np.array(flares.loc[flares['Host_ID'] == stars, 'Periastron_Epoch'])[0]
    flare_epochs = np.array(flares.loc[flares['Host_ID'] == stars, 'Flare_Epoch_TESS'])
    peri_bool = np.array(targets.loc[targets['Host_ID'] == star, 'Periastron?'])
    peri = peri_bool.item()
    peri_upper = np.array(flares.loc[flares['Host_ID'] == stars, 'Periastron_Upper'])
    peri_upper = peri_upper[0]
    peri_lower = np.array(flares.loc[flares['Host_ID'] == stars, 'Periastron_Lower'])
    peri_lower = peri_lower[0]
    dlogZ = np.array(targets.loc[targets['Host_ID'] == star, 'dlog(Z)'])
    dlogZ = dlogZ.item()
    KS = np.array(targets.loc[targets['Host_ID'] == star, 'p_KS'])[0]
    KS = KS.item()
    KU = np.array(targets.loc[targets['Host_ID'] == star, 'p_KU'])[0]
    KU = KU.item()
    AD = np.array(targets.loc[targets['Host_ID'] == star, 'p_AD'])[0]
    AD = AD.item()
    K = [KS, KU, AD]
    K_val = np.argmin(np.array(K))
    K_str = ['KS', 'KU', 'AD']
    if (K[K_val] < 5e-2 and dlogZ > 5) and len(flare_epochs) > 5:
        if peri == False:
            try:
                flare_epochs2 = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == int(stars), 'Flare_Epoch'])
            except:
                flare_epochs2 = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == stars, 'Flare_Epoch'])
            flare_epochs = np.array(flares.loc[flares['Host_ID'] == stars, 'Flare_Epoch_TESS'])
            flare_phases = np.array(flares.loc[flares['Host_ID'] == stars, 'Transit_Phase'])
        elif peri == True:
            continue
        files = os.listdir("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Hosts/" + str(stars))
        time = np.linspace(1250, 4000, num=100000)
        phases = np.mod((np.linspace(2457000 + 1250, 2457000 + 4000, num=100000) - (periastron_epoch + period/2)), period)/period
        check = 0        
        for data in files:
            try:
                data_dir = "C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Hosts/" + str(stars) + '/' + str(data)
                lc = tier0(data_dir)
                if np.min(time) < min_time:
                    min_time = np.min(time)
                if np.max(time) > max_time:
                    max_time = np.max(time)
                ax.scatter(lc.time, (((0.04)*(lc.detrended_flux-np.min(lc.detrended_flux))/(np.max(lc.detrended_flux)-np.min(lc.detrended_flux))) + shift+.955), s=0.5, c = 'black')
            except:
                continue
        ax.hlines(1 + shift, xmin =1250, xmax=4500, color = 'black')
        
        for index, epochs in enumerate(flare_epochs):
            # if (flare_phases[index] < peri_upper + 0.5 and flare_phases[index] > peri_lower + 0.5) or (flare_phases[index] < 0.6 and flare_phases[index] < 0.4):
            #     ax.vlines(epochs, colors = 'red', ymin=(0.95 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--')
            ax.vlines(epochs, colors = 'blue', ymin=(0.95 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--')

        # ax.text(3500, .9625 + shift,str(proper[index]) +'\n' + str(len(flare_epochs)) +' Flares' + '\n' + 'KU: ' + str("{:.2E}".format(round(KU, 3))) + ' KS: ' + str("{:.2E}".format(round(KS, 4))), font=font,fontsize = 12)
        if star == 'TOI-1062':
            ax.text(3600, .9575 + shift,str(star) +'\n$N_{flare} = ' + str(len(flare_epochs)) + '$\n$p_{' + K_str[K_val] + '}:\,' + str(round(K[K_val],6)) + '$\n' + "$\overline{\Delta\log(Z)}=$" + str(round(dlogZ,1)) , font=font,fontsize =9, horizontalalignment='left')
        else:
            ax.text(3600, .9575 + shift,str(star) +'\n$N_{flare} = ' + str(len(flare_epochs)) + '$\n$p_{' + K_str[K_val] + '}:\,' + str(round(K[K_val], 3)) + '$\n' + "$\overline{\Delta\log(Z)}=$" + str(round(dlogZ,1)) , font=font,fontsize =9, horizontalalignment='left')
        
        counter3 -= 1
        shift += 0.05
        c_index += 1
        shift2 += 1/float(num)
    counter += 1

ax.hlines(1, 1500, 1500, linestyle = '--', color = 'blue', label = 'Flares')
# ax.hlines(1, 1500, 1500, linestyle = '--', color = 'red', label = 'Flares Near Periastron')
ax.set_xticklabels(ax.get_xticklabels(), font=font,fontsize = 11)
ax.set_ylim(0.95, 0.95+ shift)
ax.set_xlim(1250, 4000)
ax.legend(loc = 'lower left', prop={"family":"serif", "size":11})
plt.savefig('CandidateLCs_Peri.png', dpi=300, bbox_inches='tight' )    
plt.show()


################## TOI LOOP ###################################
# for index, star in enumerate(hosts_TOI):
#     if star in hosts_TOI:
#         try:
#             counter = 1
#             stars = star
#             period = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == float(stars), 'Period'])[0]
#             # periastron_epoch = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == stars, 'Transit_Epoch'])[0]
#             flare_epochs = np.array(flares.loc[flares_TOI['Host_ID'] == float(stars), 'Flare_Epoch_TESS'])
#             peri_bool = np.array(targets.loc[targets_TOI['Host_ID'] == float(stars), 'Periastron?'])
#             peri = False
#             # peri_upper = np.array(flares.loc[flares_TOI['Host_ID'] == stars, 'Periastron_Upper'])
#             # peri_upper = peri_upper[0]
#             # peri_lower = np.array(flares.loc[flares_TOI['Host_ID'] == stars, 'Periastron_Lower'])
#             # peri_lower = peri_lower[0]
#             dlogZ = np.array(targets.loc[targets_TOI['Host_ID'] == float(stars), 'dlog(Z)'])
#             dlogZ = dlogZ.item()
#             KS = np.array(targets.loc[targets_TOI['Host_ID'] == float(stars), 'p_KS'])[0]
#             KS = KS.item()
#             KU = np.array(targets.loc[targets_TOI['Host_ID'] == float(stars), 'p_KU'])[0]
#             KU = KU.item()
#             AD = np.array(targets.loc[targets_TOI['Host_ID'] == float(stars), 'p_AD'])[0]
#             AD = AD.item()
#             K = [KS, KU, AD]
#             K_val = np.argmin(np.array(K))
#             K_str = ['KS', 'KU', 'AD']
#             if (K[K_val] < 5e-2 and dlogZ > 3) and len(flare_epochs) > 5:
#                 if peri == False:
#                     try:
#                         flare_epochs2 = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == int(stars), 'Flare_Epoch'])
#                     except:
#                         flare_epochs2 = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == stars, 'Flare_Epoch'])
#                     flare_epochs = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == stars, 'Flare_Epoch_TESS'])
#                     flare_phases = np.array(flares_TOI.loc[flares_TOI['Host_ID'] == stars, 'Periastron_Phase'])
#                 elif peri == False:
#                     continue
#                 files = os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs/' + str(stars))
#                 time = np.linspace(1250, 4000, num=100000)
#                 # phases = np.mod((np.linspace(2457000 + 1250, 2457000 + 4000, num=100000) - (periastron_epoch + period/2)), period)/period
#                 check = 0        
#                 for data in files:
#                     try:
#                         data_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs/' + str(stars) + '/' + str(data)
#                         lc = tier0(data_dir)
#                         if np.min(time) < min_time:
#                             min_time = np.min(time)
#                         if np.max(time) > max_time:
#                             max_time = np.max(time)
#                         ax.scatter(lc.time, (((0.04)*(lc.detrended_flux-np.min(lc.detrended_flux))/(np.max(lc.detrended_flux)-np.min(lc.detrended_flux))) + shift+.955), s=0.5, c = 'black')
#                     except:
#                         continue
#                 ax.hlines(1 + shift, xmin =1250, xmax=4500, color = 'black')
                
#                 for index, epochs in enumerate(flare_epochs):
#                     # if (flare_phases[index] < peri_upper + 0.5 and flare_phases[index] > peri_lower + 0.5) or (flare_phases[index] < 0.6 and flare_phases[index] < 0.4):
#                     #     ax.vlines(epochs, colors = 'red', ymin=(0.95 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--')
#                     ax.vlines(epochs, colors = 'blue', ymin=(0.95 + shift), ymax = (1 + shift), alpha=0.75, linestyle='--')
        
#                 # ax.text(3500, .9625 + shift,str(proper[index]) +'\n' + str(len(flare_epochs)) +' Flares' + '\n' + 'KU: ' + str("{:.2E}".format(round(KU, 3))) + ' KS: ' + str("{:.2E}".format(round(KS, 4))), font=font,fontsize = 12)
#                 ax.text(3600, .9575 + shift,str(star) +'\n$N_{flare} = ' + str(len(flare_epochs)) + '$\n$p_{' + K_str[K_val] + '}:\,' + str(round(K[K_val], 3)) + '$\n' + "$\overline{\Delta\log(Z)}=$" + str(round(dlogZ,1)) , font=font,fontsize = 9, horizontalalignment='left')
                
#                 counter3 -= 1
#                 shift += 0.05
#                 c_index += 1
#                 shift2 += 1/float(num)
#             counter += 1
#         except:
#             continue

# ax.hlines(1, 1500, 1500, linestyle = '--', color = 'blue', label = 'Flares')
# # ax.hlines(1, 1500, 1500, linestyle = '--', color = 'red', label = 'Flares Near Periastron')
# ax.set_xticklabels(ax.get_xticklabels(), font=font,fontsize = 11)
# ax.set_ylim(0.95, 0.95+ shift)
# ax.set_xlim(1250, 4000)
# ax.legend(loc = 'lower left', prop={"family":"serif", "size":11})
# plt.savefig('CandidateLCs_TOI.png', dpi=300, bbox_inches='tight' )    
# plt.show()