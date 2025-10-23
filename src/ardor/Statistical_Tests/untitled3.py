# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 16:54:28 2025

@author: whitsett.n
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 13:22:11 2025

@author: whitsett.n
"""

import numpy as np
import matplotlib.cm
import pandas as pd
from ardor.Utils.Utils import df_return, find_absolute_minimum_2d, find_absolute_maximum_2d
from ardor.Statistical_Tests.K_Tests import K_Tests
from ardor.Statistical_Tests.MLE import VM_Unbinned_likelihood
from matplotlib import pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
for fontpath in font_manager.findSystemFonts():
    if 'lmroman10-regular' in fontpath.lower():
        path = fontpath
    if 'lmroman10-italic' in fontpath.lower():
        italicpath = fontpath

# Register and get name
font_manager.fontManager.addfont(path)
font_manager.fontManager.addfont(italicpath)

font = FontProperties(fname=path)
font_name = font.get_name()

rcParams["font.family"] = font_name
rcParams["mathtext.fontset"] = "cm"

colors = matplotlib.cm.tab20(range(20))
data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")
Exo_data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Stat_Tests/Exo_VM.csv")
TOI_data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_News.csv")
TOI_K = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/Stat_Tests/TOI_VM.csv")
hosts = set(data['Host_ID'])
TOIs = set(data['Host_ID'])
TOI_pku = TOI_K['p_KU_TOI'].to_list()
TOI_TS = TOI_K['TS_TOI'].to_list()
host_pku = Exo_data['p_KU_Trans'].to_list()
peri_pku = Exo_data['p_KU_Peri'].to_list()
host_TS = Exo_data['TS_Trans'].to_list()
peri_TS = Exo_data['TS_Peri'].to_list()
data_pku = host_pku + peri_pku + TOI_pku
data_TS = host_TS + peri_TS + TOI_TS
base = []
twosig = 0
threesig = 0
foursig = 0
fivesig = 0
total = 0
for datas in data_pku:
    if datas < 0.05 and datas > 0.003:
        twosig += 1
    if datas < 0.003 and datas > 0.00006:
        threesig += 1
    if datas < 0.00006 and datas > 0.0000003:
        foursig += 1
    if datas < 0.0000003:
        fivesig += 1
    total += 1
print(twosig/total, threesig/total, foursig/total, fivesig/total)
for host in hosts:
    flares = data.loc[data['Host_ID'] == host, 'Flare_Epoch']
    N = len(flares)
    if N >= 3:
        base.append(N)
for TOI in TOIs:
    flares = data.loc[data['Host_ID']== TOI, 'Flare_Epoch']
    N = len(flares)
    if N >= 3:
        base.append(N)
p_KU_grand = []
TS_grand = []
for trials in range(100):
    TS_list = []
    p_KU = []
    for samples in base[0:len(hosts)]:
        flare_samp = np.random.uniform(0, 1, size = samples)
        result = K_Tests(flare_samp, [1], [2457000], output_message=False, sampling=False)
        mu, kappa, TS = VM_Unbinned_likelihood(flare_samp)
        TS_list.append(TS)
        p_KU.append(result[1])
    print(trials)
    p_KU_grand.append(p_KU)
    TS_grand.append(TS_list)
fig, (ax1, ax2) = plt.subplots(1,2, sharey=True, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(7, 3.5))
twosig = 0
threesig = 0
foursig = 0
fivesig = 0
total = 0
log_bins = np.logspace(np.log10(np.nanmin(data_pku)), np.log10(find_absolute_maximum_2d(p_KU_grand)), 20)
ax1.hist(data_pku, bins=log_bins, histtype='step', color='green', label='Hosts')
ax1.hist(TOI_pku, bins=log_bins, histtype='step', color='red', label='TOI PC')
for data_sets in p_KU_grand:
    
    log_bins = np.logspace(np.log10(find_absolute_minimum_2d(p_KU_grand)), np.log10(find_absolute_maximum_2d(p_KU_grand)), 20)
    ax1.hist(data_sets, bins=log_bins, alpha=0.15, log = True,color='black')
    for data in data_sets:
        if data < 0.05 and data > 0.003:
            twosig += 1
        if data < 0.003 and data > 0.00006:
            threesig += 1
        if data < 0.00006 and data > 0.0000003:
            foursig += 1
        if data < 0.0000003:
            fivesig += 1
        total += 1
    print(twosig, threesig, foursig, fivesig, total)
print(twosig, threesig, foursig)
ax1.set_ylim(0.5, 500)
ax1.axvline(0.05,linestyle = '--', color = colors[0], label = r'$2\sigma$')
ax1.axvline(0.003,linestyle = '--', color = colors[1], label = r'$3\sigma$')
ax1.axvline(0.00006,linestyle = '--', color = colors[2], label = r'$4\sigma$')
ax1.legend(loc = 'upper left', ncol=2, handlelength=1, columnspacing=1)
ax1.set_xlabel(r'$p_{KU}$')
ax1.set_ylabel('Count')
ax1.set_xscale('log')
twosig = 0
threesig = 0
foursig = 0
fivesig = 0
total = 0
log_bins = np.linspace(find_absolute_minimum_2d(TS_grand), 5.1, 15)
ax2.set_xlim(-0.1, 5.25)
ax2.hist(data_TS, bins=log_bins, histtype='step', color='green', label='Hosts')
ax2.hist(TOI_TS, bins=log_bins, histtype='step', color='red', label='TOI PC')
for data_sets in TS_grand:
    log_bins = np.linspace(find_absolute_minimum_2d(TS_grand), find_absolute_maximum_2d(TS_grand), 15)
    ax2.hist(data_sets, bins=log_bins, alpha=0.15, log = True, color='black')
    for data in data_sets:
        if data > 2 and data < 3:
            twosig += 1
        if data > 3 and data < 4:
            threesig += 1
        if data > 4 and data < 5:
            foursig += 1
        if data > 5:
            
            fivesig += 1
        total += 1
    print(twosig, threesig, foursig, fivesig, total)
ax2.set_ylim(0.5, 500)
ax2.axvline(2,linestyle = '--', color = colors[0], label = r'$2\sigma$')
ax2.axvline(3,linestyle = '--', color = colors[1], label = r'$3\sigma$')
ax2.axvline(4,linestyle = '--', color = colors[2], label = r'$4\sigma$')
ax2.axvline(5, linestyle = '--', color = colors[3], label = r'$5\sigma$')
ax2.set_xlabel(r'$\sqrt{TS}$')
plt.subplots_adjust(wspace=0,hspace=0)
plt.tight_layout()
plt.savefig("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Induced_Flares/TS_K_Samp.png", dpi=400)
plt.show()