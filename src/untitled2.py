# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 13:22:11 2025

@author: whitsett.n
"""

import numpy as np
import pandas as pd
from ardor.Utils.Utils import df_return, find_absolute_minimum_2d, find_absolute_maximum_2d
from ardor.Statistical_Tests.K_Tests import K_Tests
from ardor.Statistical_Tests.MLE import VM_Unbinned_likelihood
from matplotlib import pyplot as plt
data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")
Exo_data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Stat_Tests/Exo_VM.csv")
hosts = set(data['Host_ID'])
host_pku = Exo_data['p_KU_Trans'].to_list()
peri_pku = Exo_data['p_KU_Peri'].to_list()
host_TS = Exo_data['TS_Trans'].to_list()
peri_TS = Exo_data['TS_Peri'].to_list()
data_pku = host_pku + peri_pku
data_TS = host_TS + peri_TS
base = []
twosig = 0
threesig = 0
foursig = 0
fivesig = 0
total = 0
for data in data_pku:
    if data < 0.05 and data > 0.003:
        twosig += 1
    if data < 0.003 and data > 0.00006:
        threesig += 1
    if data < 0.00006 and data > 0.0000003:
        foursig += 1
    if data < 0.0000003:
        fivesig += 1
    total += 1
print(twosig/total, threesig/total, foursig/total, fivesig/total)
# for host in hosts:
#     flares = data.loc[data['Host_ID'] == host, 'Flare_Epoch']
#     N = len(flares)
#     if N >= 3:
#         base.append(N)

# p_KU_grand = []
# TS_grand = []
# for trials in range(100):
#     TS_list = []
#     p_KU = []
#     for samples in base:
#         flare_samp = np.random.uniform(0, 1, size = samples)
#         result = K_Tests(flare_samp, [1], [2457000], output_message=False, sampling=False)
#         mu, kappa, TS = VM_Unbinned_likelihood(flare_samp)
#         TS_list.append(TS)
#         p_KU.append(result[1])
#     p_KU_grand.append(p_KU)
#     TS_grand.append(TS_list)
# fig, (ax1, ax2) = plt.subplots(1,2, sharey=True, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(7, 3.5))
# twosig = 0
# threesig = 0
# foursig = 0
# fivesig = 0
# total = 0
# log_bins = np.logspace(np.log10(np.nanmin(data_pku)), np.log10(find_absolute_maximum_2d(p_KU_grand)), 20)
# ax1.hist(data_pku, bins=log_bins, histtype='step', color='black', label='Observed')
# for data_sets in p_KU_grand:
    
#     log_bins = np.logspace(np.log10(find_absolute_minimum_2d(p_KU_grand)), np.log10(find_absolute_maximum_2d(p_KU_grand)), 20)
#     ax1.hist(data_sets, bins=log_bins, alpha=0.05, log = True)
#     for data in data_sets:
#         if data < 0.05 and data > 0.003:
#             twosig += 1
#         if data < 0.003 and data > 0.00006:
#             threesig += 1
#         if data < 0.00006 and data > 0.0000003:
#             foursig += 1
#         if data < 0.0000003:
#             fivesig += 1
#         total += 1
#     print(twosig, threesig, foursig, fivesig, total)
# print(twosig, threesig, foursig)
# ax1.axvline(0.68, linestyle='--', label=r'$1\sigma$', color='red')
# ax1.axvline(0.05, linestyle='--',label=r'$2\sigma$', color='blue')
# ax1.axvline(0.003, linestyle='--', label=r'$3\sigma$', color='magenta')
# ax1.axvline(0.00006, linestyle='--',label=r'$4\sigma$', color='cyan')
# ax1.legend(loc = 'upper left', ncol=2, handlelength=1, columnspacing=1)
# ax1.set_xlabel(r'$p_{KU}$')
# ax1.set_ylabel('Count')
# ax1.set_xscale('log')
# twosig = 0
# threesig = 0
# foursig = 0
# fivesig = 0
# total = 0
# log_bins = np.linspace(find_absolute_minimum_2d(TS_grand), 5.1, 15)
# ax2.set_xlim(-0.1, 5.25)
# ax2.hist(data_TS, bins=log_bins, histtype='step', color='black', label='Observed')
# for data_sets in TS_grand:
#     log_bins = np.linspace(find_absolute_minimum_2d(TS_grand), find_absolute_maximum_2d(TS_grand), 15)
#     ax2.hist(data_sets, bins=log_bins, alpha=0.05, log = True)
#     for data in data_sets:
#         if data > 2 and data < 3:
#             twosig += 1
#         if data > 3 and data < 4:
#             threesig += 1
#         if data > 4 and data < 5:
#             foursig += 1
#         if data > 5:
#             fivesig += 1
#         total += 1
#     print(twosig, threesig, foursig, fivesig, total)
# ax2.axvline(1, linestyle='--', label=r'$1\sigma$', color='red')
# ax2.axvline(2, linestyle='--',label=r'$2\sigma$', color='blue')
# ax2.axvline(3, linestyle='--', label=r'$3\sigma$', color='magenta')
# ax2.axvline(4, linestyle='--',label=r'$4\sigma$', color='cyan')
# ax2.set_xlabel(r'$\sqrt{TS}$')
# plt.subplots_adjust(wspace=0,hspace=0)
# plt.tight_layout()
# plt.savefig("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Induced_Flares/TS_K_Samp.png", dpi=400)
# plt.show()