from astropy.io import fits
from astropy.timeseries import LombScargle as LS
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import numpy as np
import pandas as pd
import statistics as st
import os
import warnings
import planck_law as pl
import aflare
import time as timer
import copy
import allesfitter_priors
import shutil
from scipy.interpolate import splrep

# fig = plt.figure(figsize=(15,10))
# gs = GridSpec(2,2, figure = fig)
# ax1 = fig.add_subplot(gs[0,:])
# ax2 = fig.add_subplot(gs[1, 0])
# ax3 = fig.add_subplot(gs[1,1])
# TESS_Folder_ID = [x[1] for x in os.walk('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/Periastron_Hosts/')]
# TOI_Catalog = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Periastron_hosts.csv')
# time, flux, error = TESS_data_extract('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/TESS Data/Periastron_Hosts/Proxima Cen/tess2019140104343-s0012-0000000388857263-0144-s_lc.fits', PDCSAP_ERR=True)
# time, flux = delete_nans(time, flux)
# detrend_flux, mov_average = SMA_detrend(time, flux, 8)
# flares, lengths = flare_ID(detrend_flux, 3)
# ax1.plot(time, np.array(flux)- 145000, linewidth=2)
# ax1.plot(time, detrend_flux, linewidth=2, c='green')
# ax1.plot(time, np.array(mov_average) - 145000, c='red', linestyle='--', linewidth=1.75)
# ax1.set_xlim(1632.1, 1633)
# ax1.set_ylim(-500, 3500)
# print(len(flares))
# for index in range(len(flares)):
#     ax1.scatter(time[flares[index]:flares[index]+lengths[index]-1], np.array(detrend_flux[flares[index]:flares[index]+lengths[index]-1]), s=15, c='r')
# # plt.axhline(y=np.median(detrend_flux), c = 'r')
# # ax2.plot(time[4300:4400], detrend_flux[4300:4400], linewidth=1, c ='green')

# count = 0
# for flare_events in flares[40], flares[48]:
#     if flare_events >= 100 and len(flux) - flare_events > 100:
#         new_time = time[flare_events-100:flare_events+100]
#         new_data = flux[flare_events-100:flare_events+100]
#     elif flare_events < 100:
#         new_time = time[0+flare_events:flare_events+100]
#         new_data = flux[0+flare_events:flare_events+100]
#     elif len(flux) - flare_events < 100:
#         new_time = time[flare_events:]
#         new_data = flux[flare_events:]
#     recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
#     norm_time = time[flare_events]
#     events = np.where(new_data == recenter)[0][0]    
#     criteria1 = False
#     # if recenter > np.mean(new_data)+3*(np.std(new_data)):
#         # criteria1 = True
#     # if criteria1 == True and new_data[events+1] > np.mean(new_data)+2*(np.std(new_data)):
#     new_time = (new_time - new_time[events])*24*60
#     if lengths[index] >= 25:
#         alles_data = new_data/np.median(new_data)
#         popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
#         squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
#         chi2_cutoff = 20.843
#     elif lengths[index] >= 15 and lengths[index] < 25:
#         alles_data = new_data/np.median(new_data)
#         popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
#         squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
#         chi2_cutoff = 11.912
#     elif lengths[index] > 5 and lengths[index] < 15:
#         alles_data = new_data/np.median(new_data)
#         popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
#         squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
#         chi2_cutoff = 3.455
#     elif lengths[index] <= 5:
#         alles_data = new_data/np.median(new_data)
#         popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
#         squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
#         chi2_cutoff = 1.3
#     chi_squared = np.sum(squares)
#     if count == 0:
#         ax2.plot(new_time[events:events+30], exp_decay(new_time[events:events+30], popt[0], popt[1], popt[2]), linewidth=2, linestyle='--', c='r')
#         ax2.plot(new_time[70:150], alles_data[70:150], c = 'g')
#         print(chi_squared)
#     if count == 1:
#         ax3.plot(new_time[events:events+10], exp_decay(new_time[events:events+10], popt[0], popt[1], popt[2]), linewidth=2, linestyle='--', c='r')
#         ax3.plot(new_time[70:150], alles_data[70:150], c = 'g')
#         print(chi_squared)
#     count += 1
# plt.savefig('')
# plt.ylim(40, 100)


