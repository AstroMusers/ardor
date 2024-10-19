# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 21:47:16 2023

@author: Nate Whitsett
"""

from astropy.io import fits
from astropy.timeseries import LombScargle as LS
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import numpy as np
import pandas as pd
import statistics as st
import os
import Flare
time, flux = Flare.TESS_data_extract('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/G Type TESS Data/G 2.24312e+18/mastDownload/TESS/tess2019226182529-s0015-0000000236965747-0151-s/tess2019226182529-s0015-0000000236965747-0151-s_lc.fits')

TESS_Folder_ID = [x[1] for x in os.walk('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/M Dwarf Hosts/')]
total_flares = 0
total_possible_flares = 0
total_observation_time = 0
flare_lengths = []

TOI_ID_list = []
flare_number = []
peak_time = []
amplitude = []
time_scale = []
Teff = []
radius = []
flare_phase = []
total_flare_energies = []
total_flare_phases = []
item_count = 0
# for M_dwarves in TESS_Folder_ID[0]:
#     if M_dwarves.endswith('.01') == True:
        
#         ##Iteration Scheme
#         TOI_ID = float(M_dwarves[3::])
#         a = os.listdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/M Dwarf Hosts/' + M_dwarves)
#         # os.mkdir('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flares_New/' + M_dwarves + '/')
#         print(item_count, M_dwarves)
        
#         ##Relevant parameters from TOI catalog
#         flare_count = 1
        
#         ##Trackable values per star
phase_folded_flare_list = []
flare_amplitude = []
flare_time_scale = []
flare_energy = []
accepted_flare_index = []
flare_phase = []
accepted_flare_number = []
observation_time = 0
possible_flares = 0
#         for folders in a:
b, pdcsap_flux, pdcsap_error = Flare.TESS_data_extract('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/G Type TESS Data/G 1.8084e+18/mastDownload/TESS/tess2022190063128-s0054-0000000050880583-0227-s/tess2022190063128-s0054-0000000050880583-0227-s_lc.fits', PDCSAP_ERR=True)
time, flux = Flare.delete_nans(b, pdcsap_flux)
detrend_flux = Flare.SMA_detrend(time, flux, 80, LS_Iterations=5)
flares, lengths = Flare.flare_ID(detrend_flux, 3)
plt.scatter(time, flux, s=1)
# if folders.endswith('a_fast-lc.fits') == True:
#     observation_time += len(time)*(0.33333333)
#     total_observation_time += len(time)*(0.33333333)
# elif folders.endswith('a_fast-lc.fits') == False:  
observation_time += len(time)*2
total_observation_time += len(time)*2
index = 0
possible_flares += len(flares)
total_possible_flares += len(flares)
flare_count = 0
for flare_events in flares:
    if flare_events >= 100 and len(flux) - flare_events > 100:
        new_time = time[flare_events-100:flare_events+100]
        new_data = flux[flare_events-100:flare_events+100]
    elif flare_events < 100:
        new_time = time[0+flare_events:flare_events+100]
        new_data = flux[0+flare_events:flare_events+100]
    elif len(flux) - flare_events < 100:
        new_time = time[flare_events:]
        new_data = flux[flare_events:]
    recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
    c, d = Flare.flare_ID(np.array(new_data), 3)
    norm_time = time[flare_events]
    events = np.where(new_data == recenter)[0][0]
    criteria1 = False
    if recenter > np.mean(new_data)+3*(np.std(new_data)):
        criteria1 = True
    if criteria1 == True and new_data[events+1] > np.mean(new_data)+2*(np.std(new_data)) and len(c) > 0:
        new_time = (new_time - new_time[events])*24*60
        if lengths[index] >= 25:
            # new_time = np.array(new_time[events-10:events+50])*24*60
            # new_data = np.array(new_data[events-10:events+50])
            alles_data = new_data/np.median(new_data)
            BJD_time = time[events-10:events+50]
            popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
            squares = (alles_data[events:events+30] - Flare.exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
            chi2_cutoff = 18
            plt.clf()
        elif lengths[index] >= 15 and lengths[index] < 25:
            # new_time = np.array(new_time[events-10:events+30])*24*60
            # new_data = np.array(new_data[events-10:events+30])
            alles_data = new_data/np.median(new_data)
            BJD_time = time[events-10:events+30]
            popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
            squares = (alles_data[events:events+20] - Flare.exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
            chi2_cutoff = 9.5
        elif lengths[index] > 5 and lengths[index] < 15:
            # new_time = np.array(new_time[events-10:events+20])*24*60
            # new_data = np.array(new_data[events-10:events+20])
            alles_data = new_data/np.median(new_data)
            BJD_time = time[events-10:events+20]
            popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
            squares = (alles_data[events:events+10] - Flare.exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
            chi2_cutoff = 2.167
        elif lengths[index] <= 5:
            # new_time = np.array(new_time[events-10:events+14])*24*60
            # new_data = np.array(new_data[(events-10):(events+14)])
            alles_data = new_data/np.median(new_data)
            BJD_time = time[events-10:events+14]
            popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
            squares = (alles_data[events:events+7] - Flare.exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
            chi2_cutoff = 1.2
            plt.clf()
        chi_squared = np.sum(squares)
        if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
            half_max = (alles_data[8:15].max()-np.median(alles_data[0:8]))/2
            time_scale.append(popt[1])
            flare_time_scale.append(popt[1])
            amplitude.append(popt[0])
            flare_amplitude.append(popt[0])
            peak_time.append(norm_time)           
            flare_number.append(flare_count)
            X = np.column_stack((new_time[events-30:events+40], alles_data[events-30:events+40]))
            baseline = st.median(new_data)*(lengths[index])*2
            median = st.median(new_data)
            accepted_flare_index.append(flares[index])
            accepted_flare_number.append(flare_count)
            print(flare_count)
            if lengths[index] > 5:
                print('Flare ' + str(flare_count) + ' length: ' + str(lengths[index]))
            # np.savetxt('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flares_New/' + M_dwarves + '/Flare' + str(flare_count) + '.csv', X, delimiter=',')
            flare_count += 1
            total_flares += 1
            
                
    index += 1
print(total_flares)