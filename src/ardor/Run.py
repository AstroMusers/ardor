# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 21:13:01 2024

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
import warnings
import planck_law as pl
import aflare
import time as timer
import copy
import allesfitter_priors
import shutil
import Flare
# import allesfitter

warnings.filterwarnings("ignore")
data_dir = '/ugrad/whitsett.n/Flare_Data/Tier_2/Hosts'
output_dir = '/ugrad/whitsett.n/Flare_Data/Tier_2'
TESS_Folder_ID = os.listdir(data_dir)
TOI_Catalog = pd.read_csv('/ugrad/whitsett.n/Flare_Data/Reference_Files/Alfven_Catalog.csv')
total_flares = 0
total_possible_flares = 0
total_observation_time = 0
TOI_ID_list = []
flare_number = []
peak_time = []
amplitude = []
time_scale = []
Teff = []
radius = []
flare_phase = []
TOI_period = []
total_flare_energies = []
total_flare_phases = []
item_count = 0
total_periastron_list = []
total_periastron_epoch_list = []
total_epoch_list = []
list_index = 0
e_list = []
a_list = []
tier0_tau = []
tier1_tau = []
tier2_tau = []

for M_dwarves in TESS_Folder_ID[100:]:
    print(M_dwarves)
    ##Trackable values per star
    phase_folded_flare_list = []
    flare_amplitude = []
    flare_time_scale = []
    flare_energy = []
    accepted_flare_index = []
    flare_phase = []
    accepted_flare_number = []
    observation_time = 0
    possible_flares = 0
    transit_epoch_list = []
    periastron_epoch_list = []
    epoch_list = []
    
    #Total Trackable Lists for all Data
    e_list = []
    a_list = []
    min_approach_list = []
    TOI_ID_list = []
    flare_number = []
    peak_time = []
    amplitude = []
    time_scale = []
    Teff = []
    radius = []
    flare_phase = []
    TOI_period = []
    total_flare_energies = []
    total_flare_phases = []
    item_count = 0
    total_transit_epoch_list = []
    total_periastron_epoch_list = []
    total_epoch_list = []
    
    ##Iteration Scheme
    TOI_ID = str(M_dwarves).replace(' ', '')
    b = os.listdir(data_dir + '/' + M_dwarves)
    try:
        os.mkdir(output_dir + '/' + M_dwarves + '/')
    except:
        print('Already analyzed')
        continue
    ##Relevant parameters from TOI catalog
    period = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'pl_orbper'])[0]
    epoch = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'pl_tranmid'])[0]
    stellar_radius = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'st_rad'])[0]
    T = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'st_teff'])[0]
    e = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'pl_orbeccen'])[0]
    a = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'pl_orbsmax'])[0]
    min_approach = (1-e)*a
    transit_epoch = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'pl_tranmid'])[0]
    epoch_periastron = np.array(TOI_Catalog.loc[TOI_Catalog['Host_Name'] == TOI_ID, 'pl_orbtper'])[0]
    if T == '':
        T = np.NAN
    if stellar_radius == '':
        stellar_radius = np.NAN
    if epoch_periastron == '':
        epoch_periastron = np.NAN
    if transit_epoch == '':
        transit_epoch = np.NAN
    if epoch == '':
        epoch = np.NAN
    flare_count = 0


    for indexs, folders in enumerate(b):
        t1 = timer.time()
        try:
            b, pdcsap_flux, pdcsap_error = Flare.TESS_data_extract(data_dir + '/' + M_dwarves + '/' + folders, PDCSAP_ERR=True)
        except:
            continue
        if folders.endswith('a_fast-lc.fits') == True:
            observation_time += len(pdcsap_flux)*(0.33333333)
            total_observation_time += len(pdcsap_flux)*(0.33333333)
            fast = True
        elif folders.endswith('a_fast-lc.fits') == False:  
            observation_time += len(pdcsap_flux)*2
            total_observation_time += len(pdcsap_flux)*2
            fast = False
        time, flux, pdcsap_error = Flare.delete_nans(b, pdcsap_flux, pdcsap_error)
        detrend_flux = Flare.SMA_detrend(time, flux, pdcsap_error, 80)
        t2 = timer.time()
        tier0_tau.append(t2-t1)
        flares, lengths = Flare.flare_ID(detrend_flux, 3, fast=fast)
        index = 0
        possible_flares += len(flares)
        total_possible_flares += len(flares)
        print('LC #: ' + str(indexs + 1),len(flares))
        for flare_events in flares:         
            if flare_events >= 50 and len(flux) - flare_events > 50:
                new_time = time[flare_events-50:flare_events+50]
                new_data = flux[flare_events-50:flare_events+50]
                new_error = pdcsap_error[flare_events-50:flare_events+50]
                recenter = np.max(new_data[45:50+lengths[index]])
            elif flare_events < 50:
                new_time = time[0:flare_events+50]
                new_data = flux[0:flare_events+50]
                new_error = pdcsap_error[flare_events:flare_events+50]
                try:
                    recenter = np.max(new_data[flare_events-5:flare_events+lengths[index]])
                except:
                    continue
            elif len(flux) - flare_events < 50:
                new_time = time[flare_events:]
                new_data = flux[flare_events:]
                new_error = pdcsap_error[flare_events:]
                try:
                    recenter = np.max(new_data[flare_events-5:flare_events+lengths[index]])
                except:
                    continue
            norm_time = time[flare_events]
            events = np.where(new_data == recenter)[0][0]
            if period != '':
                phase = Flare.phase_folder(new_time, period, epoch)
                flare_phase_value = phase[events]
            elif period == '':
                phase = []
                flare_phase_value = np.NAN
                period = np.NAN
            t3 = timer.time()
            tier1_tau.append(t3-t2)
            new_time = (new_time - new_time[events])*24*60
            try:
                if lengths[index] >= 25:
                    alles_data = new_data/np.median(new_data)
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000, sigma = error[events:events+30], absolute_sigma=True)
                    squares =((alles_data[events:events+30] - Flare.exp_decay(new_time[events:events+30], *popt))/(error[events:events+30]))**2
                    chi_squared = np.sum(squares)/27
                    cutoff = 15
                elif lengths[index] >= 15 and lengths[index] < 25:
                    alles_data = new_data/np.median(new_data)
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+20],alles_data[events:events+20], maxfev=5000, sigma = error[events:events+20], absolute_sigma=True)
                    squares = ((alles_data[events:events+20] - Flare.exp_decay(new_time[events:events+20], *popt))/(error[events:events+20]))**2
                    chi_squared = np.sum(squares)/17
                    cutoff = 7.5
                elif lengths[index] > 5 and lengths[index] < 15:
                    alles_data = new_data/np.median(new_data)
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000, sigma = error[events:events+10], absolute_sigma=True)
                    squares = ((alles_data[events:events+10] - Flare.exp_decay(new_time[events:events+10], *popt))/(error[events:events+10]))**2
                    chi_squared = np.sum(squares)/7
                    cutoff = 5
                elif lengths[index] <= 5:
                    alles_data = new_data/np.median(new_data)
                    error = new_error/np.median(new_data)
                    popt, pcov = curve_fit(Flare.exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000, sigma = error[events:events+7], absolute_sigma=True)
                    squares = ((alles_data[events:events+7] - Flare.exp_decay(new_time[events:events+7], *popt))/(error[events:events+7]))**2
                    chi_squared = np.sum(squares)/4
                    cutoff = 5
                if chi_squared < cutoff and popt[0] > 0 and popt[1] > 0:
                    print(index, lengths[index])
                    plt.plot(new_time[events:events+20], Flare.exp_decay(new_time[events:events+20], *popt))
                    plt.plot(new_time, alles_data)
                    plt.show()
                    a = input()
                    flare_count += 1
                    time_scale.append(popt[1])
                    flare_time_scale.append(popt[1])
                    amplitude.append(popt[0])
                    flare_amplitude.append(popt[0])
                    peak_time.append(norm_time)
                    
                    TOI_ID_list.append(TOI_ID)
                    flare_number.append(flare_count)
                    try:
                        energy = Flare.bolo_flare_energy(popt, T, stellar_radius, pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000), t_unit='days')
                    except:
                        energy = np.NaN
                    flare_energy.append(energy)
                    total_flare_energies.append(energy)
                    Teff.append(T)
                    radius.append(stellar_radius)
                    try:
                        X = np.column_stack((new_time, alles_data, error))
                    except:
                        X = np.column_stack(([0],[0],[0]))
                    baseline = st.median(new_data)*(lengths[index])*2
                    median = st.median(new_data)
                    flare_phase.append(flare_phase_value)
                    accepted_flare_index.append(flares[index])
                    accepted_flare_number.append(flare_count)
                    total_flare_phases.append(flare_phase_value)
                    TOI_period.append(period)
                    total_transit_epoch_list.append(transit_epoch)
                    total_periastron_epoch_list.append(epoch_periastron)
                    transit_epoch_list.append(transit_epoch)
                    periastron_epoch_list.append(epoch_periastron)
                    total_epoch_list.append(epoch)
                    epoch_list.append(epoch)
                    a_list.append(a)
                    e_list.append(e)
                    min_approach_list.append(min_approach)
                    if lengths[index] > 5:
                        print('Flare ' + str(flare_count) + ' length: ' + str(lengths[index]))
                    np.savetxt(output_dir + '/' + M_dwarves + '/Flare' + str(flare_count) + '.csv', X, delimiter=',')
                    total_flares += 1
                    t4 = timer.time()
                    tier2_tau.append(t4-t3)
                    index += 1
                else:
                    print(index, chi_squared, lengths[index])
                    plt.plot(new_time[events:events+20], Flare.exp_decay(new_time[events:events+20], *popt))
                    plt.plot(new_time, alles_data)
                    plt.show()
                    a = input()
                    index += 1
            except TypeError:
                print('Flare at end of LC')
                index += 1
                continue
                
        ZZ = np.column_stack((np.array(np.array(tier0_tau).mean()), np.array(np.array(tier1_tau).mean()/len(flares)), np.array(np.array(tier2_tau).mean())))
        with open("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/Time_Stats/Time_Stats.csv", "a") as f:
            np.savetxt(f, ZZ, delimiter=",", fmt='%s')
            f.close()
        tier0_tau = []
        tier1_tau = []
        tier2_tau = []
    # if flare_count == 0:
    #     os.rmdir(output_dir + '/' + M_dwarves + '/')
    #     continue
    if observation_time == 0:
        observation_time = -1
    if possible_flares == 0 :
        possible_flares = -1
    # Y = np.column_stack((flare_amplitude, flare_time_scale, flare_energy))
    # np.savetxt(output_dir + '/' + M_dwarves + '/All_Flares.csv', Y, delimiter=',')
    # f = open(output_dir + '/' + M_dwarves + '/Host_Statistics.txt', 'w')
    # f.write('Flares/Day: ' + str(flare_count*60*24/observation_time) + '\n' + 'Possible Flares: ' + str(possible_flares) + '\n' +
    #         'Accepted Flares: ' + str(flare_count) + '\n' + 'Accepted/Possible Flares: ' + str(flare_count/possible_flares) + '\n' +
    #         'Observation Time (min): ' + str(observation_time))
    f.close()
    item_count += 1
    observation_time = np.ones(len(Teff))*observation_time
    # ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, total_flare_energies, Teff, radius, TOI_period, observation_time, transit_epoch_list, periastron_epoch_list, a_list, e_list, min_approach_list))
    # with open('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_Flares.csv', "a") as f:
    #     np.savetxt(f, ZZ, delimiter=",", fmt='%s')
    #     f.close()
# g = open(output_dir + '/All_Exoplanet_Stats.txt', 'w')
# g.write('Total Flares: ' + str(total_flares) + '\n' + 'Net Observation Time (Days): ' + str(total_observation_time/(60*24)))
# g.close()