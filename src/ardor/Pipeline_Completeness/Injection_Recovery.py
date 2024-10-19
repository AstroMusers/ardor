# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 09:35:59 2023

@author: Nate Whitsett
"""

import aflare as af
import Flare
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os
import planck_law
from scipy.integrate import simpson
import time as timer
def amp_log_normal():
    return np.random.lognormal(np.log(0.020744), sigma=np.log(4.33924339))

def FWHM_uniform():
    return np.random.uniform(0.001388888, 0.041)
def Flare_injection(light_curve, sp_type = 'M', flare_type='Flaring', fast=False, detrend = True, rates = True):
    location_list = []
    ## Approximate flare rate per 2 minute cadence of flaring M/F stars (~0.5 flares/day)
    if flare_type == 'Flaring':
        if sp_type == 'M':
            rate = 2.8e-4
        if sp_type == 'F':
            rate = 6e-05
        if sp_type == 'G':
            rate = 1.18e-4
        if sp_type == 'K':
            rate = 1.19e-4
    ## Poor statistics on this, but G type stars flare ~2e-5 per 2 minute cadence
    elif flare_type == 'Not Flaring':
        rate = 2.78e-8
    ## Adjust times for 20s cadence
    if fast == True:
        rate /= 6
    time, data, error = Flare.TESS_data_extract(light_curve, PDCSAP_ERR=True)
    if detrend == True:
        data = Flare.SMA_detrend(time, data, error, 300)
    ## Iterate over the time scale of the light curve
    flares = 0
    location_list = []
    if rates == False:
        counter = 0
        for flares in range(15):
            location = np.random.randint(200, len(data)-200)
            for locations in location_list:
                while location > locations - 100 and location < locations + 100 and counter < 10000:
                    location = np.random.randint(50, len(data)-50)
                    counter += 1
            location_list.append(location)
    inject_location_index = []
    flare_inject_dict_T1 = dict()
    flare_inject_dict_T2 = dict()
    if rates == True:
        for interval in range(len(time) - 200):
            flare_check = np.random.random()
            flare_rate = rate
            if flare_rate >= flare_check:
                location = interval
                counter = 0 
                for locations in location_list:
                    while location > locations - 100 and location < locations + 100 and counter < 10000:
                        location = np.random.randint(50, len(data)-50)
                        counter += 1
                sample_baseline = data[location-300:location+300]
                baseline_error = np.std(sample_baseline)
                normalized_sample = sample_baseline
                FWHM = FWHM_uniform()
                amp = amp_log_normal()
                flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
                normalized_sample_inject = normalized_sample + flare_inject
                data[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
                location_list.append(location)
                flares += 1
                integral = simpson(af.aflare1(time, time[location], FWHM, amp), x=time)
                inject_location_index.append(location)
                flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, False, True, integral]
                flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, False, True, integral]
    if rates== False:
        for locations in location_list:
            location = locations
            counter = 0 
            sample_baseline = data[location-300:location+300]
            baseline_error = np.std(sample_baseline)
            normalized_sample = sample_baseline
            FWHM = FWHM_uniform()
            amp = amp_log_normal()
            flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            data[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
            flares += 1
            integral = simpson(af.aflare1(time, time[location], FWHM, amp), x=time)
            inject_location_index.append(location)
            flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, False, True, integral]
            flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, False, True, integral]
    return data, time, error, flare_inject_dict_T1, flare_inject_dict_T2



TESS_Folder_ID = os.listdir('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G Type TESS Data')


TOI_ID_list = []
test_flare_list_T1 = []
test_flare_list_T2 = []
lc_num = 0
item_count = 0
random_sample = []
check = 0


tier0_tau = []
tier1_tau = []
tier2_tau = []

for stars in range(5):
    sample = np.random.randint(0, len(TESS_Folder_ID))
    if TESS_Folder_ID[sample] in random_sample:
        sample = np.random.randint(0, len(TESS_Folder_ID))
    random_sample.append(TESS_Folder_ID[sample])
    
for M_dwarves in random_sample:
    ##Iteration Scheme
    a = os.listdir('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G Type TESS Data/' + M_dwarves)
    print(item_count, M_dwarves) 
    
    ##Relevant parameters from TOI catalog
    flare_count = 1
    file = 0
    ##Trackable values per star
    for folders in a:
        t1 = timer.time()
        if file == 0:
            location_list = []
            inject_location_index = []
            flux, time, error, flare_inject_dict_T1, flare_inject_dict_T2 = Flare_injection('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G Type TESS Data/' + M_dwarves + '/' + folders, detrend = True, rates = False)
            if folders.endswith('a_fast-lc.fits') == True:
                fast = True
            elif folders.endswith('a_fast-lc.fits') == False:  
                fast = False
            t2 = timer.time()
            tier0_tau.append(t2-t1)
            flare_events_T1, lengths = Flare.flare_ID(flux, 3, fast = fast)
            t3 = timer.time()
            tier1_tau.append(t3-t2)
            if len(flare_events_T1) != 0:
                for flares in flare_events_T1:
                    for keys in flare_inject_dict_T1.keys():
                        if keys - 50 < flares and keys + 50 > flares:
                            flare_inject_dict_T1[keys][3] = True        
            flare_events_T2 = Flare.tier2(time, flux, error, flare_events_T1, lengths, chi_square_cutoff = 5, host_name = 'My_Host', csv = False,Sim = False, injection= True)
            if len(flare_events_T2) != 0:
                for flares in flare_events_T2:
                    for keys in flare_inject_dict_T2.keys():
                        if keys - 50 < flares and keys + 50 > flares:
                            flare_inject_dict_T2[keys][3] = True
            t4 = timer.time()
            tier2_tau.append(t4-t3)
            XX = np.column_stack((np.array(np.array(tier0_tau).mean()), np.array(np.array(tier1_tau).mean()),  np.array(np.array(tier2_tau).mean())))
            with open("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/Time_Stats/Time_Stats.csv", "a") as f:
                np.savetxt(f, XX, delimiter=",", fmt='%s')
                f.close()
            tier0_tau = []
            tier1_tau = []
            tier2_tau = []
            test_flare_list_T1.append(flare_inject_dict_T1)
            test_flare_list_T2.append(flare_inject_dict_T2)
            p = 0
            p1 = 0
            for keys in flare_inject_dict_T1.keys():
                if flare_inject_dict_T1[keys][3] == 1:
                    p += 1
            for keys in flare_inject_dict_T2.keys():
                if flare_inject_dict_T1[keys][3] == 1:
                    p1 += 1
            lc_num += 1
            print('LC #: ' + str(lc_num))
            file += 1
        else:
            continue
        
amp_list_T1=[]
FWHM_list_T1 = []
result_T1 = []
error_T1 = []
integral_T1 = []
injected_T1 = []

amp_list_T2=[]
FWHM_list_T2 = []
result_T2 = []
error_T2 = []
integral_T2 = []
injected_T2 = []
for d in test_flare_list_T1:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list_T1.append(value[0])
        FWHM_list_T1.append(value[1])
        error_T1.append(value[2])
        result_T1.append(value[3])
        injected_T1.append(value[4])
        integral_T1.append(value[5])
for d in test_flare_list_T2:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list_T2.append(value[0])
        FWHM_list_T2.append(value[1])
        error_T2.append(value[2])
        result_T2.append(value[3])
        injected_T2.append(value[4])
        integral_T2.append(value[5])
        
ZZ = np.stack((amp_list_T1, FWHM_list_T1, error_T1, integral_T1, result_T1, injected_T1))
ZZ = np.transpose(ZZ)
np.savetxt('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/Params/Injection_Recovery_G_Type_T1.csv', ZZ, delimiter=',')

Z = np.stack((amp_list_T2, FWHM_list_T2, error_T2, integral_T2, result_T2, injected_T2))
Z = np.transpose(Z)
np.savetxt('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/Params/Injection_Recovery_G_Type_T2.csv', Z, delimiter=',')