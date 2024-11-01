# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 09:35:59 2023

@author: Nate Whitsett
"""


import ardor.Flares.aflare as af
import ardor.Flares.Flare as Flare
import numpy as np
import os
from scipy.integrate import simpson
import time as timer
import warnings
import pandas as pd
import csv
from matplotlib import pyplot as plt
warnings.filterwarnings("ignore")
def amp_log_normal():
    return np.random.lognormal(np.log(0.025), sigma=np.log(4.33924339))

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
            rate = 5e-5
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
        lc = Flare.tier0(light_curve, scale=401)
        time = lc.time
        data = lc.flux
        error = lc.error
    ## Iterate over the time scale of the light curve
    flares = 0
    location_list = [500*i + x for i, x in enumerate(sorted(np.random.choice(range(200, len(data)-7500), 15)))]
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
                    while location > locations - 800 and location < locations + 800 and counter < 10000:
                        location = np.random.randint(50, len(data)-50)
                        counter += 1
                sample_baseline = data[location-300:location+300]
                baseline_error = np.std(sample_baseline)
                if np.isnan(baseline_error) == True:
                    print('Nan!')
                    continue
                normalized_sample = sample_baseline
                FWHM = FWHM_uniform()
                amp = amp_log_normal()
                flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
                normalized_sample_inject = normalized_sample + flare_inject
                data[location-300:location+300] = normalized_sample_inject
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
            if np.isnan(baseline_error) == True:
                print('Nan!')
                continue
            flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            data[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
            flares += 1
            integral = simpson(af.aflare1(time, time[location], FWHM, amp), x=time)
            inject_location_index.append(location)
            flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, False, True, integral]
            flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, False, True, integral]
    return data, time, error, flare_inject_dict_T1, flare_inject_dict_T2


def Injection_Recovery(input_dir, output_dir, star_type = 'G', rate = False):
    lc_num = 0
    item_count = 0
    tier0_tau = []
    tier1_tau = []
    tier2_tau = []
    for trials in range(10):
        for M_dwarves in os.listdir(input_dir):
            ##Iteration Scheme
            a = os.listdir(input_dir + '/' + M_dwarves)
            print(item_count, M_dwarves) 
            
            ##Relevant parameters from TOI catalog
            file = 0
            ##Trackable values per star
            try:
                for folders in a:
                    t1 = timer.time()
                    flux, time, error, flare_inject_dict_T1, flare_inject_dict_T2 = Flare_injection(input_dir + '/' + M_dwarves + '/' + folders, detrend = True, rates = rate, sp_type = star_type)
                    if folders.endswith('a_fast-lc.fits') == True:
                        continue
                        fast = True
                    elif folders.endswith('a_fast-lc.fits') == False:  
                        fast = False
                    t2 = timer.time()
                    tier0_tau.append(t2-t1)
                    flare_events_T1, lengths = Flare.flare_ID(flux, 2.5, fast = fast)
                    t3 = timer.time()
                    tier1_tau.append(t3-t2)
                    if len(flare_events_T1) != 0:
                        for flares in flare_events_T1:
                            for keys in flare_inject_dict_T1.keys():
                                if keys - 50 < flares and keys + 50 > flares:
                                    flare_inject_dict_T1[keys][3] = True        
                    flare_events_T2 = Flare.tier2(time, flux, error, flare_events_T1, lengths, chi_square_cutoff = 20, host_name = 'My_Host', csv = False,Sim = False, injection= True)
                    if len(flare_events_T2) != 0:
                        for flares in flare_events_T2:
                            for keys in flare_inject_dict_T2.keys():
                                if keys - 100 < flares and keys + 100 > flares:
                                    flare_inject_dict_T2[keys][3] = True
                    t4 = timer.time()
                    tier2_tau.append(t4-t3)
                    XX = np.column_stack((np.array(np.array(tier0_tau).mean()), np.array(np.array(tier1_tau).mean()),  np.array(np.array(tier2_tau).mean())))
                    with open(output_dir + '/Time_Stat.csv', "a") as f:
                        np.savetxt(f, XX, delimiter=",", fmt='%s')
                        f.close()
                    tier0_tau = []
                    tier1_tau = []
                    tier2_tau = []
                    test_flare_list_T1 = []
                    test_flare_list_T2 = []
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
                    with open(output_dir + '/Injection_Recovery_T1.csv', "a") as f:
                        np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                        f.close()
                    Z = np.stack((amp_list_T2, FWHM_list_T2, error_T2, integral_T2, result_T2, injected_T2))
                    Z = np.transpose(Z)
                    with open(output_dir + '/Injection_Recovery_T2.csv', "a") as f:
                        np.savetxt(f, Z, delimiter=",", fmt='%s')
                        f.close()
            except:
                continue
    

def Injection_Recovery_Grid(data_dir, grid_dir, label = 'T1'):
    
    data = pd.read_csv(data_dir)
    data.columns.values[0] = 'Amplitude'
    data.columns.values[1] = 'FWHM'
    data.columns.values[2] = 'Error'
    data.columns.values[3] = 'Integral'
    data.columns.values[4] = 'Accepted?'
    amp = list(data['Amplitude'])
    FWHM = list(data['FWHM'])
    error = list(data['Error'])
    Bool = list(data['Accepted?'])
    integral = list(data['Integral'])
    amp_bins = np.logspace(-3, 0, num=17)
    FWHM_bins = np.linspace(0.001388888, 0.041, num=17)
    error_bins = np.linspace(0.0001, 0.25, num=16)
    x = []
    y = []

    for i in range(len(amp_bins)-1):
        tmp = []
        for index in range(len(amp)):
            if amp[index] > amp_bins[i] and amp[index] < amp_bins[i+1]:
                tmp.append([amp[index], FWHM[index], Bool[index]])
        x.append(tmp)
    total = 0
    for i in range(len(FWHM_bins)):
        tmp = []
        for cells in x:
            count = 0
            pos = 0
            for flares in cells:
                if flares[1] > FWHM_bins[i] and flares[1] < FWHM_bins[i+1]:
                    count += 1
                    total += 1
                    if flares[2] == 1:
                        pos += 1
            if i == 16:
                continue
            if count == 0:
                tmp.append(np.nan)
            elif count != 0:
                tmp.append(pos/count * 100)
        y.append(tmp)

    with open(grid_dir + '/Injection_Recover_Grid_' + label + '.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(y)
Injection_Recovery('C:/Users/Nate Whitsett/Desktop/M_Type_LC/', 'C:/Users/Nate Whitsett/Desktop/M_Type_LC/Output/', rate = False)
# Injection_Recovery_Grid('C:/Users/Nate Whitsett/Desktop/M_Type_LC/Output/Injection_Recovery_T1.csv', 'C:/Users/Nate Whitsett/Desktop/M_Type_LC/Output/', label = 'T1')