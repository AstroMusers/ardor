# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 21:47:16 2023

@author: Nate Whitsett
"""

import aflare as af
import Flare
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os
import time as t

def exp_decay(x, a, b, c):
    '''
    

    Parameters
    ----------
    x : numpy array
        Time axis data
    a : float
        multiplicative parameter of the exponential decay function
    b : float
        time constant parameter of the exponential decay function
    c : float
        Time-offset for the exponential decay function

    Returns
    -------
    numpy array
        Gives exponential decay model output with the given parameters

    '''
    return a * np.exp(-b * x) + c
def amp_log_normal():
    return np.random.lognormal(np.log(0.020744), sigma=np.log(4.33924339))


def FWHM_uniform():
    return np.random.uniform(0.001388888, 0.041)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

TESS_Folder_ID = os.listdir('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/TESS Data/G Type TESS Data')

TOI_ID_list = []
test_flare_list = []
lc_num = 0
item_count = 0
random_sample = []
check = 0
test_time = 0
files = 0

total_flares_detected_T1 = 0
total_flares_detected_T2 = 0
inject_flare_count = 0

inject_flare_recovered_T1 = 0
inject_flare_recovered_T2 = 0

false_positive_T1 = 0
false_positive_T2 = 0

manual_inspect_flare = 0


test_flare_list_T1 = []
test_flare_list_T2 = []
tier_1_time = 0
tier_2_time = 0
sigmas = [2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5]
# sigmas = [3]
for cutoff in sigmas:
    for G_Stars in range(20):
        G_Stars = TESS_Folder_ID[20]
        if files >= 20:
            files = 0
            break
        ##Iteration Scheme
        a = os.listdir('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/TESS Data/G Type TESS Data/' + G_Stars)
        
        ##Trackable values per star
        location_list = []
        for folders in range(1):
            folders = a[0]
            if files >= 20:
                break
            flare_inject_dict_T1 = dict()
            flare_inject_dict_T2 = dict()
            inject_location_index = []
            b, pdcsap_flux, pdcsap_error = Flare.TESS_data_extract('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/TESS Data/G Type TESS Data/' + G_Stars + '/' + folders, PDCSAP_ERR=True)
            time, flux, pdcsap_error = Flare.delete_nans(b, pdcsap_flux, pdcsap_error)
            if folders.endswith('a_fast-lc.fits') == True:
                test_time += len(time)*(0.33333333)
                fast = True
            elif folders.endswith('a_fast-lc.fits') == False:  
                test_time += len(time)*2
                fast = False
            
            for flares in range(10):
                location = np.random.randint(50, len(flux)-50)
                counter = 0 
                for locations in location_list:
                    while location > locations - 100 and location < locations + 100 and counter < 10000:
                        location = np.random.randint(50, len(flux)-50)
                        counter += 1
                inject_location_index.append(location)
                sample_baseline = flux[location-300:location+300]
                normalized_sample = sample_baseline/np.median(sample_baseline)
                FWHM = 0
                while FWHM 
                FWHM = FWHM_uniform()
                amp = amp_log_normal()
                flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
                normalized_sample_inject = normalized_sample + flare_inject
                flux[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
                baseline_error = np.std(normalized_sample)
                flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, False, False, cutoff]
                flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, False, False, cutoff]
                inject_flare_count += 1
                location_list.append(location)
                
            t0 = t.time()
            print(cutoff)
            # detrend_flux = Flare.SMA_detrend(time, flux, pdcsap_error, time_scale = 80)
            plt.plot(time, flux)
            plt.show()
            flares, lengths = Flare.flare_ID(flux, cutoff, fast=fast)
            index = 0
            for flare_events in flares:
                if flare_events >= 50 and len(flux) - flare_events > 50:
                    new_time = time[flare_events-50:flare_events+50]
                    new_data = flux[flare_events-50:flare_events+50]
                    new_error = pdcsap_error[flare_events-50:flare_events+50]
                elif flare_events < 50:
                    new_time = time[0+flare_events:flare_events+50]
                    new_data = flux[0+flare_events:flare_events+50]
                    new_error = pdcsap_error[0+flare_events:flare_events+50]
                elif len(flux) - flare_events < 50:
                    new_time = time[flare_events:]
                    new_data = flux[flare_events:]
                    new_error = pdcsap_error[flare_events:]
                recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
                norm_time = time[flare_events]
                events = np.where(new_data == recenter)[0][0]
                peak_centered_event = flare_events + (events-50)
                criteria1 = False
                flare_location = find_nearest(inject_location_index, flare_events)
                if np.abs(flare_location - flare_events) < 10 and flare_inject_dict_T1[flare_location][4] == False:
                    inject_flare_recovered_T1 += 1
                    flare_inject_dict_T1[flare_location][3] = True
                else:
                    false_positive_T1 += 1
                    flare_inject_dict_T1[flare_events] = [new_data[0]-np.median(new_data), 0.1, baseline_error, False, True, cutoff]
                total_flares_detected_T1 += 1
                t1 = t.time()
                tier_1_time += t1-t0
                t2 = t.time()
                new_time = (new_time - new_time[events])*24*60
                if len(new_data)  == 100 :
                    if lengths[index] >= 25:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000, sigma = error[events:events+30], absolute_sigma=True)
                        squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(error[events:events+30])
                        chi_squared = np.sum(squares)/27
                    elif lengths[index] >= 15 and lengths[index] < 25:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000, sigma = error[events:events+20],absolute_sigma=True)
                        squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(error[events:events+20])
                        chi_squared = np.sum(squares)/17
                    elif lengths[index] > 5 and lengths[index] < 15:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000, sigma = error[events:events+10],absolute_sigma=True)
                        squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(error[events:events+10])
                        chi_squared = np.sum(squares)/7
                    elif lengths[index] <= 5:
                        alles_data = new_data/np.median(new_data)
                        error = new_error/np.median(new_data)
                        popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
                        squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(error[events:events+7]**2)
                        chi_squared = np.sum(squares)/4

                    emp_FWHM = -np.log(0.5)/(popt[1])/(24*60)
                    t3 = t.time()
                    tier_2_time += t3 - t2
                    if chi_squared < 1.5 and popt[1] > 0 and popt[0] > 0:
                        flare_location = find_nearest(inject_location_index, flare_events)
                        if np.abs(flare_location - flare_events) < 10 and flare_inject_dict_T2[flare_location][4] == False:
                            inject_flare_recovered_T2 += 1
                            flare_inject_dict_T2[flare_location][3] = True
                        else:
                            # plt.scatter(new_time, alles_data)
                            # plt.show()
                            # a = input('Continue?')
                            # plt.clf()
                            # if a == 'y':
                            #     manual_inspect_flare += 1
                            #     flare_inject_dict_T1[peak_centered_event][3] = True
                            #     flare_inject_dict_T1[peak_centered_event][4] = False
                            #     flare_inject_dict_T2[peak_centered_event] = [popt[0], emp_FWHM, baseline_error, True, False]
                            # else:
                            false_positive_T2 += 1
                            flare_inject_dict_T2[peak_centered_event] = [popt[0], emp_FWHM, baseline_error, False, True, cutoff]
                        total_flares_detected_T2 += 1
                        index += 1
                
            test_flare_list_T1.append(flare_inject_dict_T1)
            test_flare_list_T2.append(flare_inject_dict_T2)
            lc_num += 1
            print('LC #: ' + str(lc_num), files)
            print(inject_flare_recovered_T2)
            files += 1
        

print('Result of test:')
print('Files Tested: ' + str(files))
print('Total Time Analyzed: ' + str(test_time) + 'minutes')

print('Tier 1 Detections: ' + str(total_flares_detected_T1))
print('Tier 2 Detections: ' + str(total_flares_detected_T2))

print('Tier 1 False Detections: '+ str(false_positive_T1))
print('Tier 2 False Detections: '+ str(false_positive_T2))

print('Tier 1 True Flares Recovered: ' + str(inject_flare_recovered_T1))
print('Tier 2 True Flares Recovered: ' + str(inject_flare_recovered_T2))

print('Tier 1 True Flares Missed: ' + str(inject_flare_count - inject_flare_recovered_T1))
print('Tier 2 True Flares Missed: ' + str(inject_flare_count - inject_flare_recovered_T2))

print('Tier 1 %-True Flares Recovered: ' + str(inject_flare_recovered_T1/inject_flare_count))
print('Tier 2 %-True Flares Recovered: ' + str(inject_flare_recovered_T2/inject_flare_count))

print('Tier 1 False Positive Rate: ' + str(false_positive_T1/test_time))
print('Tier 2 False Positive Rate: ' + str(false_positive_T2/test_time))

print('Tier 2 Manually Inspected, Non Injected Flares: ' + str(manual_inspect_flare))

amp_list_T1=[]
FWHM_list_T1 = []
result_T1 = []
error_T1 = []
false_positive_T1 = []
sigma_T1 = []

amp_list_T2=[]
FWHM_list_T2 = []
result_T2 = []
error_T2 = []
false_positive_T2 = []
sigma_T2 = []
for d in test_flare_list_T1:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list_T1.append(value[0])
        FWHM_list_T1.append(value[1])
        error_T1.append(value[2])
        result_T1.append(value[3])
        false_positive_T1.append(value[4])
        sigma_T1.append(value[5])
        

for d in test_flare_list_T2:
    for key,value in d.items():
        # notice the difference here, instead of appending a nested list
        # we just append the key and value
        # this will make temp_list something like: [a0, 0, a1, 1, etc...]
        amp_list_T2.append(value[0])
        FWHM_list_T2.append(value[1])
        error_T2.append(value[2])
        result_T2.append(value[3])
        false_positive_T2.append(value[4])
        sigma_T2.append(value[5])
        
ZZ = np.stack((amp_list_T1, FWHM_list_T1, error_T1, result_T1, false_positive_T1, sigma_T1))
ZZ = np.transpose(ZZ)
np.savetxt('C:/Users/natha/OneDrive - Washington University in St. Louis/Injection_Test_T1.csv', ZZ, delimiter=',')

Z = np.stack((amp_list_T2, FWHM_list_T2, error_T2, result_T2, false_positive_T2, sigma_T2))
Z = np.transpose(Z)
np.savetxt('C:/Users/natha/OneDrive - Washington University in St. Louis/Injection_Test_T2.csv', Z, delimiter=',')