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
import collections as c
from matplotlib import pyplot as plt
from scipy.stats import poisson
warnings.filterwarnings("ignore")

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 11,
        }
def amp_log_normal():
    return np.random.lognormal(np.log(0.025), sigma=np.log(4.33924339))

def FWHM_uniform():
    return np.random.uniform(0.0005, 0.041)


def Flare_Injection(light_curve, sp_type = 'M', flare_type='Flaring', fast=False, detrend = True, rates = True,PR = False, sigma = 3, chi_sq = 1):
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
                    continue
                normalized_sample = sample_baseline
                FWHM = FWHM_uniform()
                amp = amp_log_normal()
                flare_inject = af.aflare1(time[location-1000:location+1000], time[location], FWHM, amp)
                normalized_sample_inject = normalized_sample + flare_inject
                data[location-1000:location+1000] = normalized_sample_inject
                location_list.append(location)
                flares += 1
                integral = simpson(af.aflare1(time, time[location], FWHM, amp), x=time)
                inject_location_index.append(location)
                if PR == False:
                    flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, True, False, integral]
                    flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, True, False, integral]
                if PR == True:
                    flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, True, False, sigma]
                    flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, True, False, sigma]
    if rates== False:
        for locations in location_list:
            location = locations
            counter = 0 
            sample_baseline = data[location-1000:location+1000]
            baseline_error = np.std(sample_baseline)
            normalized_sample = sample_baseline
            FWHM = FWHM_uniform()
            amp = amp_log_normal()
            if np.isnan(baseline_error) == True:
                continue
            flare_inject = af.aflare1(time[location-1000:location+1000], time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            data[location-1000:location+1000] = np.median(sample_baseline)*normalized_sample_inject
            flares += 1
            integral = simpson(af.aflare1(time, time[location], FWHM, amp), x=time)
            inject_location_index.append(location)
            if PR == False:
                flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, True, False, integral]
                flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, True, False, integral]
            if PR == True:
                flare_inject_dict_T1[location] = [amp, FWHM, baseline_error, True, False, sigma, chi_sq]
                flare_inject_dict_T2[location] = [amp, FWHM, baseline_error, True, False, sigma, chi_sq]
    LightCurve = c.namedtuple('LightCurve', ['time', 'flux', 'error'])
    lc = LightCurve(time, data, error)
    return lc, flare_inject_dict_T1, flare_inject_dict_T2


def Injection_Recovery(input_dir, output_dir, star_type = 'G', rate = False, old = False):
    lc_num = 0
    item_count = 0
    tier0_tau = []
    tier1_tau = []
    tier2_tau = []
    for iterations in range(10):
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
                    lc, flare_inject_dict_T1, flare_inject_dict_T2 = Flare_Injection(input_dir + '/' + M_dwarves + '/' + folders, detrend = True, rates = rate, sp_type = star_type)
                    if folders.endswith('a_fast-lc.fits') == True:
                        continue
                        fast = True
                    elif folders.endswith('a_fast-lc.fits') == False:  
                        fast = False
                    t2 = timer.time()
                    tier0_tau.append(t2-t1)
                    flare_events_T1, lengths = Flare.flare_ID(lc.flux, 3, fast = fast, injection = True, old = True)
                    t3 = timer.time()
                    tier1_tau.append(t3-t2)                
                    if len(flare_events_T1) != 0:
                        for flares in flare_events_T1:
                            for keys in flare_inject_dict_T1.keys():
                                if keys - 50 < flares and keys + 50 > flares:
                                    flare_inject_dict_T1[keys][4] = True
                    flare_events_T2 = Flare.tier2(lc.time, lc.flux, lc.error, flare_events_T1, lengths, chi_square_cutoff = 10, host_name = 'My_Host', csv = False,Sim = False, injection= True)
                    if len(flare_events_T2) != 0:
                        for flares in flare_events_T2:
                            for keys in flare_inject_dict_T2.keys():
                                if keys - 50 < flares and keys + 50 > flares:
                                    flare_inject_dict_T2[keys][4] = True
                    for keys in flare_inject_dict_T1.keys():
                        if flare_inject_dict_T1[keys][4] == False:
                            del flare_inject_dict_T2[keys]
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
                            amp_list_T1.append(value[0])
                            FWHM_list_T1.append(value[1])
                            error_T1.append(value[2])
                            result_T1.append(value[3])
                            injected_T1.append(value[4])
                            integral_T1.append(value[5])
                    for d in test_flare_list_T2:
                        for key,value in d.items():
                            amp_list_T2.append(value[0])
                            FWHM_list_T2.append(value[1])
                            error_T2.append(value[2])
                            result_T2.append(value[3])
                            injected_T2.append(value[4])
                            integral_T2.append(value[5])
                    ZZ = np.stack((amp_list_T1, FWHM_list_T1, error_T1, integral_T1,injected_T1, result_T1))
                    ZZ = np.transpose(ZZ)
                    output_file_T1 = '/Injection_Recovery_T1.csv'
                    output_file_T2 = '/Injection_Recovery_T2.csv'
                    with open(output_dir + output_file_T1, "a") as f:
                        np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                        f.close()
                    Z = np.stack((amp_list_T2, FWHM_list_T2, error_T2, integral_T2,injected_T2, result_T2))
                    Z = np.transpose(Z)
                    with open(output_dir + output_file_T2, "a") as f:
                        np.savetxt(f, Z, delimiter=",", fmt='%s')
                        f.close()
            except:
                continue
    

def Injection_Recovery_Grid(data_dir, grid_dir, label = 'T1', energy = False):
    
    data = pd.read_csv(data_dir)
    data.columns.values[0] = 'Amplitude'
    data.columns.values[1] = 'FWHM'
    data.columns.values[2] = 'Error'
    data.columns.values[3] = 'Integral'
    data.columns.values[5] = 'Injected?'
    data.columns.values[4] = 'Accepted?'
    amp = list(data['Amplitude'])
    FWHM = list(data['FWHM'])
    injection = list(data['Injected?'])
    Bool = list(data['Accepted?'])
    integral = list(data['Integral'])
    amp_bins = np.logspace(-3, 0, num=17)
    FWHM_bins = np.linspace(0.001388888, 0.041, num=17)
    int_bins = np.logspace(-5, -1, num=17)
    x = []
    y = []
    if energy == False:
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
        with open(grid_dir + '/Injection_Recovery_Grid_' + label + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(y)
    elif energy == True:
        for i in range(len(int_bins)-1):
            tmp = []
            for index in range(len(integral)):
                if integral[index] > int_bins[i] and integral[index] < int_bins[i+1]:
                    tmp.append([integral[index], FWHM[index], Bool[index]])
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

        with open(grid_dir + '/Injection_Recovery_Grid_Energy' + label + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(y)
        
def Precision_Recall(input_dir, output_dir, sigma_list, chi_set, star_type = 'G', rate = False):
    lc_num = 0
    item_count = 0
    tier0_tau = []
    tier1_tau = []
    tier2_tau = []
    for sigma in sigma_list:
        for M_dwarves in os.listdir(input_dir)[0:20]:
            ##Iteration Scheme
            a = os.listdir(input_dir + '/' + M_dwarves)
            print(item_count, M_dwarves) 
            
            ##Relevant parameters from TOI catalog
            file = 0
            ##Trackable values per star
            try:
                for folders in a:
                    t1 = timer.time()
                    lc, flare_inject_dict_T1, flare_inject_dict_T2 = Flare_Injection(input_dir + '/' + M_dwarves + '/' + folders, detrend = True, rates = rate, sp_type = star_type, PR = True, sigma = sigma, chi_sq = 10)
                    if folders.endswith('a_fast-lc.fits') == True:
                        continue
                        fast = True
                    elif folders.endswith('a_fast-lc.fits') == False:  
                        fast = False
                    t2 = timer.time()
                    tier0_tau.append(t2-t1)
                    flare_events_T1, lengths = Flare.flare_ID(lc.flux, sigma, fast = fast, injection = True)
                    t3 = timer.time()
                    tier1_tau.append(t3-t2)                
                    if len(flare_events_T1) != 0:
                        for flares in flare_events_T1:
                            for keys in flare_inject_dict_T1.keys():
                                if keys - 100 < flares and keys + 100 > flares:
                                    flare_inject_dict_T1[keys][4] = True
                                    flare_inject_dict_T2[keys][6] = 10
                    count = 0
                    for keys in flare_inject_dict_T1.keys():
                        count += int(flare_inject_dict_T1[keys][3])
                    for diff in range(abs(len(flare_events_T1)- count)):
                        flare_inject_dict_T1[diff] = [0, 0, 0, False, True, sigma, 10]
                    flare_events_T2 = Flare.tier2(lc.time, lc.flux, lc.error, flare_events_T1, lengths, chi_square_cutoff = 10, host_name = 'My_Host', csv = False,Sim = False, injection= True)
                    if len(flare_events_T2) != 0:
                        for flares in flare_events_T2:
                            for keys in flare_inject_dict_T2.keys():
                                if keys - 100 < flares and keys + 100 > flares:
                                    flare_inject_dict_T2[keys][4] = True
                                    flare_inject_dict_T2[keys][6] = 10
                    count = 0
                    for keys in flare_inject_dict_T2.keys():
                        count += int(flare_inject_dict_T2[keys][3])
                    for diff in range(abs(len(flare_events_T2)- count)):
                        flare_inject_dict_T2[diff] = [0, 0, 0, False, True, sigma, 10]
                    t4 = timer.time()
                    tier2_tau.append(t4-t3)
                    for keys in flare_inject_dict_T1.keys():
                        if flare_inject_dict_T1[keys][4] == False:
                            del flare_inject_dict_T2[keys]
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
                    lc_num += 1
                    print('LC #: ' + str(lc_num))
                    file += 1
                    amp_list_T1=[]
                    FWHM_list_T1 = []
                    result_T1 = []
                    error_T1 = []
                    injected_T1 = []
                    sigma_T1 = []
                    chi_sq_T1 = []
                    
                    amp_list_T2=[]
                    FWHM_list_T2 = []
                    result_T2 = []
                    error_T2 = []
                    injected_T2 = []
                    sigma_T2 = []
                    chi_sq_T2 = []
                    for d in test_flare_list_T1:
                        for key,value in d.items():
                            amp_list_T1.append(value[0])
                            FWHM_list_T1.append(value[1])
                            error_T1.append(value[2])
                            injected_T1.append(value[3])
                            result_T1.append(value[4])
                            sigma_T1.append(value[5])
                            chi_sq_T1.append(value[6])
                    for d in test_flare_list_T2:
                        for key,value in d.items():
                            amp_list_T2.append(value[0])
                            FWHM_list_T2.append(value[1])
                            error_T2.append(value[2])
                            injected_T2.append(value[3])
                            result_T2.append(value[4])
                            sigma_T2.append(value[5])
                            chi_sq_T2.append(value[6])
                    ZZ = np.stack((amp_list_T1, FWHM_list_T1, error_T1,injected_T1,result_T1,sigma_T1, chi_sq_T1))
                    ZZ = np.transpose(ZZ)
                    output_file_T1 = '/Precision_Recall_T1.csv'
                    output_file_T2 = '/Precision_Recall_T2.csv'
                    with open(output_dir + output_file_T1, "a") as f:
                        np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                        f.close()
                    Z = np.stack((amp_list_T2, FWHM_list_T2, error_T2,injected_T2,result_T2,sigma_T2,chi_sq_T2))
                    Z = np.transpose(Z)
                    with open(output_dir + output_file_T2, "a") as f:
                        np.savetxt(f, Z, delimiter=",", fmt='%s')
                        f.close()
            except:
                continue

def Precision_Recall_Data_Processing(PR_T1_dir, PR_T2_dir, sigma_set, chi_set,output_dir, tier = 'T1'):
    data_T1 = pd.read_csv(PR_T1_dir, index_col=None, header=None)
    data_T2 = pd.read_csv(PR_T2_dir, index_col=None, header=None)
    Bool = data_T2[3].array
    False_Pos = data_T2[4].array
    sigma = data_T2[5].array
    chi_sq = data_T2[6].array
    Bool1 = data_T1[3].array
    False_Pos1 = data_T1[4].array
    sigma1 = data_T1[5].array
    chi_sq1 = data_T1[6].array
    FN = 0
    TP = 0
    FP = 0
    x = []
    y = []
    x1 = []
    y1 = []
    count = 0
    x_err_l1 = []
    x_err_h1 = []
    y_err_l1 = []
    y_err_h1 = []
    
    if tier == 'T1':
        counts = 0
        
        for bins in range(len(sigma_set)):
            for index, num in enumerate(sigma1):
                if sigma1[index] == sigma_set[counts]:
                    if Bool1[index] == 1 and False_Pos1[index] == 1:
                        TP += 1
                        count += 1
                    if Bool1[index] == 1 and False_Pos1[index] == 0:
                        FN += 1
                        count += 1
                    if Bool1[index] == 0 and False_Pos1[index] == 1:
                        FP += 1
                        count +=1
            
            precision_list = []
            recall_list = []
            precision = TP/(TP+FP)
            recall = TP/(TP+FN)
            for trials in range(10000):
                TP_t = poisson.rvs(TP)
                FP_t = poisson.rvs(FP)
                FN_t = poisson.rvs(FN)
                recall_list.append(((TP_t)/(TP_t+FN_t)))
                precision_list.append((TP_t)/(TP_t+FP_t))
            precision_list = np.sort(np.array(precision_list))
            recall_list = np.sort(np.array(recall_list))
            y_err_l1.append(np.abs(precision_list[9500]-precision))
            y_err_h1.append(np.abs(precision_list[500]-precision))
            x_err_l1.append(np.abs(recall_list[9500]-recall))
            x_err_h1.append(np.abs(recall_list[500]-recall))
            x1.append((TP)/(TP+FN))
            y1.append((TP)/(TP+FP))
            counts += 1
            TP = 0
            FP = 0
            FN = 0
        
        
        # upper = 100
        # lower = 94
        x_err_l = []
        x_err_h = []
        y_err_l = []
        y_err_h = []
        counts = 0
        for bins in range(len(sigma_set)):
            for index, num in enumerate(sigma):
                if sigma[index] == sigma_set[counts]:
                    if Bool[index] == 1 and False_Pos[index] == 1:
                        TP += 1
                        count += 1
                    if Bool[index] == 1 and False_Pos[index] == 0:
                        FN += 1
                        count += 1
                    if Bool[index] == 0 and False_Pos[index] == 1:
                        FP += 1
                        count +=1
            precision_list = []
            recall_list = []
            precision = TP/(TP+FP)
            recall = TP/(TP+FN)
            for trials in range(10000):
                TP_t = poisson.rvs(TP)
                FP_t = poisson.rvs(FP)
                FN_t = poisson.rvs(FN)
                recall_list.append((TP_t)/(TP_t+FN_t))
                precision_list.append((TP_t)/(TP_t+FP_t))
               
            precision_list = np.sort(np.array(precision_list))
            recall_list = np.sort(np.array(recall_list))
            
            y_err_l.append(np.abs(precision_list[9500]-precision))
            y_err_h.append(np.abs(precision_list[500]-precision))
            x_err_l.append(np.abs(recall_list[9500]-recall))
            x_err_h.append(np.abs(recall_list[500]-recall))
            x.append((TP)/(TP+FN))
            y.append((TP)/(TP+FP))
            counts +=1 
            TP = 0
            FP = 0
            FN = 0
        x_err_h1[0] = 0
        x_err_h[0] = 0
        x_err1 = [x_err_l1, x_err_h1]
        y_err1 = [y_err_l1, y_err_h1]
        x_err = [x_err_l, x_err_h]
        y_err = [y_err_l, y_err_h]
        print(simpson(y, x=x))
        # plt.xlim(0, 1.1)
        # x = x/np.max(x)
        # x1 = x1/np.max(x1)
        # plt.savefig('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Injection Tests/Precision_Recall/Precision_Recall_Curve.png', dpi=600, bbox_inches="tight")
        plt.show()
    elif tier == 'T2':
        counts = 0
        
        for bins in range(len(chi_set)):
            for index, num in enumerate(chi_sq1):
                if chi_sq1[index] == chi_set[counts]:
                    if Bool1[index] == 1 and False_Pos1[index] == 1:
                        TP += 1
                        count += 1
                    if Bool1[index] == 1 and False_Pos1[index] == 0:
                        FN += 1
                        count += 1
                    if Bool1[index] == 0 and False_Pos1[index] == 1:
                        FP += 1
                        count +=1
            precision_list = []
            recall_list = []
            precision = TP/(TP+FP)
            recall = TP/(TP+FN)
            for trials in range(10000):
                TP_t = poisson.rvs(TP)
                FP_t = poisson.rvs(FP)
                FN_t = poisson.rvs(FN)
                recall_list.append(((TP_t)/(TP_t+FN_t)))
                precision_list.append((TP_t)/(TP_t+FP_t))
            precision_list = np.sort(np.array(precision_list))
            recall_list = np.sort(np.array(recall_list))
            y_err_l1.append(np.abs(precision_list[9500]-precision))
            y_err_h1.append(np.abs(precision_list[500]-precision))
            x_err_l1.append(np.abs(recall_list[9500]-recall))
            x_err_h1.append(np.abs(recall_list[500]-recall))
            x1.append((TP)/(TP+FN))
            y1.append((TP)/(TP+FP))
            counts += 1
            TP = 0
            FP = 0
            FN = 0
        
        
        # upper = 100
        # lower = 94
        x_err_l = []
        x_err_h = []
        y_err_l = []
        y_err_h = []
        counts = 0
        for bins in range(len(chi_set)):
            for index, num in enumerate(chi_sq):
                if chi_sq[index] == chi_set[counts]:
                    if Bool[index] == 1 and False_Pos[index] == 1:
                        TP += 1
                        count += 1
                    if Bool[index] == 1 and False_Pos[index] == 0:
                        FN += 1
                        count += 1
                    if Bool[index] == 0 and False_Pos[index] == 1:
                        FP += 1
                        count +=1
            precision_list = []
            recall_list = []
            precision = TP/(TP+FP)
            recall = TP/(TP+FN)
            for trials in range(10000):
                TP_t = poisson.rvs(TP)
                FP_t = poisson.rvs(FP)
                FN_t = poisson.rvs(FN)
                recall_list.append((TP_t)/(TP_t+FN_t))
                precision_list.append((TP_t)/(TP_t+FP_t))
               
            precision_list = np.sort(np.array(precision_list))
            recall_list = np.sort(np.array(recall_list))
            
            y_err_l.append(np.abs(precision_list[9500]-precision))
            y_err_h.append(np.abs(precision_list[500]-precision))
            x_err_l.append(np.abs(recall_list[9500]-recall))
            x_err_h.append(np.abs(recall_list[500]-recall))
            x.append((TP)/(TP+FN))
            y.append((TP)/(TP+FP))
            counts +=1 
            TP = 0
            FP = 0
            FN = 0
        x_err_h1[0] = 0
        x_err_h[0] = 0
        x_err1 = [x_err_l1, x_err_h1]
        y_err1 = [y_err_l1, y_err_h1]
        x_err = [x_err_l, x_err_h]
        y_err = [y_err_l, y_err_h]
        df = pd.DataFrame()
        df['x'] = x
        df['x-'] = x_err[0]
        df['x+'] = x_err[1]
        df['y'] = y
        df['y-'] = y_err[0]
        df['y+'] = y_err[1]
        df['x1'] = x1
        df['x-1'] = x_err1[0]
        df['x+1'] = x_err1[1]
        df['y1'] = y1
        df['y-1'] = y_err1[0]
        df['y+1'] = y_err1[1]
        df.to_csv(output_dir, index = None)

# Injection_Recovery('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_LC', 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_IR')
Injection_Recovery_Grid('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_IR/Injection_Recovery_T1.csv', 
                        'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_IR',
                        label = 'T1', energy = False)
# Precision_Recall('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_LC', 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_IR', [2, 2.5, 3, 3.5, 4, 4.5, 5], [1, 2, 5, 7.5, 10, 15, 20])
# Precision_Recall('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/M_Type_LC', 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/M_Type_IR', [2, 2.5, 3, 3.5, 4, 4.5, 5], [1, 2, 5, 7.5, 10, 15, 20])

# Precision_Recall_Data_Processing('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/M_Type_IR/Precision_Recall_T1.csv', 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/G_Type_IR/Precision_Recall_T1.csv', [2, 2.5, 3, 3.5, 4, 4.5, 5], chi_set=[1, 2, 5, 7.5, 10, 15, 20], tier = 'T1')