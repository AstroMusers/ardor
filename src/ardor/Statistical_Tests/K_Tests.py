# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:42:31 2024

@author: whitsett.n
"""

import numpy as np
import pandas as pd
import os
from scipy.stats import uniform
from scipy.stats import ks_1samp
from skgof import ad_test
import csv
from itertools import zip_longest
from astropy.stats import kuiper
import ardor.Statistical_Tests.MLE as MLE
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
def K_Test_Sim_Kappa(DataFrame, output_dir, descriptor='G'):
    # interest_hosts = [['p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N']]
    data_frame = pd.DataFrame()
    set_hosts = set(DataFrame['Host_ID'])
    kappa_0KU = []
    kappa_25KU = []
    kappa_5KU = []
    kappa_1KU = []
    kappa_2KU = []
    kappa_4KU = []
    kappa_8KU = []
    kappa_0AD = []
    kappa_25AD = []
    kappa_5AD = []
    kappa_1AD = []
    kappa_2AD = []
    kappa_4AD = []
    kappa_8AD = []
    kappa_0KS = []
    kappa_25KS = []
    kappa_5KS = []
    kappa_1KS = []
    kappa_2KS = []
    kappa_4KS = []
    kappa_8KS = []
    KU_List = [kappa_0KU, kappa_25KU, kappa_5KU, kappa_1KU, kappa_2KU, kappa_4KU,kappa_8KU]
    AD_List = [kappa_0AD, kappa_25AD, kappa_5AD, kappa_1AD, kappa_2AD, kappa_4AD,kappa_8AD]
    KS_List = [kappa_0KS, kappa_25KS, kappa_5KS, kappa_1KS, kappa_2KS, kappa_4KS,kappa_8KS]
    for index, hosts in enumerate(set_hosts):

        if float(hosts[-1]) == 8 or float(hosts[-1]) == 4 or float(hosts[-1]) == 0 or float(hosts[-1]) == 2 or float(hosts[-1]) == 1:
            kappa = float(hosts[-1])
        elif float(hosts[-3:]) == 0.5:
            kappa = 0.5
        elif float(hosts[-4:]) == 0.25:
            kappa = 0.25
        peri_phases = np.array(DataFrame.loc[DataFrame['Host_ID'] == hosts, 'Phase'])
        peri_phases = np.sort(peri_phases)
        D, p_KS = ks_1samp(peri_phases, uniform.cdf, args=(0, 1))
        A, p_AD = ad_test(peri_phases, uniform(0,1), assume_sorted=True)
        a, p_KU = kuiper(peri_phases, uniform.cdf, args=(0, 1))
        if kappa == 0:
            kappa_0KU.append(p_KU)
            kappa_0KS.append(p_KS)
            kappa_0AD.append(p_AD)
        elif kappa == 0.25:
            kappa_25KU.append(p_KU)
            kappa_25KS.append(p_KS)
            kappa_25AD.append(p_AD)
        elif kappa == 0.5:
            kappa_5KU.append(p_KU)
            kappa_5KS.append(p_KS)
            kappa_5AD.append(p_AD)
        elif kappa == 1:
            kappa_1KU.append(p_KU)
            kappa_1KS.append(p_KS)
            kappa_1AD.append(p_AD)
        elif kappa == 2:
            kappa_2KU.append(p_KU)
            kappa_2KS.append(p_KS)
            kappa_2AD.append(p_AD)
        elif kappa == 4:
            kappa_4KU.append(p_KU)
            kappa_4KS.append(p_KS)
            kappa_4AD.append(p_AD)
        elif kappa == 8:
            kappa_8KU.append(p_KU)
            kappa_8KS.append(p_KS)
            kappa_8AD.append(p_AD)
    for index, kappas in enumerate([0,0.25,0.5,1,2,4,8]):
        sig_count_twoKU = 0
        sig_count_threeKU = 0
        for x in KU_List[index]:
            if x < 0.05:
                sig_count_twoKU += 1
            if x < 0.0027:
                sig_count_threeKU += 1
        KU_List[index].append(sig_count_twoKU/50)
        KU_List[index].append(sig_count_threeKU/50)
        sig_count_twoAD = 0
        sig_count_threeAD = 0
        for x in AD_List[index]:
            if x < 0.05:
                sig_count_twoAD += 1
            if x < 0.0027:
                sig_count_threeAD += 1
        AD_List[index].append(sig_count_twoAD/50)
        AD_List[index].append(sig_count_threeAD/50)
        sig_count_twoKS = 0
        sig_count_threeKS = 0
        for x in KS_List[index]:
            if x < 0.05:
                sig_count_twoKS += 1
            if x < 0.0027:
                sig_count_threeKS += 1
        KS_List[index].append(sig_count_twoKS/50)
        KS_List[index].append(sig_count_threeKS/50)
        data_frame['p_KU_' + str(kappas) + '_G'] = KU_List[index]
        data_frame['p_AD_' + str(kappas) + '_G'] = AD_List[index]
        data_frame['p_KS_' + str(kappas) + '_G'] = KS_List[index]
    data_frame.to_csv(output_dir, index = False)
    
def Peri_K_Test_Sampling(flares, N, error = True, peri_phases_upper = 0, peri_phases_lower = 0):
    AD_Sample_list = []
    KS_Sample_list = []
    if error == True:
        for samples in range(N):
            sample_flares = flares
            check = np.random.random()
            if check > 0.5:
                sample_flares = sample_flares + np.abs(np.random.normal(0, scale=np.abs(peri_phases_upper)))
            elif check < 0.5:
                sample_flares = sample_flares - np.abs(np.random.normal(0, scale=np.abs(peri_phases_lower)))
            for index in range(len(sample_flares)):
                if sample_flares[index] > 1:
                    sample_flares[index] = sample_flares[index] - 1 
                if flares[index] < 0:
                    sample_flares[index] = sample_flares[index] + 1
            sample_flares = np.sort(sample_flares)
            D, p_ks_samp = ks_1samp(sample_flares, uniform.cdf, args=(0, 1))
            A, p_ad_samp = ad_test(sample_flares, uniform(0,1), assume_sorted=True)
            KS_Sample_list.append(p_ks_samp)
            AD_Sample_list.append(p_ad_samp)
    elif error == False:
        for samples in range(N):
            sample_flares = flares
            check = np.random.random()
            if check > 0.5:
                sample_flares = sample_flares + np.abs(np.random.normal(0, scale=0.1))
            elif check < 0.5:
                sample_flares = sample_flares - np.abs(np.random.normal(0, scale=0.1))
            for index in range(len(sample_flares)):
                if sample_flares[index] > 1:
                    sample_flares[index] = sample_flares[index] - 1 
                if flares[index] < 0:
                    sample_flares[index] = sample_flares[index] + 1
            sample_flares = np.sort(sample_flares)
            D, p_ks_samp = ks_1samp(sample_flares, uniform.cdf, args=(0, 1))
            A, p_ad_samp = ad_test(sample_flares, uniform(0,1), assume_sorted=True)
            KS_Sample_list.append(p_ks_samp)
            AD_Sample_list.append(p_ad_samp)
    return np.median(KS_Sample_list), np.median(AD_Sample_list)
    
def K_Tests_Total(flare_df, output_dir, TOI = False, CDFs = False, sample = 10):
    KS_Sample_list = []
    AD_Sample_list = []
    N = sample

    # interest_hosts = [['Host_ID', 'p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N', 'dlog(Z)', 'Periastron?', 'Obs_Time', 'Period', 'Close_Approach']]
    interest_hosts = [['Host_ID', 'p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N', 'dlog(Z)', 'Periastron?']]
    All_CDF = []
    set_hosts = set(flare_df['Host_ID'])
    for index, hosts in enumerate(set_hosts):
        print(hosts)
        phases = np.array(flare_df.loc[flare_df['Host_ID'] == hosts, 'Transit_Phase'], dtype=float)
        dlogz = np.median(np.array(flare_df.loc[flare_df['Host_ID'] == hosts, 'dlogZ']))
        if TOI == False:
            peri_phases = np.array(flare_df.loc[flare_df['Host_ID'] == hosts, 'Periastron_Phase'])
            peri_phases_lower = np.abs(np.array(flare_df.loc[flare_df['Host_ID'] == hosts, 'Periastron_Lower'])[0])
            peri_phases_upper = np.array(flare_df.loc[flare_df['Host_ID'] == hosts, 'Periastron_Upper'])[0]
            if len(phases) < 5 or len(peri_phases) < 5:
                continue
            if np.isnan(np.mean(peri_phases)) == False:
                peri = True
            if np.isnan(np.mean(peri_phases)) == True:
                peri = False
            if np.isnan(np.mean(phases)) == True and peri == False:
                continue
            if np.isnan(peri_phases_lower) == True or np.isnan(peri_phases_upper) == True:
                error = False
            if np.isnan(peri_phases_lower) == False and np.isnan(peri_phases_upper) == False:
                error = True
        elif TOI == True:
            peri = False
        if peri == True:
            phases = peri_phases
            if len(phases) < 3:
                continue
            phases = np.sort(phases)
            D, p_KS = ks_1samp(phases, uniform.cdf, args=(0, 1))
            A, p_AD = ad_test(phases, uniform(0,1), assume_sorted=True)
            a, p_KU = kuiper(phases, uniform.cdf, args=(0, 1))
            for phase_sample in range(N):
                if error == False:
                    phases = phases + np.random.normal(0, scale=0.05)
                elif error == True:
                    check = np.random.random()
                    if check > 0.5:
                        phases = phases + np.abs(np.random.normal(0, scale=np.abs(peri_phases_upper)))
                    elif check < 0.5:
                        phases = phases - np.abs(np.random.normal(0, scale=np.abs(peri_phases_lower)))
                    for index in range(len(phases)):
                        if phases[index] > 1:
                            phases[index] = phases[index] - 1 
                        if phases[index] < 0:
                            phases[index] = phases[index] + 1
                    phases = np.sort(phases)
                    D, p_ks_samp = ks_1samp(phases, uniform.cdf, args=(0, 1))
                    A, p_ad_samp = ad_test(phases, uniform(0,1), assume_sorted=True)
                    KS_Sample_list.append(p_ks_samp)
                    AD_Sample_list.append(p_ad_samp)
        if peri == False:
            if len(phases) < 3:
                continue
            phases = np.sort(phases)
            D, p_KS = ks_1samp(phases, uniform.cdf, args=(0, 1))
            A, p_AD = ad_test(phases, uniform(0,1))
            a, p_KU = kuiper(phases, uniform.cdf, args=(0, 1))
            for phase_sample in range(N):
                phases += 1/N
                for index in range(len(phases)):
                    if phases[index] > 1:
                        phases[index] = phases[index] - 1 
                    if phases[index] < 0:
                        phases[index] = phases[index] + 1
                phases = np.sort(phases)
                D, p_ks_samp = ks_1samp(phases, uniform.cdf, args=(0, 1))
                A, p_ad_samp = ad_test(phases, uniform(0,1), assume_sorted=True)
                KS_Sample_list.append(p_ks_samp)
                AD_Sample_list.append(p_ad_samp)
        if p_KU < 1:
            interest_hosts.append([hosts, p_KU, p_KS, p_AD, np.median(KS_Sample_list), np.median(AD_Sample_list), len(phases), np.mean(dlogz), peri])
            x = (np.sort(phases)).tolist()
            y = (np.arange(len(x))/float(len(x))).tolist()
            x.insert(0, str(hosts) + '_x')
            y.insert(0, str(hosts) + '_y')
            x.append(1)
            y.append(1)
            All_CDF.append(x)
            
            KS_Sample_list = []
            AD_Sample_list = []
        if CDFs == True:
            with open(output_dir + "/Host_CDFs.csv","w+", newline='') as f:
                writer = csv.writer(f)
                for values in zip_longest(*All_CDF):
                    writer.writerow(values)
        with open(output_dir + "/K_Tests.csv","w+", newline='') as f:
            writer = csv.writer(f)
            for values in interest_hosts:
                writer.writerow(values)

def K_Tests(flares, periods, epoch, KS = True, KU = True, AD = True, sampling = True,
            N = 10, sample_type = ['periastron'], peri_error = None):
    '''
    Generates test statistics for different GoF tests to test if the provided
    flare epoch distribution deviates significantly from a uniform distribution.

    Parameters
    ----------
    flares : array-like, float
        Flare epochs, in BJD.
    period : array-like, float
        Period(s) of the planet(s).
    epoch : array-like, float
        Transit or periastron epoch(s), in BJD. Should be same length as 
        period array.
    KS : Bool, optional
        Set to true if you wish to compute the 1-Sample Kolmolgorov-Smirnov
        test statistic. The default is True.
    KU : Bool, optional
        Set to true if you wish to compute the 1-Sample Kuiper
        test statistic. The default is True.. The default is True.
    AD : Bool, optional
        Set to true if you wish to compute the 1-Sample Anderson-Darling
        test statistic. The default is True.
    sampling : Bool, optional
        Set to true if you want to return the average test statistic
        for either the AD or KS tests sampled at different starting
        phases. Will add an additional test statistic for each
        test sampled. The default is True.
    N: int
        If sampling is true, the number of samples to draw.
    sample_type: array-like, string ('transit' or 'periastron')
        List of whether to sample assuming a transit (uniform sampling) or assuming
        a known argument of periastron (Gaussian).
    peri_error: 2D array, floats
        The uncertainty in the epoch of periastron for sampling, in days.
        If symmetric, provide an array of length 1. If asymmetric, provide 
        as 2 element array in the format [lower,upper].

    Returns
    -------
    ND Array: numpy array
        Returns an ND array returning the test statistics

    '''
    K_tests = []
    message_str = 'Output Order\n'
    idx_counter = 0
    for index, period in enumerate(periods):
        phases = ((flares - (epoch[index] + period/2)) % period)/period
        print(phases)
        phases = np.sort(phases)
        if KS == True:
            idx_counter += 1
            D, p_KS = ks_1samp(phases, uniform.cdf, args=(0, 1))
            K_tests.append(p_KS)
            message_str += str(idx_counter) + ' KS\n'
        if KU == True:
            idx_counter += 1
            a, p_KU = kuiper(phases, uniform.cdf, args=(0, 1))
            K_tests.append(p_KU)
            message_str += str(idx_counter) + ' KU\n'
        if AD == True:
            idx_counter += 1
            A, p_AD = ad_test(phases, uniform(0,1), assume_sorted=True)
            K_tests.append(p_AD)
            message_str += str(idx_counter) + ' AD\n'
        if sampling == True:
            p_AD_Samp = []
            p_KS_Samp = []
            phase1 = phases
            for samples in range(N):
                if sample_type[index].lower() == 'transit':
                    check = np.random.random()
                    if check < 0.5:
                        phase1 += np.random.random()
                    elif check > 0.5:
                        phase1 -= np.random.random()
                if sample_type[index].lower() == 'periastron':
                    if peri_error == None:
                        peri_error = [0.1]
                    if len(peri_error) == 1:
                        phase1 += np.random.normal(scale=peri_error[0])
                    elif len(peri_error) == 2:
                        check = np.random.random()
                        if check < 0.5:
                            phase1 -= np.random.normal(scale=peri_error[0])
                        elif check > 0.5:
                            phase1 += np.random.normal(scale=peri_error[1])
                
                for idx, samps in enumerate(phase1):
                    if phase1[idx] > 1:
                        phase1[idx] -= 1
                    if phase1[idx] < 0:
                        phase1[idx] += 1
                phase1 = np.sort(phase1)
                A, p_ad_samp = ad_test(phase1, uniform(0,1), assume_sorted=True)
                D, p_ks_samp = ks_1samp(phase1, uniform.cdf, args=(0, 1))
                p_AD_Samp.append(p_ad_samp)
                p_KS_Samp.append(p_ks_samp)
            if KS == True:
                idx_counter += 1
                K_tests.append(np.median(p_KS_Samp))
                message_str += str(idx_counter) + ' KS Sampled\n'
            if AD == True:
                idx_counter += 1
                K_tests.append(np.median(p_AD_Samp))
                message_str += str(idx_counter) + ' AD Sampled\n'
        print(message_str)
    return K_tests

data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")

flares = data.loc[data['Host_ID'] == 'YZCet', 'Flare_Epoch']

K_Tests(flares, [3.8888], [2457010], N=50, sample_type=['periastron'])
    