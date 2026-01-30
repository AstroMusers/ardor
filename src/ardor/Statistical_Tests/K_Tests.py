# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:42:31 2024

@author: whitsett.n
"""

import numpy as np
import pandas as pd
from scipy.stats import uniform
from scipy.stats import ks_1samp
from skgof import ad_test
import csv
from itertools import zip_longest
from astropy.stats import kuiper
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
    '''
    Computes KU, AD, KS tests for simulations found in SPI_Simulation.py.

    Parameters
    ----------
    DataFrame : pd.df
        Dataframe of the simulates flares.
    output_dir : str
        Directory for the output.
    descriptor : str, optional
        What to name the outputfile. The default is 'G'.

    Returns
    -------
    .csv file of the output of the tests.

    '''
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
    

def K_Tests(phases, KS = True, KU = True, AD = True, sampling = True,
            N = 10, output_message = True):
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
    sample_type: array-like, int
        The epoch type:
            0: The target has no constraint on any epoch.
            1: The target has a constraint on the epoch of periastron
            2: The target has a constraint on transit epoch, but not periastron
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
    phases = np.sort(phases)
    K_sampled = []
    if len(phases) >= 3 and np.isnan(phases[0]) == False:
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
            phase_offset = np.linspace(0, 1, N+1)
            KS_samples = []
            AD_samples = []
            for offset in phase_offset:
                phases_sampled = (phases + offset) % 1
                phases_sampled = np.sort(phases_sampled)
                if KS == True:
                    D, p_KS = ks_1samp(phases_sampled, uniform.cdf, args=(0, 1))
                    KS_samples.append(p_KS)
                    message_str += str(idx_counter) + ' KS Sampled\n'
                if AD == True:
                    idx_counter += 1
                    A, p_AD = ad_test(phases_sampled, uniform(0,1), assume_sorted=True)
                    AD_samples.append(p_AD)
                    message_str += str(idx_counter) + ' AD Sampled\n'
            if KS == True:
                K_sampled.append(np.median(KS_samples))
            if AD == True:
                K_sampled.append(np.median(AD_samples))
    if output_message == True:
        print(message_str)
    if sampling == True:
        return K_tests, K_sampled
    elif sampling == False:
        return K_tests
        