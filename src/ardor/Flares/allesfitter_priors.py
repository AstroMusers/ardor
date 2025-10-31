# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:48:26 2023

@author: Nathan
"""

import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import simpson
import numpy as np
import os
from ardor.Flares import aflare
import math
from ardor.Utils.planck_law import planck_law, planck_integrator
from ardor.Utils.Utils import asymmetric_sample
import allesfitter
import astropy.io.ascii as ascii
from astropy.table import Table
import shutil
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

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

def flare_energy(params, Teff, R_stellar, uncertainty = False, param_lower = [0,0,0], param_upper = [0,0,0], Teff_l = 0, Teff_u = 0,
                 R_stellar_l = 0, R_stellar_u = 0, N = 1, N_samp = 1000):
    if len(params) != 3*N:
        raise ValueError("Parameter length does not match number of flares.")
    for flares in range(N):
        if uncertainty == False:
            x = np.linspace(0, 0.02, num = 2000)
            y = aflare.aflare(x, params)
            flare_area = simpson(y, x)
            color_factor = planck_integrator(600e-6, 1000e-6, Teff)/planck_integrator(600e-6, 1000e-6, 9000)
            energy = (5.67e-8)*(9000**4)*(flare_area)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*color_factor*(1e7)*86400
            return energy
        if uncertainty == True:
            x = np.linspace(0, 0.02, num = 2000)
            y = aflare.aflare(x, params)
            flare_area = simpson(y, x)
            color_factor = planck_integrator(600e-6, 1000e-6, Teff)/planck_integrator(600e-6, 1000e-6, 9000)
            energy = (5.67e-8)*(9000**4)*(flare_area)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*color_factor*(1e7)*86400
            samples = []
            for flare_samples in range(N_samp):
                if N == 1:
                    R_stellar_sample = asymmetric_sample(R_stellar, R_stellar_u, R_stellar_l)
                    Teff_sample = asymmetric_sample(Teff, Teff_u, Teff_l)
                    fwhm_sample = asymmetric_sample(params[1], param_upper[1], param_lower[1])
                    ampl_sample = asymmetric_sample(params[2], param_upper[2], param_lower[2])
                    y = aflare.aflare(x, [params[0], fwhm_sample, ampl_sample])
                    flare_temp = np.random.normal(loc=9000, scale=500)
                    color_factor = planck_integrator(600e-6, 1000e-6, Teff_sample)/planck_integrator(600e-6, 1000e-6, flare_temp)
                    energy_sample = (5.67e-8)*(flare_temp**4)*(flare_area)*np.pi*(R_stellar_sample*6.957e8*R_stellar_sample*6.957e8)*color_factor*(1e7)*86400
                    samples.append(energy_sample)
                if N == 2:
                    R_stellar_sample = asymmetric_sample(R_stellar, R_stellar_u, R_stellar_l)
                    Teff_sample = asymmetric_sample(Teff, Teff_u, Teff_l)
                    fwhm1_sample = asymmetric_sample(params[1], param_upper[1], param_lower[1])
                    ampl1_sample = asymmetric_sample(params[2], param_upper[2], param_lower[2])
                    fwhm2_sample = asymmetric_sample(params[4], param_upper[4], param_lower[4])
                    ampl2_sample = asymmetric_sample(params[5], param_upper[5], param_lower[5])
                    y = aflare.aflare(x, [params[0], fwhm1_sample, ampl1_sample, params[3], fwhm2_sample, ampl2_sample])
                    flare_temp = np.random.normal(loc=9000, scale=500)
                    color_factor = planck_integrator(600e-6, 1000e-6, Teff_sample)/planck_integrator(600e-6, 1000e-6, flare_temp)
                    energy_sample = (5.67e-8)*(flare_temp**4)*(flare_area)*np.pi*(R_stellar_sample*6.957e8*R_stellar_sample*6.957e8)*color_factor*(1e7)*86400
                    samples.append(energy_sample)
                if N == 3:
                    R_stellar_sample = asymmetric_sample(R_stellar, R_stellar_u, R_stellar_l)
                    Teff_sample = asymmetric_sample(Teff, Teff_u, Teff_l)
                    fwhm1_sample = asymmetric_sample(params[1], param_upper[1], param_lower[1])
                    ampl1_sample = asymmetric_sample(params[2], param_upper[2], param_lower[2])
                    fwhm2_sample = asymmetric_sample(params[4], param_upper[4], param_lower[4])
                    ampl2_sample = asymmetric_sample(params[5], param_upper[5], param_lower[5])
                    fwhm3_sample = asymmetric_sample(params[7], param_upper[7], param_lower[7])
                    ampl3_sample = asymmetric_sample(params[8], param_upper[8], param_lower[8])
                    y = aflare.aflare(x, [params[0], fwhm1_sample, ampl1_sample, params[3], fwhm2_sample, ampl2_sample, params[6], fwhm3_sample, ampl3_sample])
                    flare_temp = np.random.normal(loc=9000, scale=500)
                    color_factor = planck_integrator(600e-6, 1000e-6, Teff_sample)/planck_integrator(600e-6, 1000e-6, flare_temp)
                    energy_sample = (5.67e-8)*(flare_temp**4)*(flare_area)*np.pi*(R_stellar_sample*6.957e8*R_stellar_sample*6.957e8)*color_factor*(1e7)*86400
                    samples.append(energy_sample)
            samples = sorted(samples)
            return energy, -np.abs(energy - samples[int(N_samp*0.32)]), np.abs(energy - samples[int(N_samp*0.68)])

def csv_cleaner(flare_csv_dir):
    data = pd.read_csv(flare_csv_dir, header=None, index_col=False)
    time = np.array(data[0])
    flux = np.array(data[1])
    error = np.array(data[2])
    gap_index = 0
    error_index = []
    for index in range(len(time) - 1):
        dt = time[index + 1] - time[index]
        if dt > 10:
            gap_index = index + 1
    
    for index in range(len(time)):
        if np.isnan(error[index]) == True:
            error[index] = 1e-3
            error_index.append(index)
    for index in error_index:
        error[index] = np.average(error)
        if np.isnan(np.average(error)) == True:
            error[index] = 1e-3
    if gap_index != 0 and gap_index < len(time)/2:
        time = time[gap_index:]
        flux = flux[gap_index:]
        error = error[gap_index:]
    elif gap_index != 0 and gap_index >= len(time)/2:
        time = time[:gap_index]
        flux = flux[:gap_index]
        error = error[:gap_index]
    output = np.stack((time,flux,error)).T
    np.savetxt(flare_csv_dir, output, delimiter=',')

def multipeak_model_compare(target_file, working_dirs, template_dir, baseline='hybrid_spline', N_models = 3):
    if len(working_dirs) != N_models + 1:
        raise ValueError("Working dirs length does not match number of models.")
    csv_cleaner(target_file)
    file_num = os.path.basename(target_file).split('.csv')[0]
    log_Z1 = 0
    log_Z2 = 0
    log_Z3 = 0
    log_Z4 = 0
    dlogZ = 1000
    ## Copy the csv to each working directory
    for idx, working_dir in enumerate(working_dirs):
        ## Determine Bayes Factor for baseline model
        if idx == 0:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z1 = allesfitter.get_logZ(working_dir)
        ## Determine Bayes Factor for N=1 model
        if idx == 1:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z2 = allesfitter.get_logZ(working_dir)
            ## If change in Bayes Factor < 5, exit and return dlogZ
            dlogZ = log_Z2[0][0] - log_Z1[0][0]
            if dlogZ < 5:
                return dlogZ
        ## If change in Bayes Factor >= 5, continue to next model
        if idx == 2:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z3 = allesfitter.get_logZ(working_dir)
            dlogZ = log_Z3[0][0] - log_Z2[0][0]
        if dlogZ < 5:
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            return dlogZ
        if idx == 3:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z4 = allesfitter.get_logZ(working_dir)
            dlogZ = log_Z4[0][0] - log_Z3[0][0]
            return dlogZ

def construct_param_file(output_dir, priors = None, baseline = 'hybrid_spline', 
                         model = 'flare', N=1, name = 'Flare'):
    if model == 'flare':
        if N == 0:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1']
            if priors == None:
                best_guess = [0, 0.01, 0.1]
                priors = ['uniform -1 1', 'uniform 0 0.1', 'uniform 0 3']
                labels = ['Flare_Time', 'Flare_FWHM', 'Flare_Amp.']
            units = ['days', 'days', 'rel. flux']
        if N == 1:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1']
            if priors == None:
                best_guess = [0, 0.01, 0.1]
                priors = ['uniform -0.01 0.0025', 'uniform 0 0.1', 'uniform 0 3']
                labels = ['Flare_Time', 'Flare_FWHM', 'Flare_Amp.']
            units = ['days', 'days', 'rel. flux']
        elif N == 2:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                           'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2']
            if priors == None:
                best_guess = [0, 0.01, 0.1, 0.0035, 0.01, 0.1]
                priors = ['uniform -0.01 0.0025', 'uniform 0 0.1', 'uniform 0 3', 
                          'uniform 0.0025 0.02', 'uniform 0 0.1', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.']
            units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux']
        elif N == 3:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                           'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2',
                           'flare_tpeak_3', 'flare_fwhm_3', 'flare_ampl_3']
            if priors == None:
                best_guess = [0, 0.05, 0.1, 0.01, 0.0035, 0.1, 0.02, 0.0075, 0.1]
                priors = ['uniform -0.01 0.0025', 'uniform 0 0.5', 'uniform 0 3', 
                          'uniform 0.0025 0.02', 'uniform 0 0.5', 'uniform 0 3',
                          'uniform 0.02 0.03', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.',
                            'Flare3_Time', 'Flare3_FWHM', 'Flare3_Amp.']
            units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux',
                        'days', 'days', 'rel. flux']
    if baseline == 'sample_GP_Matern32':
        param_names += [f'ln_err_flux_{name}', f'baseline_gp_offset_flux_{name}', f'baseline_gp_matern32_lnsigma_flux_{name}', 
                        f'baseline_gp_matern32_lnrho_flux_{name}']
        best_guess += [-7,0,-5, 0]
        priors += ['uniform -10 -2', 'uniform -0.02 0.02', 'uniform -15 0', 'uniform 0 15']
        labels += [f'ln_err_flux_{name}', rf'GP_Matern32_Baseline_{name}', rf'GP_Matern32_ln$sigma$_{name}', 
                   rf'GP_Matern32_ln$rho$_{name}']
        units += ['rel. flux', 'rel. flux', '', '']
    elif baseline == 'hybrid_spline':
        param_names += [f'ln_err_flux_{name}']
        best_guess += [-7]
        priors += ['uniform -10 -2']
        labels += [f'ln_err_flux_{name}']
        units += ['rel. flux']
    table = Table()
    table.add_column(param_names, name='#name')
    table.add_column(best_guess, name='value')
    table.add_column([1]*len(param_names), name='fit')
    table.add_column(priors, name='bounds')
    table.add_column(labels, name='label')
    table.add_column(units, name='unit')
    table.write(output_dir, overwrite=True, format='csv')

def construct_settings_file(settings_file_dir, working_dir, baseline = 'hybrid_spline', 
                         N=1, name = 'Flare', multi_process_bool = True
                         , cores = 10, walkers = 15, steps = 5000, burn_in = 1000):
    """Constructs settings file for Tier 3 allesfitter MCMC run.

    Args:
        settings_file_dir (str): Directory of template settings file.
        working_dir (str): Directory to save constructed settings file.
        baseline (str, optional): Baseline you want to use. Defaults to 'hybrid_spline'. 'sample_GP_Matern32' also available.
        N (int, optional): Number of flares to fit. Defaults to 1.
        name (str, optional): Name of the fit. Defaults to 'Flare'.
        multi_process_bool (bool, optional): Whether to use multiprocessing. Defaults to True.
        cores (int, optional): Number of cores to use for multiprocessing. Defaults to 10.
        walkers (int, optional): Number of walkers for MCMC. Defaults to 15.
        steps (int, optional): Total number of steps for MCMC. Defaults to 5000.
        burn_in (int, optional): Number of burn-in steps for MCMC. Defaults to 1000.
    """
    table = ascii.read(settings_file_dir, format='csv')
    for idx, values in enumerate(table['#name']):
        if 'Flaronardo' in table['#name'][idx]:
            table['#name'][idx] = str(values).replace("Flaronardo", name)
    for idx, values in enumerate(table['value']):
        if 'Flaronardo' in table['value'][idx]:
            table['value'][idx] = str(values).replace("Flaronardo", name)

    table['value'][table['#name'] == 'inst_phot'] = name
    table['value'][table['#name'] == 'multiprocess'] = multi_process_bool
    table['value'][table['#name'] == 'multiprocess_cores'] = cores
    table['value'][table['#name'] == 'mcmc_nwalkers'] = walkers
    table['value'][table['#name'] == 'mcmc_total_steps'] = steps
    table['value'][table['#name'] == 'mcmc_burn_steps'] = burn_in
    table['value'][table['#name'] == f'baseline_flux_{name}'] = baseline
    table['value'][table['#name'] == 'N_flares'] = N
    ascii.write(table, os.path.join(working_dir, 'settings.csv'), overwrite=True, format='csv')

def find_subsequent_peaks(time, data):
    """
    Find subsequent peaks in time series data after the main peak at t=0.
    
    Parameters
    ----------
    time : array-like
        Time array with the main peak guaranteed to be at time = 0
    data : array-like
        Flux data array corresponding to the time array
        
    Returns
    -------
    list
        List of times where subsequent peaks occur, only when data is at least
        1 standard deviation above the initial baseline flux
    """
    time = np.array(time)
    data = np.array(data)
    
    # Find the main peak at t=0 (or closest to 0)
    main_peak_idx = 0
    
    # Calculate baseline flux from data before the main peak
    # Use data points that are at least 10% of the time range before the main peak
    time_range = np.max(time) - np.min(time)
    baseline_mask = time < (time[main_peak_idx] - 0.1 * time_range)
    
    if np.sum(baseline_mask) < 2:
        # If not enough baseline points, use first 10% of data points
        n_baseline = max(2, len(time) // 10)
        baseline_flux = data[:n_baseline]
    else:
        baseline_flux = data[baseline_mask]
    
    # Calculate baseline statistics
    baseline_mean = np.mean(baseline_flux)
    baseline_std = np.std(baseline_flux)*3
    threshold = baseline_mean + baseline_std
    
    # Only look for peaks after the main peak
    after_peak_mask = time > time[main_peak_idx]
    time_after = time[after_peak_mask]
    data_after = data[after_peak_mask]
    
    # Only consider data points above the threshold
    above_threshold_mask = data_after >= threshold
    
    if np.sum(above_threshold_mask) < 3:
        # Not enough points above threshold to find peaks
        return []
    
    time_valid = time_after[above_threshold_mask]
    data_valid = data_after[above_threshold_mask]
    
    # Find local maxima
    peak_times = []
    
    # Need at least 3 points to find a local maximum
    if len(data_valid) >= 3:
        for i in range(1, len(data_valid) - 1):
            # Check if current point is a local maximum
            if (data_valid[i] > data_valid[i-1] and 
                data_valid[i] > data_valid[i+1]):
                peak_times.append(time_valid[i])
    
    return peak_times
# data = pd.read_csv("/ugrad/whitsett.n/ardor_test/AUMic/Flare_28.csv", header=None)
# time = data[0]
# flux = data[1]
# subsequent_peaks = find_subsequent_peaks(time, flux)
# print(subsequent_peaks)
data_dir = "/ugrad/whitsett.n/ardor_test/AUMic/Flare_28.csv"
working_dirs = ["/ugrad/whitsett.n/ardor_test/Working_Dirs/Baseline", "/ugrad/whitsett.n/ardor_test/Working_Dirs/1Flare", "/ugrad/whitsett.n/ardor_test/Working_Dirs/2Flare", "/ugrad/whitsett.n/ardor_test/Working_Dirs/3Flare"]
template = ""
multipeak_model_compare(data_dir, working_dirs, "/ugrad/whitsett.n/ardor/templates", baseline='hybrid_spline', N_models=3)