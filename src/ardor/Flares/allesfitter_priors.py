# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:48:26 2023

@author: Nathan
"""
#%%
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import simpson
import numpy as np
import os
from ardor.Flares import aflare
import math
from ardor.Utils.planck_law import planck_law, planck_integrator
from ardor.Utils.Utils import asymmetric_sample, find_nearest, copy_output
import allesfitter
import astropy.io.ascii as ascii
from astropy.table import Table
import shutil
import pandas as pd
from matplotlib import pyplot as plt
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

def clear_workingdir(working_dir):

    files = os.listdir(working_dir)
    for file_name in files:
        if os.path.isfile(os.path.join(working_dir, file_name)):
            file_path = os.path.join(working_dir, file_name)
            os.remove(file_path)

    ## Clear the files in the results folder
    files = os.listdir(os.path.join(working_dir, 'results'))
    for file_name in files:
        file_path = os.path.join(working_dir, 'results', file_name)
        os.remove(file_path)

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


def slice_flare(time, data):
    '''
    Slices flare data to only include data around the flare event.

    Parameters
    ----------
    time : numpy array
        Time axis data.
    data : numpy array
        Flux axis data.

    Returns
    -------
    time_slice : numpy array
        Sliced time axis data.
    data_slice : numpy array
        Sliced flux axis data.

    '''
    index = find_nearest(time, 0)[1]
    return time[index:], data[index:]

def find_local_maxima(time, data):
    '''
    Finds local maxima and minima in the flare data.

    Parameters
    ----------
    time : numpy array
        Time axis data.
    data : numpy array
        Flux axis data.

    Returns
    -------
    peak_times : list
        List of times where local maxima occur.
    min_times : list
        List of times where local minima occur.
    '''
    smoothed_data = gaussian_filter1d(data, sigma=1.1)
    peak_times = [0]
    min_times = []
    for i in range(1, len(smoothed_data)-1):
        if len(peak_times) > 3:
            break
        if smoothed_data[i] > smoothed_data[i-1] and smoothed_data[i] > smoothed_data[i+1]:
            peak_times.append(time[i])
        elif smoothed_data[i] < smoothed_data[i-1] and smoothed_data[i] < smoothed_data[i+1]:
            min_times.append(time[i])
    if len(peak_times) > 3:
        peak_times = peak_times[:3]
    if len(min_times) > 3:
        min_times = min_times[:3]
    return peak_times, min_times
def return_flare_time_priors(peak_times, min_times):
    priors = []
    for idx, t in enumerate(peak_times):
        if idx == 0:
            prior_lower = -0.02
            prior_upper = 0 + min_times[0]
        elif idx != 0:
            try:
                prior_lower = min_times[idx - 1]
                prior_upper = min_times[idx]
            except IndexError:
                prior_upper = t + 0.002
                prior_lower = min_times[-1]
        priors.append([prior_lower, prior_upper])
    return priors

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
def multipeak_model_compare(target_file, working_dirs, template_dir, baseline='hybrid_spline', N_models = 3, dlogZ_threshold = 2):
    """Uses Nested Sampling to compare different number of peaks in a flare model. Begins comparing no flare
    (only noise) to a single peaked flare, then adds additional peak(s) if the Bayesian evidence supports adding
    another peak.

    Args:
        target_file (str): Path to flare csv file.
        working_dirs (str): Path to working directories for allesfitter runs.
        template_dir (str): Path to template directory for allesfitter settings.csv and parameters.csv.
        baseline (str, optional): Which allesfitter baseline to use. Defaults to 'hybrid_spline'.
        N_models (int, optional): Number of peaks to compare up to. Defaults to 3.
        dlogZ_threshold (int, optional): Threshold for change in log Bayesian evidence to add another peak. Defaults to 2.
    Raises:
        ValueError: You should have at least as many working directories as models + 1 (for baseline model).

    Returns:
        float: Change in log Bayesian evidence of the best .
    """
    if len(working_dirs) != N_models + 1:
        raise ValueError("Working dirs length does not match number of models.")
    csv_cleaner(target_file)
    file_num = os.path.basename(target_file).split('.csv')[0]
    slice_time, slice_data = slice_flare(np.array(pd.read_csv(target_file, header=None)[0]), np.array(pd.read_csv(target_file, header=None)[1]))
    peak_times, min_times = find_local_maxima(slice_time, slice_data)
    peak_time_best_guess = peak_times
    peak_time_priors = return_flare_time_priors(peak_times, min_times)
    log_Z0 = 0
    log_Z1 = 0
    log_Z2 = 0
    log_Z3 = 0
    dlogZ = 1000
    ## Copy the csv to each working directory
    for idx, working_dir in enumerate(working_dirs):
        ## Determine Bayes Factor for baseline model
        if idx == 0:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num, 
                                 peak_time_best_guess=peak_time_best_guess, peak_time_priors=peak_time_priors)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z0 = allesfitter.get_logZ(working_dir)[0][0]
        ## Determine Bayes Factor for N=1 model
        if idx == 1:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num, 
                                 peak_time_best_guess=peak_time_best_guess, peak_time_priors=peak_time_priors)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z1 = allesfitter.get_logZ(working_dir)[0][0]
            ## If change in Bayes Factor < 5, exit and return dlogZ
            dlogZ = log_Z1 - log_Z0
            if dlogZ < dlogZ_threshold:
                return print("Data is best modeled by baseline/no flare model.")
            if dlogZ > dlogZ_threshold and N_models == 1:
                print("Data is best modeled by 1 flare peak model.")
                return log_Z1 - log_Z0
        ## If change in Bayes Factor >= 5, continue to next model
        if idx == 2:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num, 
                                 peak_time_best_guess=peak_time_best_guess, peak_time_priors=peak_time_priors)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z2 = allesfitter.get_logZ(working_dir)[0][0]
            dlogZ = log_Z2 - log_Z1
            if dlogZ < dlogZ_threshold:
                print("Data is best modeled by 1 flare peak model.")
                return log_Z1 - log_Z0
            if dlogZ >= dlogZ_threshold and N_models == 2:
                print("Data is best modeled by 2 flare peak model.")
                return log_Z2 - log_Z1
        if idx == 3:
            construct_param_file(os.path.join(working_dir, 'params.csv'), baseline=baseline, N=idx, name=file_num, 
                                 peak_time_best_guess=peak_time_best_guess, peak_time_priors=peak_time_priors)
            construct_settings_file(os.path.join(template_dir, 'settings.csv'), working_dir, baseline=baseline, N=idx, name=file_num)
            shutil.copyfile(target_file, os.path.join(working_dir, os.path.basename(target_file)))
            allesfitter.ns_fit(working_dir)
            allesfitter.ns_output(working_dir)
            log_Z3 = allesfitter.get_logZ(working_dir)[0][0]
            dlogZ = log_Z3 - log_Z2
            if dlogZ < dlogZ_threshold:
                print("Data is best modeled by 2 flare peak model.")
                return log_Z2 - log_Z1
            elif dlogZ >= dlogZ_threshold:
                print("Data is best modeled by 3 flare peak model.")
                return log_Z3 - log_Z2

### Constructing parameter and settings files for allesfitter
def construct_param_file(output_dir,  peak_time_best_guess = None, peak_time_priors = None, baseline = 'hybrid_spline', 
                         model = 'flare', N=1, name = 'Flare', custom_priors = False, priors = None):
    if len(peak_time_best_guess) > 3 and N == 3:
        peak_time_best_guess = peak_time_best_guess[:3]
        peak_time_priors = peak_time_priors[:3]
        print("Number of peak time best guesses exceeds 3. Only first 3 will be used.")
    if custom_priors == False:
        if model == 'flare':
            if N == 0:
                param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1']
                best_guess = [-0.02, 0.01, 0.1]
                priors = [f'uniform -0.05 {peak_time_priors[0][1]}', 'uniform 0 0.1', 'uniform 0 3']
                labels = ['Flare_Time', 'Flare_FWHM', 'Flare_Amp.']
                units = ['days', 'days', 'rel. flux']
            if N == 1:
                param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1']
                best_guess = [peak_time_best_guess[0], 0.01, 0.1]
                priors = [f'uniform -0.01 {peak_time_priors[0][1]}', 'uniform 0 0.1', 'uniform 0 3']
                labels = ['Flare_Time', 'Flare_FWHM', 'Flare_Amp.']
                units = ['days', 'days', 'rel. flux']
            elif N == 2:
                param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                            'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2']
                best_guess = [peak_time_best_guess[0], 0.01, 0.1, peak_time_best_guess[1], 0.01, 0.1]
                priors = [f'uniform {peak_time_priors[0][0]} {peak_time_priors[0][1]}', 'uniform 0 0.1', 'uniform 0 3', 
                            f'uniform {peak_time_priors[1][0]} {peak_time_priors[1][1]}', 'uniform 0 0.1', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.']
                units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux']
            elif N == 3:
                param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                            'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2',
                            'flare_tpeak_3', 'flare_fwhm_3', 'flare_ampl_3']
                best_guess = [peak_time_best_guess[0], 0.01, 0.1, peak_time_best_guess[1], 0.01, 0.1, peak_time_best_guess[2], 0.01, 0.1]
                priors = [f'uniform {peak_time_priors[0][0]} {peak_time_priors[0][1]}', 'uniform 0 0.5', 'uniform 0 3', 
                            f'uniform {peak_time_priors[1][0]} {peak_time_priors[1][1]}', 'uniform 0 0.5', 'uniform 0 3',
                            f'uniform {peak_time_priors[2][0]} {peak_time_priors[2][1]}', 'uniform 0 0.5', 'uniform 0 3']
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

def save_params_to_csv(params, Teff, Rst, dlogZ, output_csv):
    """_summary_

    Args:
        params (_type_): _description_
        Teff (_type_): _description_
        Rst (_type_): _description_
        dlogZ (_type_): _description_
        output_csv (_type_): _description_
    """
    if len(params) == 1:
        energy = flare_energy([ params[0]['epoch'][0], params[0]['fwhm'][0],  params[0]['amplitude'][0]], Teff[0], Rst[0],
                                               [ params[0]['epoch'][1]],[ params[0]['fwhm'][1],  params[0]['amplitude'][1]],
                                               [ params[0]['epoch'][2]],[ params[0]['fwhm'][2],  params[0]['amplitude'][2]], Teff[1], Teff[2],
                                               Rst[1], Rst[2], N = len(flares), N_samp = 1000)
    if len(params) == 2:
        energy = flare_energy([params[0]['epoch'][0], params[0]['fwhm'][0],  params[0]['amplitude'][0], params[1]['epoch'][0], 
                               params[1]['fwhm'][0],  params[1]['amplitude'][0]], Teff[0], Rst[0],
                                        [params[0]['epoch'][1],params[0]['fwhm'][1],  params[0]['amplitude'][1], 
                                        params[1]['epoch'][1],params[1]['fwhm'][1],  params[1]['amplitude'][1]],
                                        [params[0]['epoch'][2], params[0]['fwhm'][2],  params[0]['amplitude'][2], 
                                         params[1]['epoch'][2], params[1]['fwhm'][2],  params[1]['amplitude'][2]], 
                                        Teff[1], Teff[2],Rst[1], Rst[2], N = len(flares), N_samp = 1000)
    if len(params) == 3:
        energy = flare_energy([params[0]['epoch'][0], params[0]['fwhm'][0],  params[0]['amplitude'][0], params[1]['epoch'][0], 
                               params[1]['fwhm'][0],  params[1]['amplitude'][0], params[2]['epoch'][0], 
                               params[2]['fwhm'][0],  params[2]['amplitude'][0]], Teff[0], Rst[0],
                                        [params[0]['epoch'][1],params[0]['fwhm'][1],  params[0]['amplitude'][1], 
                                        params[1]['epoch'][1],params[1]['fwhm'][1],  params[1]['amplitude'][1],
                                        params[2]['epoch'][1],params[2]['fwhm'][1],  params[2]['amplitude'][1]],
                                        [params[0]['epoch'][2], params[0]['fwhm'][2],  params[0]['amplitude'][2], 
                                         params[1]['epoch'][2], params[1]['fwhm'][2],  params[1]['amplitude'][2],
                                         params[2]['epoch'][2], params[2]['fwhm'][2],  params[2]['amplitude'][2]], 
                                        Teff[1], Teff[2],Rst[1], Rst[2], N = len(flares), N_samp = 1000)
    for flares in params:
        if output_csv != None:
            os.makedirs(output_csv, exist_ok=True)
            data = pd.read_csv(os.path.join(output_csv, 'Tier_3_Flare_Params.csv'), header=None)
            if data.empty:
                with open(output_csv + '/Tier_3_Flare_Params.csv', "a") as f:
                    f.write('Epoch_Best,Epoch,Epoch_Unc_L,Epoch_Unc_U,Amplitude_Best,Amplitude_Unc_L,Amplitude_Unc_U,FWHM_Best,FWHM_Unc_L,FWHM_Unc_U,Energy_Best,Energy_Unc_L,Energy_Unc_U,dlogZ\n')
            ZZ = np.column_stack((flares['epoch'][0], flares['epoch'][0] + 2457000, flares['epoch'][1], flares['epoch'][2],
                                  flares['amplitude'][0], flares['amplitude'][1], flares['amplitude'][2],
                                  flares['fwhm']['0'], flares['fwhm'][1], flares['fwhm'][2], energy[0], energy[1], energy[2],
                                  dlogZ))
            with open(output_csv + '/Tier_3_Flare_Params.csv', "a") as f:
                np.savetxt(f, ZZ, delimiter=",", fmt='%s')


working_dirs = ["/ugrad/whitsett.n/ardor_test/Working_Dirs/Baseline", "/ugrad/whitsett.n/ardor_test/Working_Dirs/1Flare", "/ugrad/whitsett.n/ardor_test/Working_Dirs/2Flare",
                  "/ugrad/whitsett.n/ardor_test/Working_Dirs/3Flare"]

# multipeak_model_compare("/ugrad/whitsett.n/ardor_test/AUMic/Flare_56.csv", working_dirs, "/ugrad/whitsett.n/ardor/templates")
# clear_workingdir("/ugrad/whitsett.n/ardor_test/Working_Dirs/2Flare")
# clear_workingdir("/ugrad/whitsett.n/ardor_test/Working_Dirs/Baseline")
# clear_workingdir("/ugrad/whitsett.n/ardor_test/Working_Dirs/1Flare")
# %%
