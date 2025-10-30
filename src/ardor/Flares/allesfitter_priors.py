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
from astropy.table import Table
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

def run_ns(target_file, working_dir, param_template_dir, settings_template_dir, baseline = False):
    #Check each folder
    # Copy relevant file to allesfitter folder
    name = os.path.basename(os.path.normpath(target_file))
    csv_cleaner(working_dir + '/' + name)
    allesfitter.ns_fit(working_dir)
    allesfitter.ns_output(working_dir)
    os.remove(working_dir + '/' + name)
    return target_file, working_dir, param_template_dir, settings_template_dir

def model_compare(target_file, model1_working_dir, model2_working_dir, model1_param_template_dir, model2_param_template_dir, model1_settings_template_dir, model2_settings_template_dir):
    run_ns(target_file, model1_working_dir, model1_param_template_dir, model1_settings_template_dir, baseline=True)
    run_ns(target_file, model2_working_dir, model2_param_template_dir, model2_settings_template_dir, baseline=False)
    model1_logz = allesfitter.get_logZ(model1_working_dir)
    model2_logz = allesfitter.get_logZ(model2_working_dir)
    d_logz = model2_logz[0][0] - model1_logz[0][0]
    for files in os.listdir(model1_working_dir + '/results'):
        os.remove(model1_working_dir + '/results/' + files)
    for files in os.listdir(model2_working_dir + '/results'):
        os.remove(model2_working_dir + '/results/' + files)
    return d_logz

def construct_param_file(param_dir, priors = None, baseline = 'hybrid_spline', 
                         model = 'flare', N=1, name = 'Flare'):
    if model == 'flare':
        if N == 1:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1']
            if priors == None:
                best_guess = [0, 0.05, 0.1]
                priors = ['uniform -1 1', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare_Time', 'Flare_FWHM', 'Flare_Amp.']
            units = ['days', 'days', 'rel. flux']
        elif N == 2:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                           'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2']
            if priors == None:
                best_guess = [0, 0.05, 0.1, 0.2, 0.05, 0.1]
                priors = ['uniform -1 1', 'uniform 0 0.5', 'uniform 0 3', 
                          'uniform -1 1', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.']
            units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux']
        elif N == 3:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                           'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2',
                           'flare_tpeak_3', 'flare_fwhm_3', 'flare_ampl_3']
            if priors == None:
                best_guess = [0, 0.05, 0.1, 0.2, 0.05, 0.1, 0.4, 0.05, 0.1]
                priors = ['uniform -0.5 0.5', 'uniform 0 0.5', 'uniform 0 3', 
                          'uniform 0 0.5', 'uniform 0 0.5', 'uniform 0 3',
                          'uniform 0 0.5', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.',
                            'Flare3_Time', 'Flare3_FWHM', 'Flare3_Amp.']
            units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux',
                        'days', 'days', 'rel. flux']
    if baseline == 'gp_matern32':
        param_names += [f'ln_err_flux_{name}', f'baseline_gp_offset_flux_{name}', f'baseline_gp_matern32_lnsigma_flux_{name}', 
                        f'baseline_gp_matern32_lnrho_flux_{name}']
        best_guess += [-7,0, -5, 0]
        priors += ['uniform -10 -2', 'uniform -0.02, 0.02', 'uniform -15 0', 'uniform -1 15']
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
    table.add_column(param_names, name='#param_name')
    table.add_column(best_guess, name='best_guess')
    table.add_column([1]*len(param_names), name='fit')
    table.add_column(priors, name='prior')
    table.add_column(labels, name='label')
    table.add_column(units, name='unit')
    ascii.write(table, os.path.join(param_dir, 'params.csv'), overwrite=True, format='csv')

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