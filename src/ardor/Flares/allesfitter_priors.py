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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx
def flare_settings(data_file_dir, settings_template_dir, output_dir, flares=1, multi_process = False, cores = 1):
    settings = pd.read_csv(settings_template_dir, index_col=False)
    name = (data_file_dir.rsplit('/', 1)[1])[:-4]
    settings.at[5, 'value'] = name
    settings.at[10, 'value'] = multi_process
    if multi_process == True:
        settings.at[11, 'value'] = cores
    settings.at[33, '#name'] = 'baseline_flux_' + str(name)
    settings.at[37, '#name'] = 'error_flux_' + str(name)
    settings.at[41, 'value'] = flares
    settings.to_csv(output_dir, index=False)
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
def allesfitter_priors(flare_csv, time_unit='min'):
    flare = pd.read_csv(flare_csv, header=None)
    flare.reset_index(inplace=True, drop=True)
    time = flare[0]
    flux = flare[1] - 1.0
    time_factor = 1140
    # try:
    if time_unit == 'min':
        params, pcov = curve_fit(exp_decay, time[int(np.where(time==0)[0]): int(np.where(time==0)[0]) + 30], flux[int(np.where(time==0)[0]): int(np.where(time==0)[0]) + 30], maxfev=5000)
        time_factor = 1140
    elif time_unit == 'days':
        params, pcov = curve_fit(exp_decay, flux[np.where(time==0)[0]: np.where(time==0)[0] + 30]*1140, flux[np.where(time==0)[0]: np.where(time==0)[0] + 30], maxfev=5000)
        time_factor = 1
    # except:
    #     params = (0, 5.0, 0.01)
    #     time_factor = 1140
    x = np.linspace(0, 10, num = 1500)/(1140)
    y = exp_decay(x, params[0], params[1]*int(time_factor), params[2])
    amp, idx = find_nearest(y, y.max()/2)
    amp = flux.max()
    tau = x[idx]
    for index in range(len(flare[2])):
        if math.isnan(flare[2][index]) == True:
            flare[2][index] = 1e-3
    if np.abs(flare[0][0]) > 1:
        flare[0] = flare[0]/time_factor
    flare.to_csv(flare_csv, index=False, header=False)
    return amp, tau

def flare_params(data_file_dir, params_template_dir, output_dir, flares=1):
    amp, tau = allesfitter_priors(data_file_dir)
    data = pd.read_csv(data_file_dir, index_col = False, header=None)
    error = data[2]
    for index in range(len(error)):
        if math.isnan(error[index]) == True:
            error[index] = 1e-3
    flux_err = np.log(np.median(error))
    params = pd.read_csv(params_template_dir, index_col=False)
    name = (data_file_dir.rsplit('/', 1)[1])[:-4]
    params.at[1, 'value'] = 0
    params.at[2, 'value'] = tau
    params.at[3, 'value'] = amp
    params.at[5, 'value'] = flux_err
    params.at[1, 'bounds'] = 'uniform ' + str(-0.001) + ' ' + str(0.001)
    params.at[2, 'bounds'] = 'uniform ' + str(0) + ' ' + str(2.5*tau)
    params.at[3, 'bounds'] = 'uniform ' + str(amp) + ' ' + str(3.5*amp)
    params.at[5, 'bounds'] = 'uniform ' + str(-1 + flux_err) + ' ' + str(1 + flux_err)
    # if flares > 1:
    #     for number in range(2, flares+1):
    #         params.at[1+3*(number-1), '#name'] = 'flare_tpeak_' + str(number)
    #         params.at[2+3*(number-1), '#name'] = 'flare_fwhm_' + str(number)
    #         params.at[3+3*(number-1), '#name'] = 'flare_ampl_' + str(number)
            
    #     params.at[5+3*(flares-1), '#name'] = '#errors (overall scaling) per instrument'
    #     params.at[6+3*(flares-1), '#name'] = 'log_err_flux_' + str(name)
    #     params.at[6+3*(flares-1), 'label'] = '$\log{\sigma_\mathrm{' + str(name) + '}}$'
    # else:
    params.at[5, '#name'] = '#errors (overall scaling) per instrument'
    params.at[5, '#name'] = 'ln_err_flux_' + str(name)
    params.at[5, 'label'] = '$\ln{\sigma_\mathrm{' + str(name) + '}}$'
    params.to_csv(output_dir, index=False)

def return_parameters(mcmc_table_dir):
    data = pd.read_csv(mcmc_table_dir)
    t_peak = data['median'][1]
    t_peak_m = data['lower_error'][1]
    t_peak_p = data['upper_error'][1]
    fwhm = data['median'][2]
    fwhm_m = data['lower_error'][2]
    fwhm_p = data['upper_error'][2]
    amp = data['median'][3]
    amp_m = data['lower_error'][3]
    amp_p = data['upper_error'][3]

    return [t_peak, t_peak_m, t_peak_p, fwhm, fwhm_m, fwhm_p, amp, amp_m, amp_p]

def flare_energy(fwhm, ampl, Teff, R_stellar, uncertainty = False, fwhm_u = 0, fwhm_l = 0, ampl_u = 0, ampl_l = 0,
                 R_stellar_l = 0, R_stellar_u = 0):
    if uncertainty == False:
        x = np.linspace(0, 0.02, num = 2000)
        y = aflare.aflare1(x, 0.01, fwhm, ampl)
        flare_area = simpson(y, x)
        color_factor = planck_integrator(600e-6, 1000e-6, Teff)/planck_integrator(600e-6, 1000e-6, 9000)
        energy = (5.67e-8)*(9000**4)*(flare_area)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*color_factor*(1e7)*86400
        return energy
    if uncertainty == True:
        x = np.linspace(0, 0.02, num = 2000)
        y = aflare.aflare1(x, 0.01, fwhm, ampl)
        flare_area = simpson(y, x)
        color_factor = planck_integrator(600e-6, 1000e-6, Teff)/planck_integrator(600e-6, 1000e-6, 9000)
        energy = (5.67e-8)*(9000**4)*(flare_area)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*color_factor*(1e7)*86400
        samples = []
        for flare_samples in range(1000):
            fwhm_sample = fwhm
            ampl_sample = ampl
            R_stellar_sample = R_stellar
            
            check1 = np.random.random()
            check2 = np.random.random()
            check3 = np.random.random()
            if check1 < 0.5:
                err = np.abs(np.random.normal(0, scale = fwhm_l))
                fwhm_sample -= err
            elif check1 > 0.5:
                err = np.abs(np.random.normal(0, scale = fwhm_u))
                fwhm_sample += err
            if check2 < 0.5:
                err = np.abs(np.random.normal(0, scale = ampl_l))
                ampl_sample -= err
            elif check2 > 0.5:
                err = np.abs(np.random.normal(0, scale = fwhm_u))
                ampl_sample += err
            if check3 < 0.5:
                err = np.abs(np.random.normal(0, scale = R_stellar_l))
                R_stellar_sample -= err
            elif check3 > 0.5:
                err = np.abs(np.random.normal(0, scale = R_stellar_u))
                R_stellar_sample += err
            
            y = aflare.aflare1(x, 0.01, fwhm_sample, ampl_sample)
            flare_temp = np.random.normal(loc=9000, scale=500)
            color_factor = planck_integrator(600e-6, 1000e-6, Teff)/planck_integrator(600e-6, 1000e-6, flare_temp)
            energy_sample = (5.67e-8)*(9000**4)*(flare_area)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*color_factor*(1e7)*86400
            samples.append(energy_sample)
        return energy, samples[320], samples[680]

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
    shutil.copyfile(target_file, working_dir + '/' + name)
    csv_cleaner(working_dir + '/' + name)
    # Run allesfitter
    if baseline == False:
        flare_params(working_dir + '/' + name, param_template_dir, working_dir + '/params.csv')
        flare_settings(working_dir + '/' + name, settings_template_dir, working_dir + '/settings.csv', multi_process=True, cores=63, flares = 1)
    elif baseline == True:
        flare_params(working_dir + '/' + name, param_template_dir, working_dir + '/params.csv')
        flare_settings(working_dir + '/' + name, settings_template_dir, working_dir + '/settings.csv', multi_process=True, cores=63, flares=0)
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
