# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:20:10 2023

@author: Nate Whitsett
"""
from astropy.io import fits
from astropy.timeseries import LombScargle as LS
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.stats import linregress
import numpy as np
import pandas as pd
import os
from ardor.Utils.planck_law import planck_law as pl
from ardor.Flares import aflare
import copy
from ardor.Flares import allesfitter_priors
import shutil
import lightkurve as lk
import collections as c
from matplotlib import pyplot as plt
# import allesfitter

def linear(x, a, b, c):
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
    return np.log(a) + b*x 
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


def TESS_FITS_csv(input_file, csv_directory, csv_name=None):
    '''

    Parameters
    ----------
    input_file : string 
        Directory to the TESS light curve fits (...lc.fits) file
    csv_directory : string
        Directory for output file.
    csv_name : string, optional
        Name for output csv file

    Returns
    -------
    .cvs
        Outputs the SAP_FLUX, PDCSAP_FLUX, and time parameters from the .fits file 
        to an easily readable .csv file

    '''
    rev_file = input_file[::-1]
    index = rev_file.index('/')
    file = ((rev_file[:index])[::-1])[:-5]
    if csv_name == None:
        directory = csv_directory + '/' + file + '.csv'
    elif csv_name != None:
        directory = csv_directory + '/' + csv_name + '.csv'
    hdul = fits.open(input_file)
    time = hdul[1].data['TIME']
    data = hdul[1].data['PDCSAP_FLUX']
    error = hdul[1].data['PDCSAP_FLUX_ERR']
    
    time, data, error = delete_nans(time, data, error)

    
    grand_list = pd.DataFrame({'time': time, 'pdcsap_flux': data, 'error': error})
    grand_list = pd.DataFrame({'time': time, 'pdcsap_flux': data, 'error': error})
    return grand_list.to_csv(directory, index=False)

def TESS_data_extract(fits_lc_file, SAP_ERR=False, PDCSAP_ERR=False):
    '''

    Parameters
    ----------
    fits_lc_file : string
        Directory of the TESS light curve fits file
    SAP_ERR : bool, optional
        True will return SAP_FLUX error. The default is False.
    PDCSAP_ERR : bool, optional
        True will return PDCSAP_FLUX error. The default is False.

    Returns
    -------
    ND np.array
        Returns an Nd array of the time, PDCSAP_FLUX, SAP_FLUX, and/or the
        SAP_FLUX and PDCSAP_FLUX error. Min 3D array, max 5D array.

    '''
    hdul = fits.open(fits_lc_file)
    time = hdul[1].data['TIME']
    sap_flux = hdul[1].data['SAP_FLUX']
    sap_flux_error = hdul[1].data['SAP_FLUX_ERR']
    pdcsap_flux = hdul[1].data['PDCSAP_FLUX']
    pdcsap_flux_error = hdul[1].data['PDCSAP_FLUX_ERR']
    time, pdcsap_flux, pdcsap_flux_error = delete_nans(time, pdcsap_flux, pdcsap_flux_error)
    #pdcsap_flux_method = hdul[1].header['PDCMETHD']
    if SAP_ERR == False and PDCSAP_ERR == False:
        return time, pdcsap_flux/np.median(pdcsap_flux)
    if SAP_ERR == False and PDCSAP_ERR == True:
        return time, pdcsap_flux/np.median(pdcsap_flux), pdcsap_flux_error/np.median(pdcsap_flux)
    if SAP_ERR == True and PDCSAP_ERR == False:
        return time, sap_flux, sap_flux_error, pdcsap_flux
    if SAP_ERR == True and PDCSAP_ERR == True:
        return time, sap_flux, pdcsap_flux/np.median(pdcsap_flux), sap_flux_error, pdcsap_flux


def phase_folder(time, period, epoch):
    '''

    Parameters
    ----------
    time : numpy array
        The time array to be phase folded
    period : float
        The period of the planet, in the same units as the time array
    epoch : float
        The time of the first transit (i.e. transit epoch), same unit as time

    Returns
    -------
    phase : numpy array
        Numpy array of the same dimension as the input time array, but phase
        folded to the period of the planet. Transit centered at 0.


    '''
    phase = (time - (epoch+period/2)) % period
    phase = phase - period/2
    return phase

def flare_ID(data, sigma, fast = False):
    '''
    

    Parameters
    ----------
    data : numpy array
        Flux data for potential flares to be identified
    sigma : float
        The detection sensitivity for flares, in standard deviations. 
        For example, a sigma value of 1.0 will count any data one standard
        deviation away from the mean as a potential flare

    Returns
    -------
    flare_indices : numpy array 
        Outputs a list of potential flares as lists of indices in the provided
        data. The index begins at the triggering data point.

    '''
    mask_data = np.ma.masked_array(copy.deepcopy(data), mask=np.zeros(len(data)))
    begin = 0
    end = 1000
    shift = 0
    for index, values in enumerate(mask_data):
        if index < 1000:
            sigma2 = np.std(mask_data[0:1000])
            median = np.median(mask_data[0:1000])
        if index > (len(mask_data) - 1000):
            sigma2 = np.std(mask_data[len(mask_data)-1000:])
            median = np.median(mask_data[len(mask_data)-1000:])       
        elif shift == 1000:
            sigma2 = np.std(mask_data[begin:end])
            median = np.median(mask_data[begin:end])
            shift = 0
            begin += 1000
            end += 1000
        if mask_data[index] > (sigma*sigma2 + median):
            mask_data.mask[index] = True
        shift += 1
    plt.plot(range(len(mask_data)), mask_data)
    plt.show()
    input('E')
    plt.clf()
    plt.plot(range(len(data)), data)
    plt.show(0)
    input('E')
    plt.clf()
    flare_indices = []
    flare_length = 0
    flare_length_list = []
    flare = False
    begin = 0
    end = 1000
    shift = 0
    peak_index = 0
    sig = sigma*mask_data[begin:end].std()
    mu = np.ma.mean(mask_data[begin:end])
    if fast == False:
        for index, flux in enumerate(data):
            try:
                if shift == 1000:
                    begin += 1000
                    end += 1000
                    sig = sigma*mask_data[begin:end].std()
                    mu = np.ma.mean(mask_data[begin:end])
                    shift = 0
                if flux > (mu + sig) and flare == False:
                    flare = True
                    flare_length += 1
                    peak_index = index
                if flare == True and data[index+1] > (mu + sig/2):
                    flare_length += 1
                elif flare == True and data[index+1] < (mu + sig/2) and flare_length < 3:
                    flare = False
                    flare_length = 0
                elif flare == True and data[index+1] < (mu + sig/2) and flare_length >= 3:
                    flare = False
                    peak_correction = np.argmax(data[peak_index:peak_index+flare_length])
                    if len(flare_indices) > 0:
                        if (flare_indices[-1] + flare_length_list[-1] + 10) > peak_index:
                            continue
                        elif (flare_indices[-1] + flare_length_list[-1] + 10) <= peak_index:
                            flare_indices.append(peak_index+peak_correction)
                            flare_length_list.append(flare_length)
                    if len(flare_indices) == 0:
                        flare_indices.append(peak_index+peak_correction)
                        flare_length_list.append(flare_length)
                    flare_length = 0
                shift += 1
            except:
                print('Flare_ID_Failed')
                continue
    elif fast == True:
        for index, flux in enumerate(data):
            try:
                if shift == 1000:
                    begin += 1000
                    end += 1000
                    sig = sigma*mask_data[begin:end].std()
                    mu = np.ma.mean(mask_data[begin:end])
                    shift = 0
                if flux > (mu + sig) and flare == False:
                    flare = True
                    flare_length += 1
                    peak_index = index
                if flare == True and data[index+1] > (mu + sig/2):
                    flare_length += 1
                elif flare == True and data[index+1] < (mu + sig/2) and flare_length < 8:
                    flare = False
                    flare_length = 0
                elif flare == True and data[index+1] < (mu + sig/2) and flare_length >= 8:
                    flare = False
                    peak_correction = np.argmax(data[peak_index:peak_index+flare_length])
                    if len(flare_indices) > 0:
                        if (flare_indices[-1] + flare_length_list[-1] + 10) > peak_index:
                            print(flare_indices[-1], flare_length_list[-1], peak_index + peak_correction)
                            continue
                        elif (flare_indices[-1] + flare_length_list[-1] + 10) <= peak_index:
                            flare_indices.append(peak_index+peak_correction)
                            flare_length_list.append(flare_length)
                    if len(flare_indices) == 0:
                        flare_indices.append(peak_index+peak_correction)
                        flare_length_list.append(flare_length)
                    flare_length = 0
                shift += 1
            except:
                print('Flare_ID_Failed')
                continue
    
    return flare_indices, flare_length_list

def delete_nans(time, data, error):
    '''
    
    Parameters
    ----------
    time : numpy array
        The time array to be cleared of NANs. Must be the same dimensionality
        and 1-to-1 with the data array.
    data : numpy array
        The data array to be cleared of NANs. Must be the same dimensionality
        and 1-to-1 with the time array

    Returns
    -------
    time1 : numpy array
        Returns the original time array, but with any NANs in the data or
        time array removed for both arrays, at the same indices, such that
        both arrays are still 1-to-1.
    data1 : numpy array
        Returns the original data array, but with any NANs in the data or
        time array removed for both arrays, at the same indices, such that
        both arrays are still 1-to-1

    '''
    time = time.tolist()
    data = data.tolist()
    error = error.tolist()
    nan_set = set()
    count_data = 0
    count_time = 0
    count_error = 0
    for indices in data:
        if np.isnan(indices) == True:
            nan_set.add(count_data)
        count_data += 1
    for indices in time:
        if np.isnan(indices) == True:
            nan_set.add(count_time)
        count_time += 1
    for indices in error:
        if np.isnan(indices) == True:
            nan_set.add(count_error)
        count_time += 1
    time1 = [i for j, i in enumerate(time) if j not in nan_set]
    data1 = [i for j, i in enumerate(data) if j not in nan_set]
    error1 = [i for j, i in enumerate(error) if j not in nan_set]
    time1 = np.array(time1)
    data1 = np.array(data1)
    error1 = np.array(error1)
    return time1, data1, error1

def SMA_detrend(time, data, error, LS_Iterations=3, time_scale=100, model=False):
    '''
    
    This applies a windowed, Single Moving Average to detrend potentially
    periodic data.
    Parameters
    ----------
    data : numpy array
        The flux data to be detrended
    time_scale : int
        The time-scale to apply the detrending.

    Returns
    -------
    numpy array
        Returns the detrended data array.
    '''
    
    time, data, error = delete_nans(time, data, error)
    count = 0
    ls = LS(time,data,error)
    frequency,power = ls.autopower(minimum_frequency=0.1, maximum_frequency=100)
    cutoff = ls.false_alarm_probability(power.max()) 
    if cutoff < 0.3 and count < LS_Iterations:
        ls = LS(time,data,error, nterms = 3)
        frequency, power = ls.autopower(minimum_frequency=0.1, maximum_frequency=100)
        best_frequency = frequency[np.argmax(power)]
        theta = ls.model_parameters(best_frequency)
        offset = ls.offset()
        design_matrix = ls.design_matrix(best_frequency, time)
        data = data - (design_matrix.dot(theta))
        LS_model = (offset + design_matrix.dot(theta))
        count += 1
    mask_data = np.ma.masked_array(copy.deepcopy(data), mask=np.zeros(len(data)))
    begin = 0
    end = 200
    shift = 0
    for index, values in enumerate(mask_data):
        if index < 200:
            sigma2 = np.std(mask_data[0:200])
            median = np.median(mask_data[0:200])
        if index > (len(mask_data) - 200):
            sigma2 = np.std(mask_data[len(mask_data)-200:])
            median = np.median(mask_data[len(mask_data)-200:])       
        elif shift == 200:
            sigma2 = np.std(mask_data[begin:end])
            median = np.median(mask_data[begin:end])
            shift = 0
            begin += 200
            end += 200
        if np.abs(mask_data[index]) > (3*sigma2 + median):
            mask_data.mask[index] = True
        shift += 1
    mov_average = []
    j = 0
    i = 0
    for a in range(time_scale - 1):
        window = mask_data.data[a : 1 + a + j]
        mov_average.append(sum(window)/(j+1))
        j += 1
    while i < len(data) - time_scale + 1:
        window = mask_data.data[i : i + time_scale]
        window_average = np.ma.round(np.ma.sum(window) / time_scale, 2)
        mov_average.append(window_average)
        i += 1
    SMA = mask_data.data - np.array(mov_average)
    if model == False:
        return mask_data.data
    elif model == True:
        return SMA + 1, mov_average
        
def lk_detrend(data, time, scale=401, return_trend = False):
    if return_trend == False:
        lc = lk.LightCurve(flux = data, time = time).flatten(scale, sigma=3, return_trend=return_trend)
        return lc.flux
    elif return_trend == True:
        lc, trend = lk.LightCurve(flux = data, time = time).flatten(scale, sigma=3, return_trend=return_trend)
        return lc.flux, trend

def flare_phase_folded_ID(phase, flare_array, period, epoch):
    new_ID_list = []
    for indices in flare_array:
        new_ID_list.append(((phase[indices] - (epoch+period/2)) % period)-period/2)
    return np.array(new_ID_list)

def bolo_flare_energy(parameters, R_stellar, planck_ratio, t_unit='days', function=exp_decay):
    a, b, c = parameters
    if t_unit == 'days':
        multiplier = 86400
        length_cap = 0.08333
    if t_unit == 'minutes':
        multiplier = 60
        length_cap = 120
    if function == exp_decay:
        integral, err = quad(function, 0, length_cap, args=(a, b, (c-1)))
    elif function == aflare.aflare1:
        integral, err =quad(function, -length_cap, length_cap, args=(a, b, c))
    energy = (5.67e-8)*(9000**4)*(integral)*np.pi*(R_stellar*6.957e8*R_stellar*6.957e8)*planck_ratio*(1e7)*multiplier
    return energy

def tier0(TESS_fits_file, scale = 401):
    '''
    Tier 0 of ardor. This function accepts a TESS '...lc.fits' file, and returns
    a named tuple which contains a NAN free, detrended and normalized 
    light curve. Additionally returns the observation time, as well as a boolean
    denoting if it is 2 minute or 20 second cadence data, as well as the 
    derived trend in the detrending process.

    Parameters
    ----------
    TESS_fits_file : string
        The TESS light curve you wish to detrend and clean up.
    Returns
    -------
    LightCurve : named tuple
        A named tuple which has keys:
            - time: time, BJD. (array)
            - flux: normalized flux. (array)
            - detrended_flux: Detrended, normalized flux. (array)
            - error: Normalized error. (array)
            - fast_bool: Boolean denoting if the TESS data is 2 minute or 20
            second cadence. fast == True means 20 second cadence. (bool)
            - obs_time: The total observation time reflected by the data, in
            minutes. This only counts data present in the file, excludes gaps
            or NAN values. (float)
            - trend: The trend removed in the detrending step. (array)
    '''
    if TESS_fits_file.endswith('a_fast-lc.fits') == True:
        fast = True
        cadence = (1/3)
    elif TESS_fits_file.endswith('a_fast-lc.fits') == False:  
        fast = False
        cadence = 2
    b, pdcsap_flux, pdcsap_error = TESS_data_extract(TESS_fits_file, PDCSAP_ERR=True)
    time, flux, pdcsap_error = delete_nans(b, pdcsap_flux, pdcsap_error)
    detrend_flux, trend = lk_detrend(flux, time, scale=scale, return_trend= True)
    observation_time = cadence*len(time)
    LightCurve = c.namedtuple('LightCurve', ['time', 'flux', 'detrended_flux', 'error', 'fast_bool', 'obs_time', 'trend'])
    lc = LightCurve(time, flux, detrend_flux, pdcsap_error, fast, observation_time, trend.flux)
    return lc


def tier1(detrend_flux, sigma, fast=False):
    '''
    

    Parameters
    ----------
    detrend_flux : numpy array
        The median centered, detrended pdcsap flux data of the TESS light curve file, in units of electrons/second. 
        Intended to be from the second output of the tier 0 function
    sigma : float
        The sensitivity cutoff in which you wish to search for flares. Typically, this is 3 sigma, though some pipelines
        go as low as 2.5. Will NOT work well below 2 sigma.
    Returns
    -------
    named tuple, the the following keys:
        flares : numpy array
            Returns an array of the **indices** of the flares in the time array passed. To get the BJD time of a given flare, 
            pass 'time[flares[0]]', etc. Used in tier 2.
        lengths : numpy array
            A 1:1 numpy array with flares, giving the approximate duration of the flare in units of **indices** of the time 
            axis. This will either be two minutes/index or 20 seconds an index, depending on what type of light curve is used.
            Used in tier 2.

    '''
    flares, lengths = flare_ID(detrend_flux, sigma, fast=fast)
    Flare = c.namedtuple('Flares', ['index', 'length'])
    flare = Flare(flares, lengths)
    return flare

def tier2(time, flux, pdcsap_error, flares, lengths, chi_square_cutoff = 1,
          output_dir = 'Output.csv', host_name = 'My_Host', T = 4000, 
          host_radius = 1, csv = True, planet_period = 5, planet_epoch = 1000, 
          Sim = False, injection = False, const = 0, obs_time=0, extract_window = 50):
    '''
    
    Parameters
    ----------
    time : numpy array
        The time data of the TESS light curve file. Intended to be from the first output of the tier 0 function.
        Units of BJD - 255700 (Days)
    detrend_flux : numpy array
        The median centered, detrended pdcsap flux data of the TESS light curve file, in units of electrons/second. 
        Intended to be from the second output of the tier 0 function
    pdcsap_error : numpy array
        The uncertainty in the detrended flux. Intended to be from the third output of the tier 0 function.
    flares : numpy array
        Returns an array of the **indices** of the flares in the time array passed. To get the BJD time of a given flare, 
        pass 'time[flares[0]]', etc. Intended as the first output from tier 2.
    lengths : numpy array
        A 1:1 numpy array with flares, giving the approximate duration of the flare in units of **indices** of the time 
        axis. This will either be two minutes/index or 20 seconds an index, depending on what type of light curve is used.
        Intended as the first output from tier 2.
    output_dir : string
        The directory in which you wish to save the flare csv files to.
    Teff : float
        The effective stellar temperature of the host star. If not passed, 4000 K will be assumed. Used to estimate 
        flare energies.
    host_radiu : float
        The radius of the host star. If not passed, a value of 1 solar radii is assumed. Used to estimate the flare energies.

    Returns
    -------
    Snippets of the flares found from 

    '''
    param = 1
    TOI_ID_list = []
    flare_number = []
    peak_time = []
    amplitude = []
    time_scale = []
    Teff = []
    radius = []
    flare_amplitude = []
    flare_time_scale = []
    accepted_flare_index = []
    accepted_flare_number = []
    param_list = []
    phase_list = []
    event_list = []
    chi_square_list = []
    obs_time_list = []
    flare_count = 0
    total_flares = 0
    flux = np.array(flux)
    if csv == True:
        os.makedirs(output_dir + '/' + str(host_name), exist_ok=True)
    for index, flare_events in enumerate(flares):
        if flare_events >= extract_window and len(flux) - flare_events > extract_window:
            new_time = time[flare_events-extract_window:flare_events+extract_window]
            new_data = flux[flare_events-extract_window:flare_events+extract_window]
            new_error = pdcsap_error[flare_events-extract_window:flare_events+extract_window]
            check = 1
        elif flare_events < extract_window:
            new_time = time[:flare_events+extract_window]
            new_data = flux[:flare_events+extract_window]
            new_error = pdcsap_error[0+flare_events:flare_events+extract_window]
            check = 2
        elif len(flux) - flare_events < extract_window:
            new_time = time[flare_events:]
            new_data = flux[flare_events:]
            new_error = pdcsap_error[flare_events:]
            check = 3
        norm_time = time[flare_events]
        new_time = (new_time - norm_time)*24*60
        if check == 1:
            events = extract_window
        elif check == 2:
            events = flare_events
        elif check == 3:
            events = len(flux) - flare_events
        r_sq = 0
        try:
            if lengths[index] >= 15:
                popt, pcov = curve_fit(exp_decay, new_time[events:events+15],new_data[events:events+15], maxfev=5000, sigma = new_error[events:events+15], absolute_sigma=True)
                squares = ((new_data[events:events+15] - exp_decay(new_time[events:events+15], *popt))/((new_error[events:events+15])))**2
                chi_squared = np.sum(squares)/13
                x = linregress(new_time[events:events+10], (np.log(new_data[events:events+10] - new_data[events+10])))
                r_sq = x.rvalue**2
            elif lengths[index] > 5 and lengths[index] < 15:
                popt, pcov = curve_fit(exp_decay, new_time[events:events+10], new_data[events:events+10], maxfev=5000, sigma = new_error[events:events+10], absolute_sigma=True)
                squares = ((new_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))/(new_error[events:events+10]))**2
                chi_squared = np.sum(squares)/8
                x = linregress(new_time[events:events+10], (np.log(new_data[events:events+10] - new_data[events+10])))
                r_sq = x.rvalue**2
            elif lengths[index] == 5:
                popt, pcov = curve_fit(exp_decay, new_time[events:events+5], new_data[events:events+5], maxfev=5000, sigma = new_error[events:events+5], absolute_sigma=True)
                squares = ((new_data[events:events+5] - exp_decay(new_time[events:events+5], *popt))/(new_error[events:events+5]))**2
                chi_squared = np.sum(squares)/(3)
            elif lengths[index] == 4:
                popt, pcov = curve_fit(exp_decay, new_time[events:events+4], new_data[events:events+4], maxfev=5000, sigma = new_error[events:events+4], absolute_sigma=True)
                squares = ((new_data[events:events+4] - exp_decay(new_time[events:events+4], *popt))/(new_error[events:events+4]))**2
                chi_squared = np.sum(squares)/(2)
            elif lengths[index] == 3:
                popt, pcov = curve_fit(exp_decay, new_time[events:events+3], new_data[events:events+3], maxfev=5000, sigma = new_error[events:events+3], absolute_sigma=True)
                squares = ((new_data[events:events+3] - exp_decay(new_time[events:events+3], *popt))/(new_error[events:events+3]))**2
                chi_squared = np.sum(squares)/(1)
        except:
            print('Flare ID Error')
            chi_squared = 100
        if (chi_squared < chi_square_cutoff and popt[0] > 0 and popt[1] > 0) or r_sq > 0.8:
            event_list.append(flare_events)
            flare_count += 1
            if Sim == True:
                param_list.append(param)
            time_scale.append(popt[1])
            flare_time_scale.append(popt[1])
            amplitude.append(popt[0])
            flare_amplitude.append(popt[0])
            peak_time.append(norm_time)
            TOI_ID_list.append(host_name)
            flare_number.append(flare_count)
            chi_square_list.append(chi_squared)
            Teff.append(T)
            radius.append(host_radius)
            obs_time_list.append(obs_time)
            phase_list.append(((time[flare_events] - (planet_epoch+planet_period/2)) % planet_period)/planet_period)
            try:
                X = np.column_stack((new_time, new_data, new_error))
            except:
                X = np.column_stack(([0],[0],[0]))
            accepted_flare_index.append(flares[index])
            accepted_flare_number.append(flare_count)
            if csv == True:
                np.savetxt(output_dir + '/' + str(host_name) + '/Flare_' + str(flare_count + const) + '.csv', X, delimiter=',')
            total_flares += 1
        
        else:
            continue
    if Sim == False and injection == False:
        ZZ = np.column_stack((TOI_ID_list, np.array(flare_number) + const, peak_time, amplitude, time_scale, Teff, radius, phase_list,obs_time_list, chi_square_list))
        return ZZ, flare_count
    if Sim == True:
        ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, Teff, radius, phase_list,obs_time_list, chi_square_list))
        return ZZ
    if csv == True:
        ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, Teff, radius,phase_list, obs_time_list, chi_square_list))
        with open(output_dir + '/' + str(host_name) + '/All_Flare_Parameters.csv', "a") as f:
            np.savetxt(f, ZZ, delimiter=",", fmt='%s')
            f.close()
        return ZZ
    if injection == True and Sim == False:
        return event_list
                    
def tier3(tier_2_output_dir, tier_3_working_dir, tier_3_output_dir, settings_template_dir, params_template_dir, host_name = 'My_Host', T = 4000, host_radius = 1, MCMC_CPUS = 1):
    #Check each folder
    flare_files = os.listdir(tier_2_output_dir)
    for csvs in flare_files:
        parameter_list_list = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
        # try:
        #Make output directory to store plots and flare data for each target
        #Omit the statistics files
        if csvs == 'All_Flare_Parameters.csv' or csvs == 'Host_Statistics.txt' or csvs == 'Flare_Phase.csv':
            continue
        # Copy relevant file to allesfitter folder
        shutil.copyfile(tier_2_output_dir + '/' + csvs, tier_3_working_dir+ '/' +  csvs)
        allesfitter_priors.csv_cleaner(tier_3_working_dir + '/' + csvs)
        # Run allesfitter
        allesfitter_priors.flare_params(tier_3_working_dir + '/' +  csvs, params_template_dir, tier_3_working_dir + '/params.csv')
        allesfitter_priors.flare_settings(tier_3_working_dir + '/' +  csvs, settings_template_dir, tier_3_working_dir + '/settings.csv', multi_process=True, cores=MCMC_CPUS)
        allesfitter.mcmc_fit(tier_3_working_dir)
        allesfitter.mcmc_output(tier_3_working_dir)

        #Extract relevant parameters and append them to relevant lists
        # Full Data csv
        parameter_list = allesfitter_priors.return_parameters(tier_3_working_dir + '/results/mcmc_table.csv')
        if len(parameter_list) == 0:
            continue
        energy = allesfitter_priors.flare_energy(parameter_list[3], parameter_list[6], T,  host_radius)
        parameter_list_list[0].append(host_name)
        parameter_list_list[1].append(csvs[:-4])
        for index in range(9):
            parameter_list_list[index+2].append(parameter_list[index])
        parameter_list_list[11].append(energy)
        parameter_list_list[12].append(T)
        parameter_list_list[13].append(host_radius)


        #Per target csv
        #Remove current csv
        os.remove(tier_3_working_dir + '/' +  csvs)
        #Copy relevant graphs and log files
        shutil.copyfile(tier_3_working_dir + '/results/mcmc_corner.pdf', tier_3_output_dir + '/mcmc_corner_' + csvs[:-4] + '.pdf')
        shutil.copyfile(tier_3_working_dir + '/results/mcmc_fit_b.pdf', tier_3_output_dir + '/mcmc_fit_' + csvs[:-4] + '.pdf')
        result_dir = os.listdir(tier_3_working_dir + '/results')
        for dirs in result_dir:
            if dirs[-3] == 'l':
                shutil.copyfile(tier_3_working_dir + '/results/' + dirs, tier_3_output_dir + '/mcmc_' + csvs[:-4] + '.log')
        #Delete results folder
        for files in os.listdir(tier_3_working_dir + '/results'):
            os.remove(tier_3_working_dir + '/results/' + files)
        # except:
        #     print(1)
        #     parameter_list_list[0].append(host_name)
        #     parameter_list_list[1].append(csvs[:-4])
        #     for index in range(9):
        #         parameter_list_list[index+2].append(np.nan)
        #     parameter_list_list[11].append(np.nan)
        #     parameter_list_list[12].append(T)
        #     parameter_list_list[13].append(host_radius)
        ZZ = np.column_stack((parameter_list_list[0], parameter_list_list[1], parameter_list_list[2], parameter_list_list[3], parameter_list_list[4], parameter_list_list[5], parameter_list_list[6], parameter_list_list[7], parameter_list_list[8], parameter_list_list[9], parameter_list_list[10], parameter_list_list[11], parameter_list_list[12], parameter_list_list[13]))
        with open(tier_3_output_dir + '/All_TOI_MCMC_Flares.csv', "a") as f:
            np.savetxt(f, ZZ, delimiter=",", fmt='%s')
            f.close()


# font = {'family': 'serif',
#         'color':  'black',
#         'weight': 'normal',
#         'size': 14,
#         }
# plt.plot(time, data)



# ### Pipeline Graphics for Publication
# L = ax1.legend(prop={'size': 9.3}, ncol=2)
# plt.setp(L.texts, family='Serif')
# x = np.linspace(0, 6, num = 50)
# popt, pcov = curve_fit(exp_decay, x, detrend[8856:8906])
# x = np.linspace(0,0.02, num = 50)
# ax2 = fig.add_subplot(gs[1, 0])
# ax2.plot(x + time[8856], exp_decay(x, popt[0], popt[1], popt[2]), zorder=10, linewidth=2, linestyle='--', c='blue', label = 'Tier 2 Curve Fit')
# ax2.plot(time[8830:8880], detrend[8830:8880],c='darkorange', linewidth=1.25)
# ax2.set_xticks([])
# ax2.set_yticks([])
# ax2.set_ylabel('Flux', fontdict = font)
# ax2.set_xlabel('Time (BJD)', fontdict = font)
# L1 = ax2.legend(prop={'size':10})
# plt.setp(L1.texts, family='Serif')



# ax3 = fig.add_subplot(gs[1, 1])
# x = np.linspace(0, 6, num = 55)
# popt, pcov = curve_fit(exp_decay, x, detrend[9155:9210])
# x = np.linspace(0,0.03, num = 50)
# ax3.plot(x + time[9154], exp_decay(x, popt[0], popt[1], popt[2]), zorder=10, linestyle='--', linewidth=2, c='blue', label = 'Tier 2 Curve Fit')
# ax3.plot(time[9120:9190], detrend[9120:9190],c='darkorange', linewidth=1.25, label = 'Detrended Flux')
# ax3.set_xticks([])
# ax3.set_yticks([])
# ax3.set_ylabel('Flux', fontdict = font)
# ax3.set_xlabel('Time (BJD)', fontdict = font)
# plt.setp(L1.texts, family='Serif')
# plt.savefig('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Ardor.png', dpi=600,bbox_inches='tight')
# plt.show()


# TESS_FITS_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/TOI_Hosts/1410/tess2019253231442-s0016-0000000199444169-0152-s_lc.fits', 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop')