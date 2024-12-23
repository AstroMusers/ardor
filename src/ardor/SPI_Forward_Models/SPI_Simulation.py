# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 20:58:28 2024

@author: Nate Whitsett
"""

import ardor.Flares.Flare as Flare
import ardor.Flares.aflare as aflare
import numpy as np
from scipy.integrate import simpson
import warnings
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
from matplotlib import pyplot as plt
warnings.filterwarnings("ignore", category=DeprecationWarning)  
warnings.filterwarnings("ignore", category=UserWarning)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }
def gaussian(x, mu, sig):
    return (
        1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )

def amp_log_normal():
    value = 1000
    while value > 0.2 or value < 0.01:
        value = np.random.lognormal(mean = 0.020744, sigma = 4.33924339)
        
    return value

def FWHM_uniform():
    return np.random.uniform(0.001388888, 0.041)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def SPI_sigma(sigma, num, duration):
    base = np.ones(num)
    phase = 0
    for index, points in enumerate(base):
        phase += 1/num
        if phase > 0.5 - duration/2 and phase < 0.5 + duration/2:
            base[index] = base[index] + gaussian(phase, 0.5, sigma)
    model = (base) / simpson(base, np.linspace(0,1,num=num))
    return model, np.linspace(0, 1, num=num)
    

def SPI_sigma_flare_injection(light_curve, SPI_parameter, SPI_duration, pl_period, sp_type = 'M', flare_type='Flaring', fast=False, theta_param = 0, phi_param = 0):
    location_list = []
    ## Approximate flare rate per 2 minute cadence of flaring M/F stars (~0.5 flares/day)
    if flare_type == 'Flaring':
        if sp_type == 'M':
            rate = 6e-4
        if sp_type == 'F':
            rate = 5e-05
        if sp_type == 'G':
            rate = 1e-5
        if sp_type == 'K':
            rate = 1e-5
    ## Poor statistics on this, but G type stars flare ~2e-5 per 2 minute cadence
    elif flare_type == 'Not Flaring':
        rate = 2.78e-8
    ## Adjust times for 20s cadence
    if fast == True:
        rate /= 6
    lc = Flare.tier0(light_curve)
    length = len(lc.time)
    phase_array = np.linspace(0, 1, length)
    model = SPI_sigma(phase_array, SPI_parameter, length, SPI_duration)
    ## Assume phase is random to begin each light curve. Periastron at 0.5
    phase = np.random.random()
    ## Iterate over the time scale of the light curve
    flares = 0
    for interval in range(length - 200):
        flare_check = np.random.random()
        flare_rate = model[find_nearest(phase_array, phase)]*rate
        if flare_rate >= flare_check:
            location = interval
            counter = 0 
            for locations in location_list:
                while location > locations - 100 and location < locations + 100 and counter < 10000:
                    location = np.random.randint(50, length-50)
                    counter += 1
            sample_baseline = lc.data[location-300:location+300]
            normalized_sample = sample_baseline/np.median(sample_baseline)
            FWHM = FWHM_uniform()
            amp = amp_log_normal()
            flare_inject = aflare.aflare1(lc.time[location-300:location+300], lc.time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            lc.data[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
            location_list.append(location)
            flares += 1
            if phase >= 0.5 - SPI_duration/2 and phase <= 0.5 + SPI_duration/2:
                print('Induced Flare')
        phase += (lc.time[interval + 1] - lc.time[interval])/pl_period 
        if phase > 1:
            phase -= 1
    return lc.data, lc.time, lc.error


def SPI_cubic_flare_injection(light_curve, star, planet, sp_type = 'M', fast=False, theta_param = 0, phi_param = 0, prior_model = True,model = None, phases = np.linspace(0, 1, num=100)):
    location_list = []
    ## Approximate flare rate per 2 minute cadence of flaring M/F stars (~0.5 flares/day)
    if sp_type == 'M':
        rate = 1.4e-4
    if sp_type == 'F':
        rate = 1e-4
    if sp_type == 'G':
        rate = 5e-5
    if sp_type == 'K':
        rate = 5e-5
    ## Adjust times for 20s cadence
    if fast == True:
        rate /= 6
    lc = Flare.tier0(light_curve)
    length = len(lc.time)
    if prior_model == False:
        model, x = OML.probability_density(star, planet, periastron_index = 1)

    ## Assume phase is random to begin each light curve. Periastron at 0.5
    random_time = np.random.random()*(planet.period)
    phase_array = ((lc.time - (random_time+ planet.period/2)) % planet.period)/planet.period
    ## Iterate over the time scale of the light curve
    flares = 0
    for interval in range(length - 200):
        phase = phase_array[interval]
        flare_check = np.random.random()
        flare_rate = model[find_nearest(phases, phase)]*rate
        if flare_rate >= flare_check:
            location = interval
            sample_baseline = lc.flux[location-300:location+300]
            normalized_sample = sample_baseline/np.median(sample_baseline)
            FWHM = FWHM_uniform()
            amp = amp_log_normal()
            flare_inject = aflare.aflare1(lc.time[location-300:location+300], lc.time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            lc.flux[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
            location_list.append(location)
            flares += 1
    return lc, random_time

def SPI_Cubic(star, planet, periastron_index, length = 10):
    probability_density = []
    bool_list = []
    for distance in planet.position:
        if distance > star.Alfven:
            bool_list.append(0)
        if distance <= star.Alfven:
            bool_list.append(1)
        probability = 1/(distance)**3
        probability_density.append(probability)
    percent = np.sum(bool_list)/len(bool_list)
    bool_list = np.array(bool_list)
    phase = planet.time/planet.period
    planet.periastron = periastron_index
    index = int(len(planet.orbit)*periastron_index/2)
    rotated_probability = probability_density[-index:] + probability_density[:-index]
    rotated_bool = np.concatenate((bool_list[int(len(bool_list)/2):],  bool_list[:int(len(bool_list)/2)]))
    uniform = np.ones(len(phase))
    rotated_probability = np.array(rotated_probability)
    rotated_probability = (rotated_probability/simpson(rotated_probability))
    probability_dist = rotated_probability*(star.Alfven/min(planet.position))*rotated_bool
    if probability_dist.sum() == 0:
        probability_dist = np.ones(len(probability_dist))
    integral = simpson(probability_dist, x= phase)
    probability_dist = probability_dist/integral
    probability_dist = (probability_dist*(planet.B/star.B)*percent + uniform)/simpson(probability_dist*(planet.B/star.B)*percent+ uniform, x= phase)
    return probability_dist, planet.phase/(2*np.pi)
