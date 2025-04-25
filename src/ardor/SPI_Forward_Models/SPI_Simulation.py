# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 20:58:28 2024

@author: Nate Whitsett
"""

import ardor.Flares.Flare as Flare
import ardor.Flares.aflare as aflare
import numpy as np
from scipy.integrate import simpson
from scipy.stats import uniform
import warnings
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
from scipy.stats import vonmises
import ardor.Utils.Utils as UT
from scipy.integrate import simpson
warnings.filterwarnings("ignore", category=DeprecationWarning)  
warnings.filterwarnings("ignore", category=UserWarning)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }
#### FORWARD MODEL PDFS ####
def SPI_kappa(kappa,loc, num):
    phase = np.linspace(0,1,num=num)              
    base = vonmises.pdf(UT.range_shift(phase, 0, 1, -np.pi, np.pi), kappa = kappa, loc = UT.range_shift(loc, 0, 1, -np.pi, np.pi))*2*np.pi
    return base, phase
def SPI_Cubic(ratio, loc, e, a, period, length = 200):
    probability_density = []
    # bool_list = []
    time, position, rot_time = OML.orbit_pos_v_time(period, e, a, orbit_length = length,  arg_periastron=loc)
    for distance in position:
    #     if distance > star.Alfven:
    #         bool_list.append(0)
    #     if distance <= star.Alfven:
    #         bool_list.append(1)
        probability = ((distance)/(a*(1-e)))**-3
        probability_density.append(probability)
    phase = (time/period)
    rot_phase = rot_time/period
    uni = np.ones(len(probability_density))
    normalized_dist = probability_density/simpson(probability_density, x = phase, dx = 0.01)
    probability_density = np.array((np.array(normalized_dist)*ratio + np.array(uni))/simpson(np.array(normalized_dist)*ratio + np.array(uni),x=phase, dx=0.01))
    print(simpson(probability_density, x=phase, dx=0.01))
    return rot_phase, probability_density


##### FLARE INJECTION IMPLEMENTATIONS #######
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


def SPI_kappa_flare_injection(light_curve, kappa, loc, pl_period, sp_type = 'M', flare_type='Flaring', fast=False, theta_param = 0, phi_param = 0):
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
    ## Poor statistics on this, but G type stars flare ~2e-5 per 2 minute cadence
    elif flare_type == 'Not Flaring':
        rate = 2.78e-8
    ## Adjust times for 20s cadence
    if fast == True:
        rate /= 6
    lc = Flare.TESS_data_extract(light_curve)
    length = len(lc.time)
    model, a = SPI_kappa(kappa, loc, length)
    ## Assume phase is random to begin each light curve. Periastron at 0.5
    random_time = lc.time[int(len(lc.time)/2)]
    phase_array = ((lc.time - (random_time+ pl_period/2)) % pl_period)/pl_period
    ## Iterate over the time scale of the light curve
    flares = 0
    for interval in range(length - 200):
        phase = phase_array[interval]
        flare_check = np.random.random()
        flare_rate = model[find_nearest(a, phase)]*rate
        if flare_rate >= flare_check:
            location = interval
            counter = 0 
            for locations in location_list:
                while location > locations - 100 and location < locations + 100 and counter < 10000:
                    location = np.random.randint(50, length-50)
                    counter += 1
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
