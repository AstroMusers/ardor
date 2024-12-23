# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:29:39 2024

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import numpy as np
def Inverse_Cubic_Function(star, planet, periastron_index):
    model, phase = SPI.SPI_Cubic(star, planet, periastron_index, length = 1000)
    f = interp1d(phase ,model, kind='linear')
    return f

def Sigma_Function(sigma, length, duration):
    model, phase = SPI.SPI_sigma(sigma, length, duration)
    f = interp1d(phase ,model, kind='linear')
    return f
def find_2d_index_of_max(arr):
    max_value = np.max(arr)  # Find the maximum value in the array
    max_index = np.where(arr == max_value)  # Get the indices of the max value
    return max_index[0][0], max_index[1][0]  # Return the first occurrence (row, col)

def SPI_Unbinned_Liklihood(flare_phases, iterations, null_model = 'Uniform', model = 'Sigma'):
    log_m_null_list_list = []
    log_m_model_list_list = []
    sigma_space = np.logspace(0, -2, num=10)
    phase_space = np.linspace(0, 1, num = 100)
    for sigmas in sigma_space:
        log_m_null_list = []
        log_m_model_list = []
        for phase_diff in phase_space:
            model = Sigma_Function(sigmas, 100, 1)
            log_m_null = 0
            log_m_model = 0
            for phases in flare_phases:
                log_m_null += np.log(1/len(flare_phases))
                log_m_model += np.log(model(phases)/len(flare_phases))
            log_m_null_list.append(log_m_null)
            log_m_model_list.append(log_m_model)
            flare_phases += 1/len(phase_space)
            for index, phase in enumerate(flare_phases):
                if flare_phases[index] > 1:
                    flare_phases[index] -= 1
                elif flare_phases[index] < 0:
                    flare_phases[index] += 1
        log_m_model_list_list.append(log_m_model_list)
        log_m_null_list_list.append(log_m_null_list)
    log_m_model_list_list = np.array(log_m_model_list_list)
    log_m_null_list_list = np.array(log_m_null_list_list)
    sigma_index, phase_index = np.unravel_index(log_m_model_list_list.argmax(), log_m_model_list_list.shape)
    optimal_simga = sigma_space[sigma_index]
    optimal_phase = 0.5 - phase_space[phase_index]
    if optimal_phase < 0:
        optimal_phase += 1
    return -2*(np.max(log_m_null_list_list) - np.max(log_m_model_list_list)), optimal_simga, optimal_phase

phases = np.random.normal(0.8, scale = 0.05, size = 1000)
results = SPI_Unbinned_Liklihood(phases, 100)
print('100 Samples from a Gaussian centered at 0.8, sigma = 0.05')
print('Unbinned Maxmimum Likeliehood: ' + str(round(results[0], 4)),
      'Optimal Sigma: ' + str(round(results[1], 4)),
      'Optimal Phase: ' + str(round(results[2], 4))
      )

print('\n100 Random Samples')
phases = np.random.random(100)
results = SPI_Unbinned_Liklihood(phases, 100)
print('Unbinned Maxmimum Likeliehood: ' + str(round(results[0], 4)),
      'Optimal Sigma: ' + str(round(results[1], 4)),
      'Optimal Phase: ' + str(round(results[2], 4))
      )