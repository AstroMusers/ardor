# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:29:39 2024

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.stats import vonmises
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
def Inverse_Cubic_Function(ratio, loc, e, star, planet):
    model, phase = SPI.SPI_Cubic(ratio, loc,e, star, planet, length=10)
    f = interp1d(np.linspace(0, 1, num=len(model)) ,model, kind='linear')
    return f

def Sigma_Function(sigma):
    model, phase = SPI.SPI_sigma(sigma, 100)
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

# results = SPI_Unbinned_Liklihood(phases, 100)
# print('100 Samples from a Gaussian centered at 0.8, sigma = 0.05')
# print('Unbinned Maxmimum Likeliehood: ' + str(round(results[0], 4)),
#       'Optimal Sigma: ' + str(round(results[1], 4)),
#       'Optimal Phase: ' + str(round(results[2], 4))
#       )

# print('\n100 Random Samples')
# phases = np.random.random(100)
# results = SPI_Unbinned_Liklihood(phases, 100)
# print('Unbinned Maxmimum Likeliehood: ' + str(round(results[0], 4)),
#       'Optimal Sigma: ' + str(round(results[1], 4)),
#       'Optimal Phase: ' + str(round(results[2], 4))
#       )

data = np.random.normal(loc=0.5, scale=1, size=1000)

# Define the model (Gaussian PDF)
def sigma_pdf(x, sigma):
    return Sigma_Function(sigma)(x)
def cubic_pdf(x, ratio, loc, e, star, planet):
    return Inverse_Cubic_Function(ratio, loc,e, star, planet)(x)

# Define the likelihood function
def sigma_likelihood(params, data, func):
    sigma = params
    print(-np.log(np.prod(func(data, sigma))), sigma)
    return -np.log(np.prod(func(data, sigma)))

def cubic_likelihood(params, data, star, planet):
    ratio, loc,e = params
    print(-np.log(np.prod(cubic_pdf(data, ratio, loc, e, star, planet))), ratio, loc, e)
    return -np.log(np.prod(cubic_pdf(data, ratio, loc, e, star, planet)))

star = OML.Star(1, 1, 1, radius=1, age=1e9, alfven=0.1)
planet = OML.Planet(1, 6.7, a=0.0649, e=0.01, B=10)
# model, phase= SPI.SPI_Cubic(5, star, planet) 
# plt.plot(phase, model)
# plt.show()
# Maximize the likelihood
initial_guess = [10, 0.7, 0.9]  # Initial guess for mean and standard deviation
# result = minimize(cubic_likelihood, initial_guess, args=(data, star, planet), bounds=[(1e-3,1e8)])

# Extract fitted parameters

data = pd.read_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares_New.csv")
flares = np.array(data.loc[data['Host_ID'] == 'HSPsc', ['Transit_Phase']])
# flares = (((flares - 0) * (2*np.pi)) / 1) + -np.pi
results = minimize(cubic_likelihood, initial_guess, args=(flares, star, planet), bounds=[(1e-3,1e8), (None, None), (0,0.99)])

plt.plot(np.linspace(0, 1), Inverse_Cubic_Function(results.x[0], results.x[1], results.x[2], star, planet)(np.linspace(0,1)))
print(results)
# plt.show()
# sigma_fit = result.x
# a = vonmises.fit(flares)
# plt.plot(np.linspace(-np.pi,np.pi, num=1000), vonmises.pdf(np.linspace(-np.pi,np.pi, num=1000), kappa=a[0], loc=a[1]))

# print(a)
plt.hist(flares, density=True)
plt.show()
# print("Fitted standard deviation:", sigma_fit)