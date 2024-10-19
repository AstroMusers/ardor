# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 20:58:28 2024

@author: Nate Whitsett
"""

import Flare
import aflare as af
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os
import time as t
import csv 
from astroquery.mast import Catalogs
from scipy.stats import gaussian_kde
from scipy.integrate import simpson
import warnings
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

def SPI_model(phase, sigma, length, duration):
    base = np.ones(length)
    model = gaussian(phase, 0.5, sigma)
    phase = 0
    for index, points in enumerate(base):
        phase += 1/length
        if phase > 0.5 - duration/2 and phase < 0.5 + duration/2:
            base[index] = base[index] + model[index]
    model = (base) / simpson(base, np.linspace(0,1,num=length))
    return model
    

def SPI_flare_injection(light_curve, SPI_parameter, SPI_duration, pl_period, sp_type = 'M', flare_type='Flaring', fast=False, theta_param = 0, phi_param = 0):
    location_list = []
    ## Approximate flare rate per 2 minute cadence of flaring M/F stars (~0.5 flares/day)
    if flare_type == 'Flaring':
        if sp_type == 'M':
            rate = 2.8e-4
        if sp_type == 'F':
            rate = 6e-05
        if sp_type == 'G':
            rate = 1.18e-4
        if sp_type == 'K':
            rate = 1.19e-4
    ## Poor statistics on this, but G type stars flare ~2e-5 per 2 minute cadence
    elif flare_type == 'Not Flaring':
        rate = 2.78e-8
    ## Adjust times for 20s cadence
    if fast == True:
        rate /= 6
    data, time, error = Flare.TESS_data_extract(light_curve, PDCSAP_ERR=True)
    data, time, error = Flare.delete_nans(time, data, error)
    phase_array = np.linspace(0, 1, len(data))
    model = SPI_model(phase_array, SPI_parameter, len(data), SPI_duration)
    ## Assume phase is random to begin each light curve. Periastron at 0.5
    phase = np.random.random()
    phaseb = phase
    ## Iterate over the time scale of the light curve
    flares = 0
    for interval in range(len(time) - 200):
        flare_check = np.random.random()
        flare_rate = model[find_nearest(phase_array, phase)]*rate
        if flare_rate >= flare_check:
            location = interval
            counter = 0 
            for locations in location_list:
                while location > locations - 100 and location < locations + 100 and counter < 10000:
                    location = np.random.randint(50, len(data)-50)
                    counter += 1
            sample_baseline = data[location-300:location+300]
            normalized_sample = sample_baseline/np.median(sample_baseline)
            FWHM = FWHM_uniform()
            amp = amp_log_normal()
            flare_inject = af.aflare1(time[location-300:location+300], time[location], FWHM, amp)
            normalized_sample_inject = normalized_sample + flare_inject
            data[location-300:location+300] = np.median(sample_baseline)*normalized_sample_inject
            location_list.append(location)
            flares += 1
            if phase >= 0.5 - SPI_duration/2 and phase <= 0.5 + SPI_duration/2:
                print('Induced Flare')
        phase += (time[interval + 1] - time[interval])/pl_period 
        if phase > 1:
            phase -= 1
    return data, time, error, phaseb, flares
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }
SPI_Params = [0.5, 0.1, 0.05]
SPI_Params2 =  [1, 5, 1000]
# for params in SPI_Params:
#     model = SPI_model(np.linspace(0, 1, num=1000), params, 1000, 0.4)
#     plt.plot(np.linspace(0,1,num=1000), model, linestyle= '-', alpha = 1, label = '$\sigma = $' + str(params))
# plt.plot(np.linspace(0,1, num=1000), np.ones(1000), label='Uniform', linestyle=':')
# L = plt.legend(prop={'size': 10})
# plt.text(0.16, 2.5, '$\phi_{\mathrm{Lower}}$', size =13)
# plt.text(0.725, 2.5, '$\phi_{\mathrm{Upper}}$', size =13)
# plt.axvline(0.3, color='black', linestyle='--', alpha=0.75)
# plt.axvline(0.7, color='black', linestyle='--', alpha=0.75)
# plt.setp(L.texts, family='Serif')
# plt.xlim(0,1)

# plt.xlabel('Phase', fontdict = font)
# plt.ylabel('Probability Density', fontdict = font)
# plt.savefig('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Sim_Prob_Density.png', dpi=400, bbox_inches='tight')
fileM = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/M Type TESS Data/Gaia DR3 1118597710821365248/tess2022164095748-s0053-0000000219698348-0226-s_lc.fits'
file = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/G Type TESS Data/G 5.51838e+18/tess2019032160000-s0008-0000000143196764-0136-s_lc.fits'
# data, time, error, phase = flare_injection(file, 0.5, 0.1, 7,0.1, sp_type = 'M')
# plt.plot(time, data)
# x = np.linspace(time[0], time[len(time)-1], num=10000)
# plt.axvspan(1521.04, 1521.74, color='red', alpha=0.4, label='Sub-Alfvenic')
# plt.axvspan(1528.04, 1528.74, color='red', alpha=0.4)
# plt.axvspan(1535.04, 1535.74, color = 'red', alpha=0.4)
# plt.title('Simulated Phase Curve (Quiescent G Type)')
# plt.xlabel('Time (BJD - 2457000)')
# plt.ylabel('PDCSAP Flux (electrons $s^{-1}$)')
# plt.legend()
# plt.plot(time, data, linewidth = 1)
# # plt.savefig('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Publication Documents/Simulated_LC.png', dpi=1000, bbox_inches = 'tight')
# plt.show()

for params in SPI_Params2:
    for simulations in range(20):
        print('Star: ' + str(simulations))
        for G_Stars in range(30):
            print('LC ' + str(G_Stars))
            data, time, error, phase, flare_count = SPI_flare_injection(fileM, params, 0.1, 7, sp_type = 'M')
            print(len(time))
            # detrend_data = Flare.SMA_detrend(time, data, error, time_scale = 60)
            flares, lengths = Flare.tier1(data, sigma=3, fast=False)
            if phase < 0.5:
                epoch = 0.5 + phase
            if phase > 0.5:
                epoch = phase - 0.5
            output_dir = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_type_SPI_Sim'
            ZZ = Flare.tier2(time, data, error, flares, lengths, chi_square_cutoff=3, host_name = 'M_Sim_' + str(simulations),T = 5500, host_radius = 1, csv=False, planet_period = 7, planet_epoch=time[0]-epoch*7,Sim=True, param = params)
            if flare_count != 0:
                print(str(len(ZZ)/flare_count*100) + '% Flares Recovered')
            try:
                with open(output_dir + '/M_Flaring_SPI_' +  str(params) + '.csv', "a") as f:
                    np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                    f.close()
            except:
                continue

# file = 'C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/Literature_Catalogs/Pietras_2020.csv'
# file2 = 'C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/Literature_Catalogs/TIC_Spec_Flares.csv'
# data = pd.read_csv(file)
# data2 = pd.read_csv(file2)
# flares = list(data2['Sector'])
# time = list(data2['Obs Time'])
# spectral = list(data2['Spectral Type'])
# OBA = []
# G = []
# F = []
# K = []
# M = []
# for index, flares in enumerate(flares):
#     if spectral[index] == 'A/B/O':
#         OBA.append(flares/time[index])
#     if spectral[index] == 'F':
#         F.append(flares/time[index])
#     if spectral[index] == 'G':
#         G.append(flares/time[index])
#     if spectral[index] == 'K':
#         K.append(flares/time[index])
#     if spectral[index] == 'M':
#         M.append(flares/time[index])
# OBA = np.array(OBA)
# G = np.array(G)
# F = np.array(F)
# K = np.array(K)
# M = np.array(M)

# fig = plt.figure(figsize=(5,5))
# gs = fig.add_gridspec(1)
# # ax1 = fig.add_subplot(gs[0])
# ax2 = fig.add_subplot(gs[0])
# # count, bins, ignored = ax1.hist(data['Amp'], range =(0, 0.5), bins = 500, density=True, align='mid', color='blue', label = 'Pietras 2022')
# # mu = np.log(0.020744)
# # sigma = np.log(4.33924339)
# # x = np.linspace(min(bins), max(bins), 20000)
# # pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
# #        / (x * sigma * np.sqrt(2 * np.pi)))
# # ax1.plot(x, pdf, linewidth = 2, linestyle = '--', c='red', label = 'Log-Normal Distribution')
# # ax1.set_xlabel('Relative Amplitude', fontdict = font)
# # ax1.set_ylabel('Probability Density', fontdict=font)
# # L = ax1.legend(prop={'size':11})
# # plt.setp(L.texts, family='Serif')
# # ax1.text(0.15, 28, r'$p(x)= \frac{1}{\sigma x \sqrt{2\pi}}\, \exp\left[-\frac{(\log(x)-\mu)^{2}}{2\sigma^{2}}\right]$', size=12)
# # ax1.text(0.27, 22, r'$\sigma = 1.5$', size = 12)
# # ax1.text(0.27, 24, r'$\mu = -3.9$', size = 12)

# ax2.hist(np.log10(M), alpha = 0.5, bins = 15, label = 'M Type')
# ax2.hist(np.log10(K), alpha = 0.75, bins = 15, label = 'K Type')
# ax2.hist(np.log10(G), alpha = 0.75, bins = 15, label = 'G Type')
# ax2.hist(np.log10(F), alpha = 0.75, bins = 15, label = 'F Type')
# ax2.hist(np.log10(OBA), alpha = 1, bins = 15, label = 'O/B/A Types')
# ax2.set_xlabel('$\log_{10} (\mathrm{Flare\;Rate\;(flares\cdot min^{-1}}))$', fontdict = font)
# ax2.set_ylabel('Count', fontdict = font)
# L2 = ax2.legend(prop={'size':10}, ncol=2)
# plt.setp(L2.texts, family='Serif')
# fig.tight_layout()
# fig.savefig('C:/Users/Nathan/OneDrive - Washington University in St. Louis/Desktop/Flare_Sim_Rate_Stats.png', dpi=1000)
# fig.show()