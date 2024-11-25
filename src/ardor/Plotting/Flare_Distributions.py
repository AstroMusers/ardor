# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 10:29:55 2024

@author: Nate Whitsett
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }
font_small ={'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
data = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Literature_Catalogs/All_Exoplanet_MCMC_Flares.csv')
data2 = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Literature_Catalogs/All_TOI_MCMC_Flares.csv")
data3 = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Literature_Catalogs/Gunther_2020.csv')
data4 = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Literature_Catalogs/Illin_2024.csv')
data5 = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Literature_Catalogs/Pietras_2020.csv')

amp = np.array(data['Amplitude'])
amp = amp[~np.isnan(amp)]
amp2 = np.array(data2['Amplitude'])
amp2 = amp2[~np.isnan(amp2)]
amp3 = np.array(data3['Amp'])
amp3 = amp3[~np.isnan(amp3)]
amp4 = np.array(data4['Amp_TESS'])
amp4 = amp4[~np.isnan(amp4)]
amp5 = np.array(data5['Amp'])
amp5 = amp5[~np.isnan(amp5)]
amp5 = [i for i in amp5 if i != 0]

FWHM = np.array(data['FWHM'])*24*60
FWHM = FWHM[~np.isnan(FWHM)]
FWHM2 = np.array(data2['FWHM'])*24*60
FWHM2 = FWHM2[~np.isnan(FWHM2)]
FWHM3 = np.array(data3['FWHM'])*24*60
FWHM3 = FWHM3[~np.isnan(FWHM3)]
FWHM4 = (np.array(data4['tstop_TESS']) - np.array(data4['tstart_TESS']))*24*60*.68

fig, ax= plt.subplots(nrows = 1, ncols = 2)
fig.set_size_inches(10,5)

count, bins, ignored = ax[0].hist(amp, bins = 100, range=(0,0.1), density = True, align = 'mid', alpha = 0.25, color= 'r')
count, bins, ignored = ax[0].hist(amp2, bins = 75, range=(0,0.1), density = True, align = 'mid', alpha = 0.25, color= 'blue')
# count, bins, ignored = ax[0].hist(amp3, bins = 300, range=(0,0.5), density = True, align = 'mid', alpha = 0.4, color= 'black', label = 'Gunther 2020')
# count, bins, ignored = ax[0].hist(amp5, bins = 300, range=(0,0.5), density = True, align = 'mid',alpha = 0.25, color= 'orange', label = 'Pietras 2022')
# count, bins, ignored = ax[0].hist(amp4, bins = 300, range=(0,0.5), density = True, align = 'mid',alpha = 0.25, color= 'blue', label = 'Illin 2024')

count, bins, ignored = ax[1].hist(FWHM, bins = 50, range=(0,20), density = True, align = 'mid', alpha = 0.75, color= 'r', label = 'Exoplanet Hosts')
count, bins, ignored = ax[1].hist(FWHM2, bins = 50, range=(0,20), density = True, align = 'mid', alpha = 0.6, color= 'cyan', label = 'TOI Hosts')
count, bins, ignored = ax[1].hist(FWHM3, bins = 50, range=(0,20), density = True, align = 'mid', alpha = 0.4, color= 'black', label = 'Gunther 2020')
count, bins, ignored = ax[1].hist(FWHM4, bins = 15, range=(0,20), density = True, align = 'mid', alpha = 0.4, color= 'orange', label = 'Illin 2024*')

mu = np.mean(np.log(amp))
sigma = np.std(np.log(amp))
print(mu, sigma, ' Whitsett_Exo')
x = np.linspace(0, 0.05, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[0].plot(x, pdf, linewidth=2, c= 'r', label = 'Exoplanet Hosts')


mu = np.mean(np.log(amp2))
sigma = np.std(np.log(amp2))
print(mu, sigma, ' Whitsett_TOI')
x = np.linspace(0, 0.05, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[0].plot(x, pdf, linewidth=2,c= 'blue', label = 'TOI Hosts')


mu = np.mean(np.log(amp3))
sigma = np.std(np.log(amp3))
print(mu, sigma, ' Gunther')
x = np.linspace(0, 0.05, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[0].plot(x, pdf, linewidth=2,  linestyle = ':',c= 'black', label = 'Gunther 2020')


mu = np.mean(np.log(amp4))
sigma = np.std(np.log(amp4))
print(mu, sigma, 'Illin')
x = np.linspace(0, 0.5, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[0].plot(x, pdf, linewidth=2,  linestyle = ':',c= 'green', label = 'Illin 2024')


mu = np.mean(np.log(amp5))
sigma = np.std(np.log(amp5))
print(mu, sigma, ' Pietras')
x = np.linspace(0, 0.5, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[0].plot(x, pdf, linewidth=2,  linestyle = ':', c = 'cyan', label = 'Pietras 2022')

mu = np.mean(np.log(FWHM))
sigma = np.std(np.log(FWHM))
x = np.linspace(0, 0.5, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[1].plot(x, pdf, linewidth=2, linestyle = '--', c= 'r')


mu = np.mean(np.log(FWHM2))
sigma = np.std(np.log(FWHM2))
x = np.linspace(0, 0.5, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[1].plot(x, pdf, linewidth=2,  linestyle = '--',c= 'cyan')

mu = np.mean(np.log(FWHM3))
sigma = np.std(np.log(FWHM3))
x = np.linspace(0, 0.5, 20000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        / (x * sigma * np.sqrt(2 * np.pi)))
ax[1].plot(x, pdf, linewidth=2,  linestyle = '--',c= 'black')



ax[0].set_xlabel('Rel. Amplitude', fontdict = font)
ax[0].set_ylabel('Density', fontdict = font)
ax[1].set_xlabel('FWHM (Min.)', fontdict = font)
ax[1].set_ylabel('Density', fontdict = font)
ax[0].set_xlim(-.002, 0.05)
ax[0].legend(prop={"family":"serif"})
ax[1].legend(prop={"family":"serif"})
ax[1].set_xlim(0,20)
ax[0].set_xticklabels(ax[0].get_xticklabels(), fontdict = font_small)
ax[0].set_yticklabels(ax[0].get_yticklabels(), fontdict = font_small)
ax[1].set_xticklabels(ax[1].get_xticklabels(), fontdict = font_small)
ax[1].set_yticklabels(ax[1].get_yticklabels(), fontdict = font_small)
plt.show()