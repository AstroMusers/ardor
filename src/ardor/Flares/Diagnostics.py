# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 15:43:08 2024

@author: Nate Whitsett
"""

import ardor.Flares.Flare as Flares
import numpy as np
import os
from matplotlib import pyplot as plt
directory = 'C:/Users/Nate Whitsett/Desktop/Exoplanet_Hosts/Hosts/AUMic'

data = os.listdir(directory)[0]
lc = Flares.TESS_data_extract(directory + '/' + data, PDCSAP_ERR=True)
lc = Flares.tier0(directory + '/' + data, scale=401)
flares, lengths = Flares.tier1(lc.detrended_flux, 3, fast=lc.fast_bool)
# for flaress in flares:
#     plt.axvline(x=time[flaress], color='red')
#     plt.plot(time[flaress-50:flaress+50], detrend_flux[flaress-50:flaress+50])
#     plt.show()
#     input('Flare # ' + str(flaress))
#     plt.clf()
a = Flares.tier2(lc.time, lc.detrended_flux, lc.error, flares, lengths, chi_square_cutoff=15, csv=False)
plt.plot(lc.time, lc.trend)
# plt.plot(time, detrend_flux)
print(len(flares) - len(a))
print(a)
# plt.plot(time[flares[0]-50:flares[0]+50], detrend_flux[flares[0]-50:flares[0]+50])