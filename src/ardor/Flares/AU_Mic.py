# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 23:31:02 2024

@author: natha
"""

import ardor.Flares.Flare as Flare
import os
from matplotlib import pyplot as plt
data = os.listdir('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Hosts/AUMic')
base = 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Hosts/AUMic/'
for files in data:
    lc = Flare.tier0(base + files)
    flares, lengths = Flare.tier1(lc.detrended_flux, 2.5)
    plt.plot(lc.time, lc.detrended_flux)
    Flare.tier2(lc.time,lc.detrended_flux, lc.error, flares, lengths, chi_square_cutoff=20, output_dir = 'C:/Users/whitsett.n/Desktop', catalog_name='AU_Mic.csv', csv=True)