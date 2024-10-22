# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 15:43:08 2024

@author: Nate Whitsett
"""

import ardor.Flares.Flare as Flares
import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
parameters = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/All_Hosts.csv')
directory = 'C:/Users/Nate Whitsett/Desktop/Exoplanet_Hosts/Hosts/'
host_list = os.listdir(directory)
flare_dir = 'C:/Users/Nate Whitsett/Desktop/Exoplanet_Hosts/Flares/'

for hosts in host_list:
    print(hosts)
    error_list = []
    try:
        input_dir = directory + str(hosts)
        output_dir = flare_dir
        try:
            os.mkdir(flare_dir + str(hosts), mode=0o777)
        except FileExistsError:
            os.mkdir(flare_dir + str(hosts) + '(1)', mode=0o777)
        file_list = os.listdir(input_dir)
        for index, files in enumerate(file_list):
            print('LC #: ' + str(index) )
            Teff = float(parameters.loc[parameters['hostname'] == hosts]['st_teff'].item())
            period = float(parameters.loc[parameters['hostname'] == hosts]['pl_orbper'].item())
            radius = float(parameters.loc[parameters['hostname'] == hosts]['st_rad'].item())
            epoch = float(parameters.loc[parameters['hostname'] == hosts]['pl_orbtper'].item())
            if np.isnan(epoch) == True:
                epoch = float(parameters.loc[parameters['hostname'] == hosts]['pl_tranmid'].item())
                if np.isnan(epoch) == True:
                    epoch = np.nan      
            lc = Flares.tier0(input_dir + '/' + str(files))
            flares = Flares.tier1(lc.detrended_flux, 2.5)
            ZZ = Flares.tier2(lc.time, lc.flux, lc.error, flares.index, flares.length,
                          chi_square_cutoff=15, output_dir=output_dir, host_name=hosts,
                          T = Teff, host_radius=radius, csv=True, planet_period=period, 
                          planet_epoch=epoch)
            with open('C:/Users/Nate Whitsett/Desktop/Exoplanet_Hosts/T2_Host_Flare_Catalog.csv', "a") as f:
                np.savetxt(f, ZZ, delimiter=",", fmt='%s')
                f.close()
    except:
        error_list.append(hosts)
        print('Error in pipeline! Consider rerunning: ' + str(hosts))
        continue
print(error_list)