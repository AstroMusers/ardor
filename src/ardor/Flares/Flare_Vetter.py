# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 23:12:53 2024

@author: Nate Whitsett
"""

import pandas as pd
import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pdf2image import convert_from_path

all_exo = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Final Data/All_Exoplanet_MCMC_Flares.csv', delimiter = ',')

hosts = all_exo['Host_Name']
dirs = os.listdir('C:/Users/natha/Desktop/All_Exoplanet')
index_track = 0

for folders in dirs:
    hosts = []
    flare_no = []
    vet = []
    counter = 1
    host = folders.replace(" ", "")
    sub_dir = os.listdir('C:/Users/natha/Desktop/All_Exoplanet/' + folders)
    print(index_track)
    for files in sub_dir:
        try:
            if files[0:14] == 'mcmc_fit_Flare' and files.endswith('.pdf'):
                images = convert_from_path('C:/Users/natha/Desktop/All_Exoplanet/' + folders + '/mcmc_fit_Flare' + str(counter) + '.pdf')
                images = images[0]
                images.save('C:/Users/natha/Desktop/All_Exoplanet/' + folders + '/mcmc_fit_Flare' + str(counter) + '.png')
                img = mpimg.imread('C:/Users/natha/Desktop/All_Exoplanet/' + folders + '/mcmc_fit_Flare' + str(counter) + '.png')
                imgplot = plt.imshow(img)
                plt.show()
                x = input('Vet? (1 good, 2 inbetween, 3 bad)')
                if x == 'skip':
                    continue
                hosts.append(folders)
                flare_no.append(counter)
                vet.append(x)
                counter += 1
        except:
            hosts.append(folders)
            flare_no.append(counter)
            vet.append(3)
            print(1)
            counter += 1
            continue
    X = np.stack((hosts,flare_no,vet)).T
    with open('C:/Users/natha/Desktop/All_Exoplanet_Vet.csv', "a") as f:
        np.savetxt(f, X, delimiter=",", fmt='%s')
        f.close()
    counter = 0
    index_track += 1