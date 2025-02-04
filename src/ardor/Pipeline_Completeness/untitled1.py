# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:05:53 2024

@author: whitsett.n
"""

from ardor.SPI_Forward_Models.SPI_Simulation import SPI_Cubic, SPI_sigma
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
B = [1, 10, 100]
e = [0, 0.05, 0.2, 0.5]
sigma = [5, 1, 0.5, 0.1, 0.05]
G_star = OML.Star(0.896, 1, 1, radius = 0.92165, age = 4.6e9, B = 1, alfven=0.1)
M_star = OML.Star(0.59, 1, 1, radius = 0.603, age = 4.6e9, B = 100, alfven=0.1)
stars = [M_star, G_star]
data = pd.DataFrame()
a = ['M', 'G']
for index, star in enumerate(stars):
    for es in e:
        for values in B:
            planet = OML.Planet(1, None, 0.05, es, values, orbit_length = 3000, Star = star)
            model, x = SPI_Cubic(M_star, planet, 1)
            plt.plot(x, model)
            # print(len(x), len(model), values, es)
            # data[a[index] + '_' + str(es) + '_' + str(values) + '_x'] = pd.Series(x)
            # data[a[index] + '_' + str(es) + '_' + str(values) + '_y'] = pd.Series(model)

for sigmas in sigma:
    model, x = SPI_sigma(sigmas, 1000, 0.3)
    # data[str(sigmas) + '_x'] = pd.Series(x)
    # data[str(sigmas) + '_y'] = pd.Series(model)
    
# data.to_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/SPI_Models.csv')