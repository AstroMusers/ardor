# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 13:33:18 2024

@author: whitsett.n
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
data = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Alfven_Parameters/Alfven_Catalog_2.csv')
r_close = data['Column1'].to_list()
r_close_p = data['r_close_+'].to_list()
r_close_m = data['r_close_-'].to_list()
alfven = data['Alfven_Rad'].to_list()
alfven_p = data['0.13'].to_list()
alfven_m = data['-0.135'].to_list()
def sampler(value1, plus1, minus1, value2, plus2, minus2):
    count = 0
    for samples in range(10000):
        rand1 = np.random.random()
        if rand1 > 0.5:
            draw = np.abs(np.random.normal(loc=0, scale = np.abs(plus1)))
            value_1_samp = value1 + draw
        elif rand1 < 0.5:
            draw = np.abs(np.random.normal(loc=0, scale = np.abs(minus1)))
            value_1_samp = value1 - draw
        rand2 = np.random.random()
        if rand2 > 0.5:
            draw = np.abs(np.random.normal(loc=0, scale = np.abs(plus2)))
            value_2_samp = value2 + draw
        elif rand2 < 0.5:
            draw = np.abs(np.random.normal(loc=0, scale = np.abs(minus2)))
            value_2_samp = value2 - draw
        if value_1_samp < value_2_samp:
            count += 1
        elif value_1_samp > value_2_samp:
            count += 0
    return count/10000
probability_list = []
for index, values in enumerate(r_close):
    probability_list.append(sampler(r_close[index], r_close_p[index], r_close_m[index], alfven[index], alfven_p[index], alfven_m[index]))

data['P(rp < RA)'] = probability_list
data.to_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Alfven_Parameters/Alfven_Catalog_3.csv')