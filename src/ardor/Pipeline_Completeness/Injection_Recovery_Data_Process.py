# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 18:52:27 2023

@author: Nate Whitsett
"""
import numpy as np
import pandas as pd
import csv
from matplotlib import pyplot as plt
data = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/Params/Injection_Recovery_M_Type_T1.csv')
amp = list(data['Amplitude'])
FWHM = list(data['FWHM'])
error = list(data['Error'])
Bool = list(data['Accepted?'])
integral = list(data['Integral'])
amp_bins = np.logspace(-3, 0, num=17)
FWHM_bins = np.linspace(0.001388888, 0.041, num=17)
error_bins = np.linspace(0.0001, 0.25, num=16)
print(amp_bins)
x = []
y = []

for i in range(len(amp_bins)-1):
    tmp = []
    for index in range(len(amp)):
        if amp[index] > amp_bins[i] and amp[index] < amp_bins[i+1]:
            tmp.append([integral[index], FWHM[index], Bool[index]])
    x.append(tmp)
total = 0
for i in range(len(FWHM_bins)):
    tmp = []
    for cells in x:
        count = 0
        pos = 0
        for flares in cells:
            if flares[1] > FWHM_bins[i] and flares[1] < FWHM_bins[i+1]:
                count += 1
                total += 1
                if flares[2] == 1:
                    pos += 1
        if i == 16:
            continue
        if count == 0:
            tmp.append(np.nan)
        elif count != 0:
            tmp.append(pos/count * 100)
    y.append(tmp)

with open('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Injection Tests/Params/Injection_Recovery_M_Type_T1_Grid.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(y)