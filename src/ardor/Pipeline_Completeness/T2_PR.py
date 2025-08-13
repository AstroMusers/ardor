# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:55:50 2024

@author: whitsett.n
"""

from ardor.Flares.Flare import TESS_data_extract
from ardor.Pipeline_Completeness.Injection_Recovery import Tier2_PR_Test
import os
import numpy as np
import pandas as pd
from scipy.integrate import simpson
files = os.listdir('C:/Users/whitsett.n/Desktop/COCONUTS-2b')
folder = 'C:/Users/whitsett.n/Desktop/COCONUTS-2b/'
epoch_list = pd.read_csv('C:/Users/whitsett.n/Desktop/COCONUTS-2A_Flare_Epochs.csv')
epoch_list = epoch_list['Epoch'].to_list()
# Tier2_PR_Test(folder, 'C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/PR_T2.csv', [0.1, 0.25, 0.5, 1, 5, 10, 20], epoch_list)

data = pd.read_csv('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/PR_T2.csv')

Precision = data['Precision']
Recall = data['Recall']

print(simpson(Precision, x = Recall, dx = 0.1))
