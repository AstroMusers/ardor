# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 10:00:26 2024

@author: natha
"""

import pandas as pd
import numpy as np
target_list = pd.read_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Spring 2024/FINESST 2024/Data_Availability.csv")

T = target_list['st_teff']
spec_type = []
for index in range(len(T)):
    temps = T[index]
    if temps < 4000 and temps > 2000:
        spec_type.append('M')
    elif temps >= 4000 and temps < 5240:
        spec_type.append('K')
    elif temps >=5240 and temps < 6050:
        spec_type.append('G')
    elif temps >= 6050 and temps < 7350:
        spec_type.append('F')
    elif temps >= 7350 and temps < 9600:
        spec_type.append('A')
    elif temps >= 9600 and temps < 29200:
        spec_type.append('B')
    elif temps >= 37800 and temps < 100000:
        spec_type.append('O')
    else:
        spec_type.append('M')

target_list['Spectral_Type'] = spec_type
target_list.to_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Spring 2024/FINESST 2024/Data_Availability.csv', index=False)