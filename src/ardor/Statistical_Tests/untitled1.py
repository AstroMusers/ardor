# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 14:12:29 2025

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import ardor.Statistical_Tests.MLE as MLE
import numpy as np



for loc in np.linspace(0, 1, num=50):
    phases = np.abs(np.random.normal(loc,scale=0.1, size=(15)))
    print(phases)
beta = MLE.SPI_Metric(3, .6, 0.03, phases, 200, 0.8)
print(beta)