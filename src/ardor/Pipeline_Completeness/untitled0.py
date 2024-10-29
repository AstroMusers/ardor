# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:49:05 2024

@author: Nate Whitsett
"""

from astroquery.mast import Catalogs
from ardor.Data_Query.Data_Query import Bulk_TESS_lc_Query
import numpy as np
results = Catalogs.query_criteria(catalog="Tic", Teff = (5000, 5010), objType="STAR", ra = (270, 272), dec=(70, 72), radius = 1)

print(np.array(results['ra']), np.array(results['dec']), np.array(results['ID']))


Bulk_TESS_lc_Query(results['ra'], results['dec'], results['ID'], 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/G_type_SPI_Sim/TESS_LC_Files', results['ID'])
