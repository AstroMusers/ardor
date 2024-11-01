# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:58:06 2024

@author: Nate Whitsett
"""

from astroquery.mast import Catalogs
from ardor.Data_Query.Data_Query import Bulk_TESS_lc_Query
import numpy as np
results = Catalogs.query_criteria(catalog="CTL", Teff = (3500, 4000), objType="STAR", ra = (281, 282), dec=(65, 66), radius = 0.5)

print(np.array(results['ra']), np.array(results['dec']), np.array(results['ID']))


Bulk_TESS_lc_Query(results['ra'], results['dec'], results['ID'], 'C:/Users/Nate Whitsett/Desktop/M_Type_LC/', results['ID'], radius=0.05)
