# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 13:58:23 2025

@author: whitsett.n
"""

import pandas as pd
import numpy as np

periods = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Host_Period.csv")
catalog = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Alfven_Catalog.csv")
period_list = []
for host in catalog['Host_ID']:
    st_per = periods.loc[periods['hostname'] == str(host), 'st_rotp']
    if st_per.empty == True:
        period_list.append(np.nan)
    else:
        period_list.append(np.array(st_per)[0])
catalog['st_rotp2'] = period_list

catalog.to_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Alfven_Catalog.csv", index = False)