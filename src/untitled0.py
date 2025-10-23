# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 14:14:16 2025

@author: natha
"""

import ardor.Utils.Utils as U


U.transfer_columns_between_files("C:/Users/natha/Downloads/PSCompPars_2025.09.11_12.13.42.csv", 
                "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Alfven_Catalog.csv", 
                "Host_ID", ['pl_orblper', 'pl_orblpererr1', 'pl_orblpererr2'], "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Alfven_Catalog2.csv")