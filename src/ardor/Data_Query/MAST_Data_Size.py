# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 18:46:37 2024

@author: natha
"""

from astroquery.mast import Observations
import pandas as pd
import numpy as np
target_list = pd.read_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Spring 2024/FINESST 2024/PSCompPars_2024.02.04_17.33.48.csv")
RA = target_list['ra']
DEC = target_list['dec']
GALEX_data = []
HST_data = []
Kepler_data = []
TESS_data= []
for index in range(len(RA)):
        kepler_count = 0
        try:
            obs_table = Observations.query_criteria(s_ra = [float(RA[index]) - 0.02, float(RA[index])+0.02], s_dec = [float(DEC[index]) - 0.02, float(DEC[index])+0.02],dataproduct_type = ["TIMESERIES"], obs_collection='TESS')
            data_products = Observations.get_product_list(obs_table)
            TESS_data.append(np.log(np.sum(data_products['size'])/(1e6)))
        except:
            TESS_data.append(0)
        try:
            obs_table = Observations.query_criteria(s_ra = [float(RA[index]) - 0.02, float(RA[index])+0.02], s_dec = [float(DEC[index]) - 0.02, float(DEC[index])+0.02],dataproduct_type = ["TIMESERIES"], obs_collection='Kepler')
            data_products = Observations.get_product_list(obs_table)
            kepler_count += np.log(np.sum(data_products['size'])/(1e6))
        except:
            kepler_count += 0
        try:
            obs_table = Observations.query_criteria(s_ra = [float(RA[index]) - 0.02, float(RA[index])+0.02], s_dec = [float(DEC[index]) - 0.02, float(DEC[index])+0.02],dataproduct_type = ["TIMESERIES"], obs_collection='K2')
            data_products = Observations.get_product_list(obs_table)
            kepler_count += np.log(np.sum(data_products['size'])/(1e6))
            Kepler_data.append(kepler_count)
        except:
            Kepler_data.append(kepler_count)
        try:
            obs_table = Observations.query_criteria(s_ra = [float(RA[index]) - 0.5, float(RA[index])+0.5], s_dec = [float(DEC[index]) - 0.5, float(DEC[index])+0.5],obs_collection='GALEX')
            data_products = Observations.get_product_list(obs_table)
            GALEX_data.append(np.log(np.sum(data_products['size'])/(1e6)))
        except:
            GALEX_data.append(0)
        try:
            obs_table = Observations.query_criteria(s_ra = [float(RA[index]) - 0.02, float(RA[index])+0.02], s_dec = [float(DEC[index]) - 0.02, float(DEC[index])+0.02], dataproduct_type = ["TIMESERIES"], obs_collection='HST')
            data_products = Observations.get_product_list(obs_table)
            HST_data.append(np.log(np.sum(data_products['size'])/(1e6)))
        except:
            HST_data.append(0)
        print(len(HST_data), len(Kepler_data), len(GALEX_data), len(TESS_data))
           
print(TESS_data, Kepler_data, GALEX_data, HST_data)
target_list['HST_data'] = HST_data
target_list['GALEX_data'] = GALEX_data
target_list['Kepler_data'] = Kepler_data
target_list['TESS_data'] = TESS_data
target_list.to_csv("C:/Users/natha/OneDrive - Washington University in St. Louis/Grad School/Spring 2024/FINESST 2024/Data_Availability.csv", index=False)