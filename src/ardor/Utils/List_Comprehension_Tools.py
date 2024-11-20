# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:34:09 2024

@author: natha
"""
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
def Data_Transfer(input_file, output_file, column_name_list=[]):
    '''
    Transfer host parameters based on host name.

    Parameters
    ----------
    input_file : str
        Input csv file directory.
    output_file : TYPE
        Output csv file directory.

    Returns
    -------
    None.

    '''
    data = pd.read_csv(input_file)
    output = pd.read_csv(output_file)
    hosts = data['hostname'].to_list()
    counter = 0
    for indices, columns in enumerate(column_name_list):
        counter = 0

        for index, names in enumerate(hosts):
            value = data.loc[(data['hostname'] == names)][columns].to_list()
            df = output.loc[output['Host_ID'] == names]
            df[columns] = np.ones(len(df['Host_ID']))*value
            if len(df) == 0:
                continue
            if counter == 0:
                output_df = df
                counter += 1
            output_df = pd.concat([output_df,df])
            output = output_df
    print(output)
    output_df.to_csv(output_file, index = False)
            
            
        
columns = ['pl_tranmid', 'pl_tranmiderr1', 'pl_tranmiderr2', 'pl_orbtper', 'pl_orbtpererr1', 'pl_orbtpererr2', 'r_close', 'pl_orbper']
Data_Transfer('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/TESS Data/All_Hosts.csv', 'C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All_Exoplanets/All_Exoplanet_MCMC_Flares.csv', columns)