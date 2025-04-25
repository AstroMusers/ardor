<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:34:09 2024

@author: natha
"""
import pandas as pd
import numpy as np
import warnings
import os
import ardor.Flares.Flare as Flare
warnings.filterwarnings("ignore")
def Data_Transfer(source_file, output_file, ID_column_header, column_headers=[], output_dir = None, specifier_column = None):
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
    source = pd.read_csv(source_file)
    output_data = pd.read_csv(output_file)
    IDs = output_data[ID_column_header]
    for header in column_headers:
        data_list = []
        for index, ID in enumerate(IDs):
            if specifier_column == None:
                data_2_transfer = source.loc[(source[ID_column_header] == str(ID)),
                                             str(header)]
                if len(data_2_transfer) > 1 and specifier_column == None:
                    print('Your list is degenerate! Additional specifier needed')
                    return None
                else:
                    try:
                        data_2_transfer = data_2_transfer.item()
                    except:
                        data_2_transfer= np.nan
            if specifier_column != None:
                data_2_transfer = source.loc[((source[ID_column_header] == str(ID)) &
                                              (source[specifier_column] == 'b')),
                                             str(header)]
                try:
                    data_2_transfer = data_2_transfer.item()
                except:
                    data_2_transfer= np.nan

            data_list.append(data_2_transfer)
            
        output_data[header] = data_list
    print(output_data)
    if output_dir == None:
        new_output = output_file[0:] + '2' + output_file[:0]
        new_output.to_csv(new_output, index=False)
    elif output_dir != None:
        output_data.to_csv(output_dir, index=False)
def range_shift(data, a, b, x, y):
    '''
    Generates a linear map to take data in range [a,b] to new range of [x,y].

    Parameters
    ----------
    data : float
        Data to linearly map.
    a : float
        Lower bound of original data.
    b : float
        Upper bound of original data.
    x : float
        Lower bound of new data.
    y : float
        Upper bound of new data.

    Returns
    -------
    new_data: float
        Linearly transformed data in new range [x,y].

    '''
    return x + (y-x)/(b-a)*(data-a)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
# source_file = "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Alfven_Parameters/Alfven_Catalog.csv"
# output_file = "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares.csv"
# data = pd.read_csv(output_file)
# hosts = set(data['Host_ID'])
# folders  = os.listdir('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs/')
# output_list = []
# for host in hosts:
#     obs_time = 0
#     count = 0
#     fast_count = 0
#     for folder in folders:
#         if folder == str(host):
#             files = os.listdir("C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs/" + folder)
#             for file in files:
#                 if file.endswith('s_lc.fits'):
#                     count += 1
#             obs_time += (count)*26
#             output_list.append(obs_time)
# output_list = np.array(output_list)
# print(np.mean(output_list))
# data['Obs_Time'] = output_list
# data.to_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_TOI_MCMC_Flares_New.csv')
                
# Data_Transfer(source_file, output_file, 'Host_ID', column_headers = ["st_met", "st_age", "st_rotp"], output_dir="C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/New_file.csv", specifier_column='pl_letter')
=======
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:34:09 2024

@author: natha
"""
import pandas as pd
import numpy as np
import warnings
import os
import ardor.Flares.Flare as Flare
warnings.filterwarnings("ignore")
def Data_Transfer(source_file, output_file, ID_column_header, column_headers=[], output_dir = None, specifier_column = None):
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
    source = pd.read_csv(source_file)
    output_data = pd.read_csv(output_file)
    IDs = output_data[ID_column_header]
    for header in column_headers:
        data_list = []
        for index, ID in enumerate(IDs):
            if specifier_column == None:
                data_2_transfer = source.loc[(source[ID_column_header] == str(ID)),
                                             str(header)]
                if len(data_2_transfer) > 1 and specifier_column == None:
                    print('Your list is degenerate! Additional specifier needed')
                    return None
                else:
                    try:
                        data_2_transfer = data_2_transfer.item()
                    except:
                        data_2_transfer= np.nan
            if specifier_column != None:
                data_2_transfer = source.loc[((source[ID_column_header] == str(ID)) &
                                              (source[specifier_column] == 'b')),
                                             str(header)]
                try:
                    data_2_transfer = data_2_transfer.item()
                except:
                    data_2_transfer= np.nan

            data_list.append(data_2_transfer)
            
        output_data[header] = data_list
    print(output_data)
    if output_dir == None:
        new_output = output_file[0:] + '2' + output_file[:0]
        new_output.to_csv(new_output, index=False)
    elif output_dir != None:
        output_data.to_csv(output_dir, index=False)
def range_shift(data, a, b, x, y):
    '''
    Generates a linear map to take data in range [a,b] to new range of [x,y].

    Parameters
    ----------
    data : float
        Data to linearly map.
    a : float
        Lower bound of original data.
    b : float
        Upper bound of original data.
    x : float
        Lower bound of new data.
    y : float
        Upper bound of new data.

    Returns
    -------
    new_data: float
        Linearly transformed data in new range [x,y].

    '''
    return x + (y-x)/(b-a)*(data-a)


# source_file = "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Alfven_Parameters/Alfven_Catalog.csv"
# output_file = "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares.csv"
# data = pd.read_csv(output_file)
# hosts = set(data['Host_ID'])
# folders  = os.listdir('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs/')
# output_list = []
# for host in hosts:
#     obs_time = 0
#     count = 0
#     fast_count = 0
#     for folder in folders:
#         if folder == str(host):
#             files = os.listdir("C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOIs/" + folder)
#             for file in files:
#                 if file.endswith('s_lc.fits'):
#                     count += 1
#             obs_time += (count)*26
#             output_list.append(obs_time)
# output_list = np.array(output_list)
# print(np.mean(output_list))
# data['Obs_Time'] = output_list
# data.to_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_TOI_MCMC_Flares_New.csv')
                
# Data_Transfer(source_file, output_file, 'Host_ID', column_headers = ["st_met", "st_age", "st_rotp"], output_dir="C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/New_file.csv", specifier_column='pl_letter')
>>>>>>> 843bcfd59f29e0297b65b05360ddcf3ed4f7ad77
    