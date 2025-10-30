# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:34:09 2024

@author: natha
"""
import pandas as pd
import numpy as np
import warnings
import math
import os
from astropy.table import Table
from astropy.io import ascii
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

def asymmetric_sample(value, upper, lower):
    upper = np.abs(upper)
    lower = np.abs(lower)
    check = np.random.random()
    if check > 0.5:
        change = np.abs(np.random.normal(loc = 0, scale = upper))
        sample = check + change
    elif check <= 0.5:
        change = np.abs(np.random.normal(loc = 0, scale = lower))
        sample = check - change
    return sample

def df_return(df, host, columns, planet = None):
    '''
    Return entries in a column based on host name.

    Parameters
    ----------
    host : str
        Name of host.
    column : str
        The columns you want to search.

    Returns
    -------
    filtered dataframe

    '''
    if planet == None:
        if len(columns) == 1:
            a = df.loc[(df['Host_ID'] == host), columns].to_numpy()
            return [x for xs in a for x in xs]
        elif len(columns) > 1:
            column_to_return = []
            for column in columns:
                column_to_return.append(df.loc[(df['Host_ID'] == host), column].to_numpy())
            return [x for xs in column_to_return for x in xs]
    else:
        if len(columns) == 1:
            return df.loc[(df['Host_ID'] == host) & (df['pl_letter'] == planet), columns].to_numpy()
        elif len(columns) > 1:
            column_to_return = []
            for column in columns:
                column_to_return.append(df.loc[(df['Host_ID'] == host) & (df['pl_letter'] == planet), column].to_numpy())
            return column_to_return

def transfer_columns_between_files(
    file_a_path,
    file_b_path,
    identifier_column,
    columns_to_transfer,
    output_file_path='file_b_updated.csv'
):
    """
    Transfers specified columns from File A to File B by matching on an identifier column.

    Args:
        file_a_path (str): Path to the source file (File A).
        file_b_path (str): Path to the destination file (File B).
        identifier_column (str): The column name used to match rows.
        columns_to_transfer (list of str): List of column names to transfer from A to B.
        output_file_path (str): Path to save the updated File B.

    Returns:
        pd.DataFrame: The updated DataFrame.
    """
    # Load files
    df_a = pd.read_csv(file_a_path)
    df_b = pd.read_csv(file_b_path)

    df_merged = df_b.copy()

    for col in columns_to_transfer:
        if col in df_b.columns:
            # Merge with suffix to avoid conflict
            temp = df_a[[identifier_column, col]]
            df_merged = df_merged.merge(
                temp,
                on=identifier_column,
                how='left',
                suffixes=('', '_from_A')
            )
            # Overwrite B's column with A's data
            df_merged[col] = df_merged[f'{col}_from_A']
            df_merged = df_merged.drop(columns=[f'{col}_from_A'])
        else:
            # Merge directly
            df_merged = df_merged.merge(
                df_a[[identifier_column, col]],
                on=identifier_column,
                how='left'
            )

    # Save result
    df_merged.to_csv(output_file_path, index=False)

    print(f'✅ Data transferred. Saved to {output_file_path}')
    return df_merged

def round_to_sig_figs(number, sig_figs):
    if number == 0:
        return 0.0
    
    # Determine the magnitude
    magnitude = int(math.floor(math.log10(abs(number))))
    
    # Calculate the scaling factor
    factor = 10**(sig_figs - 1 - magnitude)
    
    # Round and scale back
    rounded_number = round(number * factor) / factor
    return rounded_number

def find_absolute_minimum_2d(arr):
    if not arr or not arr[0]:  # Handle empty or malformed arrays
        return None

    min_value = arr[0][0]  # Initialize with the first element

    for row in arr:
        for element in row:
            if element < min_value:
                min_value = element
    return min_value

def find_absolute_maximum_2d(arr):
    if not arr or not arr[0]:  # Handle empty or malformed arrays
        return None

    max_value = arr[0][0]  # Initialize with the first element

    for row in arr:
        for element in row:
            if element > max_value:
                max_value = element
    return max_value

def boundary_conditions(val, lower, upper):
    if isinstance(val, np.ndarray) == True or isinstance(val, list) == True:
        for idx, elements in enumerate(val):
            if val[idx] < lower:
                val[idx] += upper
            if val[idx] > upper:
                val[idx] -= upper
    else:
        if val < lower:
            val += upper
        if val > upper:
            val -= upper
    return val


def add_list_to_csv(csv_path, new_column_name, values_list, output_path=None):
    """
    Add a new column to a CSV file from a list of values.
    If the CSV does not exist, create a new one.

    Parameters
    ----------
    csv_path : str
        Path to the existing CSV file (or where to create it).
    new_column_name : str
        Name of the new column to add.
    values_list : list
        The values to insert as a column.
    output_path : str, optional
        Path to save the updated CSV. If None, overwrites csv_path.

    Behavior
    --------
    - If the CSV exists:
        * If the list is shorter, pads with NaN.
        * If the list is longer, truncates.
    - If the CSV does not exist:
        * Creates a new CSV with just the given column.
    """

    if os.path.exists(csv_path):
        # Load the CSV
        df = pd.read_csv(csv_path)

        # Adjust list length to match DataFrame
        n = len(df)
        if len(values_list) < n:
            values_list = values_list + [None] * (n - len(values_list))
        elif len(values_list) > n:
            values_list = values_list[:n]

        # Add the new column
        df[new_column_name] = values_list

    else:
        # Create a new DataFrame with just this column
        df = pd.DataFrame({new_column_name: values_list})

    # Save back to file
    if output_path is None:
        output_path = csv_path
    df.to_csv(output_path, index=False)

def construct_param_file(param_dir, priors = None, baseline = 'hybrid_spline', 
                         model = 'flare', N=1, name = 'Flare'):
    if model == 'flare':
        if N == 1:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1']
            if priors == None:
                best_guess = [0, 0.05, 0.1]
                priors = ['uniform -1 1', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare_Time', 'Flare_FWHM', 'Flare_Amp.']
            units = ['days', 'days', 'rel. flux']
        elif N == 2:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                           'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2']
            if priors == None:
                best_guess = [0, 0.05, 0.1, 0.2, 0.05, 0.1]
                priors = ['uniform -1 1', 'uniform 0 0.5', 'uniform 0 3', 
                          'uniform -1 1', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.']
            units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux']
        elif N == 3:
            param_names = ['flare_tpeak_1', 'flare_fwhm_1', 'flare_ampl_1',
                           'flare_tpeak_2', 'flare_fwhm_2', 'flare_ampl_2',
                           'flare_tpeak_3', 'flare_fwhm_3', 'flare_ampl_3']
            if priors == None:
                best_guess = [0, 0.05, 0.1, 0.2, 0.05, 0.1, 0.4, 0.05, 0.1]
                priors = ['uniform -0.5 0.5', 'uniform 0 0.5', 'uniform 0 3', 
                          'uniform 0 0.5', 'uniform 0 0.5', 'uniform 0 3',
                          'uniform 0 0.5', 'uniform 0 0.5', 'uniform 0 3']
                labels = ['Flare1_Time', 'Flare1_FWHM', 'Flare1_Amp.',
                            'Flare2_Time', 'Flare2_FWHM', 'Flare2_Amp.',
                            'Flare3_Time', 'Flare3_FWHM', 'Flare3_Amp.']
            units = ['days', 'days', 'rel. flux', 'days', 'days', 'rel. flux',
                        'days', 'days', 'rel. flux']
    if baseline == 'gp_matern32':
        param_names += [f'ln_err_flux_{name}', f'baseline_gp_offset_flux_{name}', f'baseline_gp_matern32_lnsigma_flux_{name}', 
                        f'baseline_gp_matern32_lnrho_flux_{name}']
        best_guess += [-7,0, -5, 0]
        priors += ['uniform -10 -2', 'uniform -0.02, 0.02', 'uniform -15 0', 'uniform -1 15']
        labels += [f'ln_err_flux_{name}', rf'GP_Matern32_Baseline_{name}', rf'GP_Matern32_ln$sigma$_{name}', 
                   rf'GP_Matern32_ln$rho$_{name}']
        units += ['rel. flux', 'rel. flux', '', '']
    elif baseline == 'hybrid_spline':
        param_names += [f'ln_err_flux_{name}']
        best_guess += [-7]
        priors += ['uniform -10 -2']
        labels += [f'ln_err_flux_{name}']
        units += ['rel. flux']
    table = Table()
    table.add_column(param_names, name='#param_name')
    table.add_column(best_guess, name='best_guess')
    table.add_column([1]*len(param_names), name='fit')
    table.add_column(priors, name='prior')
    table.add_column(labels, name='label')
    table.add_column(units, name='unit')
    ascii.write(table, os.path.join(param_dir, 'params.csv'), overwrite=True, format='csv')

def construct_settings_file(settings_file_dir, working_dir, baseline = 'hybrid_spline', 
                         N=1, name = 'Flare', multi_process_bool = True
                         , cores = 10, walkers = 15, steps = 5000, burn_in = 1000):
    """Constructs settings file for Tier 3 allesfitter MCMC run.

    Args:
        settings_file_dir (str): Directory of template settings file.
        working_dir (str): Directory to save constructed settings file.
        baseline (str, optional): Baseline you want to use. Defaults to 'hybrid_spline'. 'sample_GP_Matern32' also available.
        N (int, optional): Number of flares to fit. Defaults to 1.
        name (str, optional): Name of the fit. Defaults to 'Flare'.
        multi_process_bool (bool, optional): Whether to use multiprocessing. Defaults to True.
        cores (int, optional): Number of cores to use for multiprocessing. Defaults to 10.
        walkers (int, optional): Number of walkers for MCMC. Defaults to 15.
        steps (int, optional): Total number of steps for MCMC. Defaults to 5000.
        burn_in (int, optional): Number of burn-in steps for MCMC. Defaults to 1000.
    """
    table = ascii.read(settings_file_dir, format='csv')
    for idx, values in enumerate(table['#name']):
        if 'Flaronardo' in table['#name'][idx]:
            table['#name'][idx] = str(values).replace("Flaronardo", name)
    for idx, values in enumerate(table['value']):
        if 'Flaronardo' in table['value'][idx]:
            table['value'][idx] = str(values).replace("Flaronardo", name)

    table['value'][table['#name'] == 'inst_phot'] = name
    table['value'][table['#name'] == 'multiprocess'] = multi_process_bool
    table['value'][table['#name'] == 'multiprocess_cores'] = cores
    table['value'][table['#name'] == 'mcmc_nwalkers'] = walkers
    table['value'][table['#name'] == 'mcmc_total_steps'] = steps
    table['value'][table['#name'] == 'mcmc_burn_steps'] = burn_in
    table['value'][table['#name'] == f'baseline_flux_{name}'] = baseline
    table['value'][table['#name'] == 'N_flares'] = N
    ascii.write(table, os.path.join(working_dir, 'settings.csv'), overwrite=True, format='csv')