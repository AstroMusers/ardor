# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:34:09 2024

@author: natha
"""
import pandas as pd
import numpy as np
import warnings
import math
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

def df_return(df, host, columns, planet = 'b'):
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

def boudnary_conditions(val, lower, upper):
    if val < lower:
        val += upper
    if val > upper:
        val -= upper
    return val