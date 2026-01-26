# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:06:28 2023

@author: Nathan
"""

from astroquery.mast import Observations
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive as nea
import os
import shutil
import numpy as np
import pandas as pd
from collections import namedtuple
from astroquery.mast import Catalogs
from ardor.SPI_Forward_Models.Orbit_Model_Library import transit_to_periastron_epoch
from ardor.Data_Query.Catalog import Catalog
import warnings
import pandas as pd
warnings.filterwarnings("ignore")
def Bulk_TESS_lc_Query(RA_list, DEC_list, TIC_ID_list, download_dir, host_name_list, radius = 0.01):
    '''
    

    Parameters
    ----------
    RA_list : list
        A list of Right Ascensions, in degrees, of your targets.
    DEC_list : list or array-like
        A list of Declinations, in degrees, of your targets.
    TIC_ID_list : list or array-like
        A list of TESS Input Catalog values of your targets.
    download_dir : list or array-like
        The primary directory where the data for each target will be.
    host_name_list : list or array-like
        The names of the target. These will correspond to the folder names 
        the data will be assigned to.
    radius : float, optional
        The radius, in degrees, in which to search around the specified RA 
        and DEC values. Typically, a value of 0.01 works well.

    Returns
    -------
    Will download TESS .fits file from the SPOC pipeline to the specified 
    directory.

    '''
    undownloaded = []
    dirs = []
    for index, TIC in enumerate(RA_list):
        count = 0
        download_dir = download_dir 
        try:
            os.mkdir(os.path.join(download_dir, str(host_name_list[index])))
        except FileExistsError:
            print('The directory already exists! Skipping the target. If you wish to rerun the data available for this target, delete or rename the current directory with its name. The run will continue.')
            dirs.append(os.path.join(download_dir, str(host_name_list[index])))
            continue
        try:
            obs_table = Observations.query_criteria(s_ra = [(float(RA_list[index]) - radius), (float(RA_list[index]) + radius)], 
                                                    s_dec = [(float(DEC_list[index]) - radius), (float(DEC_list[index]) + radius)],
                                                    calib_level = 3,
                                                    obs_collection = 'TESS',
                                                    dataproduct_type = 'TIMESERIES'
                                                    )
            data_products = Observations.get_product_list(obs_table)
        except:
            print('No products with the TIC ID and RA/DEC combination! Try a larger radius. There may be no data for this target, too.')
            continue
        
        for indices, items in enumerate(data_products['productFilename']):
            if str(items).endswith('s_lc.fits') == True or str(items).endswith('a_fast-lc.fits') == True:
                if int(items[24:40]) == int(TIC_ID_list[index]):
                    try:
                        count += 1
                        print(count)
                        Observations.download_products(data_products[indices], download_dir = download_dir)
                    except:
                        print('There appears to be a server error! This can happen if MAST does not respond in time. The potentially undownloaded file(s) will appear once the run is finished')
                        undownloaded.append(items)
                        continue  
        if count > 0:
            try:
                files = os.listdir(os.path.join(download_dir, 'mastDownload', 'TESS'))
            except:
                print('There appears to be a directory error! This can happen if your computer does not update its directories fast enough for it to recognize where to put the new file. The potentially undownloaded file(s) will appear once the run is finished')
                undownloaded.append(items)
                continue
            for data in files:
                if str(data).endswith('s') == True:
                    try:
                        os.rename(os.path.join(download_dir, 'mastDownload', 'TESS', str(data), str(data) + '_lc.fits'), os.path.join(download_dir, str(host_name_list[index]), str(data) + '_lc.fits'))
                    except:
                        print('Warning! Some files may have not downloaded. We skipped it for now, but check at the end for a list of potentially undownloaded files.')
                        undownloaded.append(data)
                        continue
                elif str(data).endswith('fast') == True:
                    try:
                        os.rename(os.path.join(download_dir, 'mastDownload', 'TESS', str(data), str(data) + '-lc.fits'), os.path.join(download_dir, str(host_name_list[index]), str(data) + '-lc.fits'))
                    except:
                        print('Warning! Some files may have not downloaded. We skipped it for now, but check at the end for a list of potentially undownloaded files.')
                        undownloaded.append(data)
                        continue
                    
            shutil.rmtree(os.path.join(download_dir, 'mastDownload'))
    for folders in os.listdir(download_dir):
        if len(os.listdir(os.path.join(download_dir, folders))) == 0:
            shutil.rmtree(os.path.join(download_dir, folders))
    print('The already existing directories are:', dirs)
    print('The potential undownloaded files are:', undownloaded)

def Query_Transit_Solution(identifier, table = "ps"):
    """
    Query the NASA Exoplanet Archive for transit solution parameters of planets.
    This function retrieves transit parameters (orbital period, transit mid-time, and transit duration)
    for planets associated with a given stellar identifier from the NASA Exoplanet Archive.
    Parameters
    ----------
    identifier : str, int, or float
        The stellar identifier. If numeric (int or float), it is treated as a TIC ID.
        If string, it is treated as a hostname and searched using a LIKE query.
    table : str, optional
        The NASA Exoplanet Archive table to query. Default is "ps".
        Supported values include "pscomppars" (Planetary Systems Composite Parameters)
        and "ps" (Planetary Systems). The aggregation method differs by table.
    Returns
    -------
    periods : numpy.ndarray
        Array of orbital periods in days for each planet with transit data.
    transit_mids : numpy.ndarray
        Array of transit mid-times (BJD) for each planet with transit data.
    durations : numpy.ndarray
        Array of transit durations in hours for each planet with transit data.
    Raises
    ------
    ValueError
        If identifier is neither a string nor numeric type.
    Notes
    -----
    - Only planets with tran_flag == 1 are included in the results.
    - For table='pscomppars', the function returns single values directly from the query.
    - For table='ps', the function returns the median period and duration, and the maximum
      transit mid-time from potentially multiple entries per planet.
    - All returned arrays have the same length, corresponding to the number of transiting
      planets found.
    """

    # Try to load from catalog first
    try:
        catalog = Catalog()
        catalog_name = 'Planet_Params.csv'
        
        # Determine hostname for catalog lookup
        if type(identifier) == float or type(identifier) == int:
            hostname = f"TIC {int(identifier)}"
        elif type(identifier) == str:
            hostname = string_checker(identifier)
        else:
            raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
        
        # Try to get data from catalog
        cached_data = catalog.get_host_data(catalog_name, hostname)
        
        if cached_data is not None and len(cached_data) > 0 and 'pl_orbper' in cached_data.columns:
            planets = cached_data['planet_letter'].unique().tolist()
            periods = cached_data.groupby('planet_letter')['pl_orbper'].first().values
            transit_mids = cached_data.groupby('planet_letter')['pl_tranmid'].first().values
            durations = cached_data.groupby('planet_letter')['pl_trandur'].first().values
            
            param_dict = {
                'periods': np.array(periods),
                'planets': planets,
                'transit_mids': np.array(transit_mids),
                'durations': np.array(durations)
            }
            return param_dict
    except (FileNotFoundError, KeyError, Exception) as e:
        print(f"ℹ️  Catalog not found or incomplete, querying astroquery... ({e})")

    if type(identifier) == float or type(identifier) == int:
        query = nea.query_criteria(table=table, select="pl_letter,tran_flag,pl_orbper,pl_tranmid,pl_trandur",
                                    where=f"tic_id='TIC {str(int(identifier))}'")
    elif type(identifier) == str:
        identifier = string_checker(identifier)
        query = nea.query_criteria(table=table, select="pl_letter,tran_flag,pl_orbper,pl_tranmid,pl_trandur",
                                    where=f"hostname like '{str(identifier)}'")
    else:
        raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
    planets = set(query['pl_letter'])
    periods = []
    transit_mids = []
    durations = []
    for planet in planets:
        tran_check = (query['tran_flag'] == 1) & (query['pl_letter'] == planet)
        if np.sum(tran_check) == 0:
            continue
        mask = (query['pl_letter'] == planet)
        query_planet = query[mask]
        if table == 'pscomppars':
            periods.append(float(query_planet['pl_orbper'].value)) 
            transit_mids.append(float(query_planet['pl_tranmid'].value))
            durations.append(float(query_planet['pl_trandur'].value))
        elif table == 'ps':
            periods.append(float(np.nanmedian(query_planet['pl_orbper'].value)))
            transit_mids.append(float(np.nanmax(query_planet['pl_tranmid'].value))) 
            durations.append(float(np.nanmedian(query_planet['pl_trandur'].value)))
    param_dict = {'periods': np.array(periods),'planets': list(planets), 'transit_mids': np.array(transit_mids), 'durations': np.array(durations)}
    if len(param_dict) ==  0:
        print(f"⚠️  No transit parameters found for {hostname}")
    # Save to catalog
    try:
        catalog = Catalog()
        catalog_name = 'Planet_Params.csv'
        catalog_path = os.path.join(catalog.catalog_dir, catalog_name)
        
        cache_data = []
        for idx, planet in enumerate(planets):
            row = {
                'hostname': hostname,
                'planet_letter': planet,
                'pl_orbper': periods[idx],
                'pl_tranmid': transit_mids[idx],
                'pl_trandur': durations[idx]
            }
            cache_data.append(row)
        
        cache_df = pd.DataFrame(cache_data)
        
        # Append or create catalog
        if os.path.exists(catalog_path):
            existing_df = pd.read_csv(catalog_path)
            # Remove old transit entries for this host
            existing_df = existing_df[~((existing_df['hostname'] == hostname) & (existing_df['planet_letter'].isin(planets)))]
            # Merge with existing data, preserving other columns
            for idx, new_row in cache_df.iterrows():
                matching_rows = existing_df[(existing_df['hostname'] == hostname) & 
                                           (existing_df['planet_letter'] == new_row['planet_letter'])]
                if len(matching_rows) > 0:
                    # Update existing row
                    for col in cache_df.columns:
                        existing_df.loc[matching_rows.index[0], col] = new_row[col]
                else:
                    # Add new row
                    existing_df = pd.concat([existing_df, pd.DataFrame([new_row])], ignore_index=True)
            cache_df = existing_df
        
        cache_df.to_csv(catalog_path, index=False)
    except Exception as e:
        print(f"⚠️  Failed to save to catalog. Error: {e}")
    
    return param_dict

def Query_Host_Params(identifier, table = "pscomppars", cache_file=None):
    """Query the NASA Exoplanet Archive for stellar parameters of a host star.
    This function retrieves stellar parameters (effective temperature, stellar radius, and stellar mass)
    for a star associated with a given identifier from the NASA Exoplanet Archive.
    
    Parameters
    ----------
    identifier : str, int, or float
        The stellar identifier. If numeric (int or float), it is treated as a TIC ID.
        If string, it is treated as a hostname and searched using a LIKE query.
    table : str, optional
        The NASA Exoplanet Archive table to query. Default is "pscomppars".
    cache_file : str, optional
        Path to a CSV file to cache the results. If provided, results will be saved to this file.
        Default is None (no caching).
    
    Returns
    -------
    Star : namedtuple
        A named tuple containing the following stellar parameters:
        - Teff : tuple of float
            Effective temperature of the star in Kelvin as (value, upper_error, lower_error).
        - Radius : tuple of float
            Stellar radius in solar radii as (value, upper_error, lower_error).
        - Mass : tuple of float
            Stellar mass in solar masses as (value, upper_error, lower_error).
    
    Notes
    -----
    - For table='ps', the function returns the median effective temperature and mass, and the maximum
      stellar radius from potentially multiple entries per star.
    """

    # Try to load from catalog first
    try:
        catalog = Catalog()
        catalog_name = 'Host_Params.csv'
        
        # Determine hostname for catalog lookup
        if type(identifier) == float or type(identifier) == int:
            hostname_lookup = f"TIC {int(identifier)}"
        elif type(identifier) == str:
            hostname_lookup = string_checker(identifier)
        else:
            raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
        
        # Try to get data from catalog
        cached_data = catalog.get_host_data(catalog_name, hostname_lookup)
        
        if cached_data is not None and len(cached_data) > 0:
            row = cached_data.iloc[0]
            Star = namedtuple("Star", ["Teff", "Radius", "Mass"])
            return Star(
                (row['st_teff'], row['st_tefferr1'], row['st_tefferr2']),
                (row['st_rad'], row['st_raderr1'], row['st_raderr2']),
                (row['st_mass'], row['st_masserr1'], row['st_masserr2'])
            )
    except (FileNotFoundError, KeyError, Exception) as e:
        print(f"ℹ️  Catalog not found or incomplete, querying astroquery... ({e})")

    # Determine input type and get hostname and TIC ID
    hostname = None
    tic_id = None
    
    if type(identifier) == float or type(identifier) == int:
        tic_id = int(identifier)
        try:
            host_query = nea.query_criteria(table=table, select="hostname",
                                            where=f"tic_id='TIC {str(tic_id)}'")
            if len(host_query) > 0:
                hostname = str(host_query['hostname'][0])
            else:
                hostname = f"TIC {tic_id}"
        except Exception as e:
            print(f"⚠️  Could not query hostname for TIC ID {tic_id}. Error: {e}")
            hostname = f"TIC {tic_id}"
        
        query = nea.query_criteria(table=table, select="st_teff,st_rad,st_mass,st_tefferr1,st_raderr1,st_masserr1,st_tefferr2,st_raderr2,st_masserr2",
                                    where=f"tic_id='TIC {str(int(identifier))}'")
    elif type(identifier) == str:
        identifier = string_checker(identifier)
        hostname = identifier
        
        # Query for TIC ID
        try:
            tic_query = nea.query_criteria(table=table, select="tic_id",
                                           where=f"hostname like '{str(identifier)}'")
            if len(tic_query) > 0:
                tic_str = str(tic_query['tic_id'][0])
                # Extract numeric part from "TIC XXXXXXX"
                if 'TIC' in tic_str:
                    tic_id = int(tic_str.split()[-1])
                else:
                    tic_id = None
            else:
                tic_id = None
        except Exception as e:
            print(f"⚠️  Could not query TIC ID for hostname {identifier}. Error: {e}")
            tic_id = None
        
        query = nea.query_criteria(table=table, select="st_teff,st_rad,st_mass,st_tefferr1,st_raderr1,st_masserr1,st_tefferr2,st_raderr2,st_masserr2",
                                    where=f"hostname like '{str(identifier)}'")
    else:
        raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
    
    if table == 'pscomppars':
        Teff = float(query['st_teff'][0].value)
        Teff_u = float(query['st_tefferr1'][0].value)
        Teff_l = float(query['st_tefferr2'][0].value)
        st_rad = float(query['st_rad'][0].value)
        st_rad_u = float(query['st_raderr1'][0].value)
        st_rad_l = float(query['st_raderr2'][0].value)
        st_mass = float(query['st_mass'][0].value)
        st_mass_u = float(query['st_masserr1'][0].value)
        st_mass_l = float(query['st_masserr2'][0].value)
    elif table == 'ps':
        Teff = float(np.nanmedian(query['st_teff'].value))
        Teff_u = float(np.nanmax(query['st_tefferr1'].value))
        Teff_l = float(np.nanmax(query['st_tefferr2'].value))
        st_rad = float(np.nanmax(query['st_rad'].value))
        st_rad_u = float(np.nanmax(query['st_raderr1'].value))
        st_rad_l = float(np.nanmax(query['st_raderr2'].value))
        st_mass = float(np.nanmedian(query['st_mass'].value))
        st_mass_u = float(np.nanmax(query['st_masserr1'].value))
        st_mass_l = float(np.nanmax(query['st_masserr2'].value))
    
    # Save to cache file if requested
    if cache_file is not None:
        try:
            cache_data = {
                'hostname': [hostname],
                'tic_id': [tic_id],
                'st_teff': [Teff],
                'st_tefferr1': [Teff_u],
                'st_tefferr2': [Teff_l],
                'st_rad': [st_rad],
                'st_raderr1': [st_rad_u],
                'st_raderr2': [st_rad_l],
                'st_mass': [st_mass],
                'st_masserr1': [st_mass_u],
                'st_masserr2': [st_mass_l]
            }
            cache_df = pd.DataFrame(cache_data)
            
            # Create directory if it doesn't exist
            cache_dir = os.path.dirname(cache_file)
            if cache_dir and not os.path.exists(cache_dir):
                os.makedirs(cache_dir)
            
            # Save to CSV
            cache_df.to_csv(cache_file, index=False)
        except Exception as e:
            print(f"⚠️  Failed to save cache file {cache_file}. Error: {e}")
    
    # Save to catalog
    try:
        catalog = Catalog()
        catalog_name = 'Host_Params.csv'
        cache_data = {
            'hostname': [hostname],
            'tic_id': [tic_id],
            'st_teff': [Teff],
            'st_tefferr1': [Teff_u],
            'st_tefferr2': [Teff_l],
            'st_rad': [st_rad],
            'st_raderr1': [st_rad_u],
            'st_raderr2': [st_rad_l],
            'st_mass': [st_mass],
            'st_masserr1': [st_mass_u],
            'st_masserr2': [st_mass_l]
        }
        cache_df = pd.DataFrame(cache_data)
        catalog_path = os.path.join(catalog.catalog_dir, catalog_name)
        
        # Append or create catalog
        if os.path.exists(catalog_path):
            existing_df = pd.read_csv(catalog_path)
            # Remove old entries for this host
            existing_df = existing_df[existing_df['hostname'] != hostname]
            cache_df = pd.concat([existing_df, cache_df], ignore_index=True)
        
        cache_df.to_csv(catalog_path, index=False)
    except Exception as e:
        print(f"⚠️  Failed to save to catalog. Error: {e}")
    
    Star = namedtuple("Star", ["Teff", "Radius", "Mass"])
    return Star((Teff,Teff_u,Teff_l), (st_rad,st_rad_u,st_rad_l), (st_mass,st_mass_u,st_mass_l))



def Query_Host_Params_TOI(identifier, cache_file=None):
    """Query the NASA Exoplanet Archive for stellar parameters of a host star using TOI identifier.

    This function retrieves stellar parameters (effective temperature and stellar radius with uncertainties)
    for a star associated with a given TOI (TESS Object of Interest) identifier from the NASA Exoplanet Archive.

    Parameters
    ----------
    identifier : int or str
        The TOI identifier number. Will be converted to integer format for the query.
    cache_file : str, optional
        Path to a CSV file to cache the results. If provided, results will be saved to this file.
        Default is None (no caching).

    Returns
    -------
        - Teff : tuple of float
            Effective temperature of the star in Kelvin as (value, upper_error, lower_error).
        - Radius : tuple of float
            Stellar radius in solar radii as (value, upper_error, lower_error).

    Raises
    ------
    ValueError
        If the identifier cannot be converted to an integer TOI ID.

    Notes
    -----
    - The function queries the NASA Exoplanet Archive using the 'toipfx' field.
    - Upper uncertainties (err1) are returned as positive values.
    - Lower uncertainties (err2) are returned as positive values (negation applied).
    - Only the first result from the query is returned.
    """
    
    # Try to load from catalog first
    toi_id = int(identifier)
    try:
        catalog = Catalog()
        catalog_name = 'Host_Params_TOI.csv'
        
        # Try to get data from catalog
        cached_data = catalog.query_catalog(catalog_name, toi_id=toi_id)
        
        if cached_data is not None and len(cached_data) > 0:
            print(f"✅ Found TOI host parameters in catalog for TOI {toi_id}")
            row = cached_data.iloc[0]
            Star = namedtuple("Star", ["Teff", "Radius"])
            return Star(
                (row['st_teff'], row['st_tefferr1'], row['st_tefferr2']),
                (row['st_rad'], row['st_raderr1'], row['st_raderr2'])
            )
    except (FileNotFoundError, KeyError, Exception) as e:
        print(f"ℹ️  Catalog not found or incomplete, querying astroquery... ({e})")
    
    try:
        query = nea.query_criteria(table="toi", select="st_teff,st_rad,st_tefferr1,st_raderr1,st_tefferr2,st_raderr2,hostname",
                                    where=f"toipfx='{str(int(identifier))}'")
    except:
        raise ValueError("Identifier must be the TOI ID as an integer.")
    
    toi_id = int(identifier)
    hostname = str(query['hostname'][0]) if 'hostname' in query.colnames and len(query) > 0 else f"TOI {toi_id}"
    
    Teff = float(query['st_teff'][0].value)
    Teff_u = float(query['st_tefferr1'][0].value)
    Teff_l = -float(query['st_tefferr2'][0].value)
    st_rad = float(query['st_rad'][0].value)
    st_rad_u = float(query['st_raderr1'][0].value)
    st_rad_l = -float(query['st_raderr2'][0].value)
    
    # Save to cache file if requested
    if cache_file is not None:
        try:
            cache_data = {
                'toi_id': [toi_id],
                'hostname': [hostname],
                'st_teff': [Teff],
                'st_tefferr1': [Teff_u],
                'st_tefferr2': [Teff_l],
                'st_rad': [st_rad],
                'st_raderr1': [st_rad_u],
                'st_raderr2': [st_rad_l]
            }
            cache_df = pd.DataFrame(cache_data)
            
            # Create directory if it doesn't exist
            cache_dir = os.path.dirname(cache_file)
            if cache_dir and not os.path.exists(cache_dir):
                os.makedirs(cache_dir)
            
            # Save to CSV
            cache_df.to_csv(cache_file, index=False)
            print(f"✅ Host parameters cached to: {cache_file}")
            print(f"   TOI ID: {toi_id}, Hostname: {hostname}")
        except Exception as e:
            print(f"⚠️  Failed to save cache file {cache_file}. Error: {e}")
    
    # Save to catalog
    try:
        catalog = Catalog()
        catalog_name = 'Host_Params_TOI.csv'
        cache_data = {
            'toi_id': [toi_id],
            'hostname': [hostname],
            'st_teff': [Teff],
            'st_tefferr1': [Teff_u],
            'st_tefferr2': [Teff_l],
            'st_rad': [st_rad],
            'st_raderr1': [st_rad_u],
            'st_raderr2': [st_rad_l]
        }
        cache_df = pd.DataFrame(cache_data)
        catalog_path = os.path.join(catalog.catalog_dir, catalog_name)
        
        # Append or create catalog
        if os.path.exists(catalog_path):
            existing_df = pd.read_csv(catalog_path)
            # Remove old entries for this TOI
            existing_df = existing_df[existing_df['toi_id'] != toi_id]
            cache_df = pd.concat([existing_df, cache_df], ignore_index=True)
        
        cache_df.to_csv(catalog_path, index=False)
    except Exception as e:
        print(f"⚠️  Failed to save to catalog. Error: {e}")
    
    Star = namedtuple("Star", ["Teff", "Radius"])
    return Star((Teff,Teff_u,Teff_l), (st_rad,st_rad_u,st_rad_l))

def Query_Transit_Solution_TOI(identifier, cache_file=None):
    """
    Query the NASA Exoplanet Archive for transit solution parameters of TOI planets.
    
    This function retrieves transit parameters (orbital period, transit mid-time, and transit duration)
    for planets associated with a given TOI (TESS Object of Interest) identifier from the NASA Exoplanet Archive.
    
    Parameters
    ----------
    identifier : int or float
        The TOI identifier number. Will be converted to integer format for the query.
    cache_file : str, optional
        Path to a CSV file to cache the results. If provided, results will be saved to this file.
        Default is None (no caching).
    
    Returns
    -------
    periods : numpy.ndarray
        Array of orbital periods in days for each planet with transit data.
    transit_mids : numpy.ndarray
        Array of transit mid-times (BJD) for each planet with transit data.
    durations : numpy.ndarray
        Array of transit durations in hours for each planet with transit data.
    
    Raises
    ------
    ValueError
        If identifier is not a numeric type (int or float).
    
    Notes
    -----
    - The function queries the 'toi' table using the 'toipfx' field.
    - The function returns the median period and duration, and the maximum
      transit mid-time from potentially multiple entries per TOI.
    - All returned arrays have the same length, corresponding to the number of
      TOIs found matching the identifier.
    """

    # Try to load from catalog first
    toi_id = int(identifier)
    try:
        catalog = Catalog()
        catalog_name = 'Transit_Solution_TOI.csv'
        
        # Try to get data from catalog
        cached_data = catalog.query_catalog(catalog_name, toi_id=toi_id)
        
        if cached_data is not None and len(cached_data) > 0:
            print(f"✅ Found TOI transit solution in catalog for TOI {toi_id}")
            planets = cached_data['toi'].unique().tolist()
            periods = cached_data.groupby('toi')['pl_orbper'].first().values
            transit_mids = cached_data.groupby('toi')['pl_tranmid'].first().values
            durations = cached_data.groupby('toi')['pl_trandurh'].first().values
            
            param_dict = {
                'periods': np.array(periods),
                'planets': planets,
                'transit_mids': np.array(transit_mids),
                'durations': np.array(durations)
            }
            return param_dict
    except (FileNotFoundError, KeyError, Exception) as e:
        print(f"ℹ️  Catalog not found or incomplete, querying astroquery... ({e})")

    if type(identifier) == float or type(identifier) == int:
        query = nea.query_criteria(table="toi", select="toi,pl_orbper,pl_tranmid,pl_trandurh,hostname",
                                    where=f"toipfx='{str(int(identifier))}'")
    else:
        raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
    
    toi_id = int(identifier)
    hostname = str(query['hostname'][0]) if 'hostname' in query.colnames and len(query) > 0 else f"TOI {toi_id}"
    
    planets = set(query['toi'])
    periods = []
    transit_mids = []
    durations = []
    for planet in planets:
        mask = (query['toi'] == planet)
        query_planet = query[mask]
        periods.append(float(np.nanmedian(query_planet['pl_orbper'].value)))
        transit_mids.append(float(np.nanmax(query_planet['pl_tranmid'].value))) 
        durations.append(float(np.nanmedian(query_planet['pl_trandurh'].value)))
    
    # Save to cache file if requested
    if cache_file is not None:
        try:
            cache_data = []
            for idx, planet in enumerate(planets):
                row = {
                    'toi_id': toi_id,
                    'hostname': hostname,
                    'toi': planet,
                    'pl_orbper': periods[idx],
                    'pl_tranmid': transit_mids[idx],
                    'pl_trandurh': durations[idx]
                }
                cache_data.append(row)
            
            cache_df = pd.DataFrame(cache_data)
            
            # Create directory if it doesn't exist
            cache_dir = os.path.dirname(cache_file)
            if cache_dir and not os.path.exists(cache_dir):
                os.makedirs(cache_dir)
            
            # Save to CSV
            cache_df.to_csv(cache_file, index=False)
            print(f"✅ Transit parameters cached to: {cache_file}")
            print(f"   TOI ID: {toi_id}, Hostname: {hostname}")
        except Exception as e:
            print(f"⚠️  Failed to save cache file {cache_file}. Error: {e}")
    
    # Save to catalog
    try:
        catalog = Catalog()
        catalog_name = 'TOI_Catalog.csv'
        cache_data = []
        for idx, planet in enumerate(planets):
            row = {
                'toi_id': toi_id,
                'hostname': hostname,
                'toi': planet,
                'pl_orbper': periods[idx],
                'pl_tranmid': transit_mids[idx],
                'pl_trandurh': durations[idx]
            }
            cache_data.append(row)
        
        cache_df = pd.DataFrame(cache_data)
        catalog_path = os.path.join(catalog.catalog_dir, catalog_name)
        
        # Append or create catalog
        if os.path.exists(catalog_path):
            existing_df = pd.read_csv(catalog_path)
            # Remove old entries for this TOI
            existing_df = existing_df[existing_df['toi_id'] != toi_id]
            cache_df = pd.concat([existing_df, cache_df], ignore_index=True)
        
        cache_df.to_csv(catalog_path, index=False)
        print(f"✅ TOI transit solution cached to catalog for TOI {toi_id}")
    except Exception as e:
        print(f"⚠️  Failed to save to catalog. Error: {e}")
    
    param_dict = {'periods': np.array(periods),'planets': list(planets), 'transit_mids': np.array(transit_mids), 'durations': np.array(durations)}
    return param_dict

def string_checker(identifier):
    """Function to standardize stellar identifiers by inserting spaces where appropriate."""
    if (identifier.startswith('BD')) and ' ' not in identifier:
        identifier = identifier[0:5] + ' ' + identifier[5:]
        return identifier
    if (identifier.startswith('Gliese')) and ' ' not in identifier:
        identifier = identifier[0:6] + ' ' + identifier[6:]
        return identifier
    if (identifier.startswith('LSPM')) and ' ' not in identifier:
        identifier = identifier[0:4] + ' ' + identifier[4:]
        return identifier
    if (identifier.startswith('L') and 'LHS' not in identifier and 'LTT' not in identifier and 'LP' not in identifier and 'Lupus' not in identifier) and ' ' not in identifier:
        identifier = identifier[0:1] + ' ' + identifier[1:]
        return identifier
    if (identifier.startswith('HD') or identifier.startswith('CD') or identifier.startswith('Gl') or identifier.startswith('WD') or identifier.startswith('HR') or identifier.startswith('GJ')) or identifier.startswith('LP') and ' ' not in identifier:
        identifier = identifier[0:2] + ' ' + identifier[2:]
        return identifier
    if (identifier.startswith('HIP') or identifier.startswith('LHS') or identifier.startswith('TAP') or (identifier.startswith('LTT')) or identifier.startswith('TIC')) and ' ' not in identifier:
        identifier = identifier[0:3] + ' ' + identifier[3:]
        return identifier
    if (identifier.startswith('EPIC') or identifier.startswith('Ross') or identifier.startswith('Wolf')) and ' ' not in identifier:
        identifier = identifier[0:4] + ' ' + identifier[4:]
        return identifier
    if (identifier[-1] == 'A' or identifier[-1] == 'B' or identifier[-1] == 'N'):
        identifier = identifier[0:-1] + ' ' + identifier[-1]
        return identifier
    if identifier.startswith('Proxima') and ' ' not in identifier:
        identifier = 'Proxima Cen'
    if identifier.startswith('AU') and ' ' not in identifier:
        identifier = 'AU Mic'
    if identifier.startswith('tau') and ' ' not in identifier:
        identifier = 'tau Boo'
    if identifier.startswith('tau') and ' ' not in identifier:
        identifier = 'tau Boo'
    if identifier.startswith('ups') and ' ' not in identifier:
        identifier = 'ups And'
    if identifier.startswith('HSPsc') and ' ' not in identifier:
        identifier = 'HS Psc'
    if identifier.startswith("Barnard") and ' ' not in identifier:
        identifier = "Barnard''s star"
    if identifier.startswith("51Peg") and ' ' not in identifier:
        identifier = "51 Peg"
        return identifier
    if identifier.startswith("61Vir") and ' ' not in identifier:
        identifier = "61 Vir"
        return identifier
    if identifier.startswith("Teegarden") and ' ' not in identifier:
        identifier = "Teegarden''s Star"
        return identifier
    if identifier.startswith("V830") and ' ' not in identifier:
        identifier = "V830 Tau"
        return identifier
    if identifier.startswith("V830") and ' ' not in identifier:
        identifier = "V830 Tau"
        return identifier
    if identifier.startswith("L1") and ' ' not in identifier:
        identifier = "L1 68"
        return identifier
    if identifier.startswith("YZ") and ' ' not in identifier:
        identifier = "YZ Cet"
        return identifier

def Query_Planet_Parameters(identifier, table = "ps", compute_epoch = False, cache_file = None):
    """
    Query the NASA Exoplanet Archive for periastron solution parameters of planets.
    This function retrieves periastron parameters (orbital period, periastron epoch, and argument of periastron)
    for planets associated with a given stellar identifier from the NASA Exoplanet Archive.
    Parameters
    ----------
    identifier : str, int, or float
        The stellar identifier. If numeric (int or float), it is treated as a TIC ID.
        If string, it is treated as a hostname and searched using a LIKE query.
    table : str, optional
        The NASA Exoplanet Archive table to query. Default is "ps".
        Supported values include "pscomppars" (Planetary Systems Composite Parameters)
        and "ps" (Planetary Systems). The aggregation method differs by table.
    compute_epoch : bool, optional
        If True, compute the periastron epoch from transit mid-time and argument of periastron, if available. Default is False.
    cache_file : str, optional
        Path to a CSV file to cache the results. If provided, results will be saved to this file.
        Default is None (no caching).
    Returns
    -------
    periods : numpy.ndarray
        Array of orbital periods in days for each planet with transit data.
    periastron_epochs : numpy.ndarray
        Array of periastron epochs (BJD) for each planet with transit data.
    omega_st : numpy.ndarray
        Array of omega_st for each planet with transit data.
    pl_letter : numpy.ndarray
        Array of planet letters for each planet with transit data.
    Raises
    ------
    ValueError
        If identifier is neither a string nor numeric type.
    Notes
    -----
    - For table='pscomppars', the function returns single values directly from the query.
    - For table='ps', the function returns the median period and duration, and the maximum
      transit mid-time from potentially multiple entries per planet.
    - All returned arrays have the same length, corresponding to the number of transiting
      planets found.
    """
    
    # Try to load from catalog first
    try:
        catalog = Catalog()
        catalog_name = 'Planet_Params.csv'
        
        # Determine hostname for catalog lookup
        if type(identifier) == float or type(identifier) == int:
            hostname_lookup = f"TIC {int(identifier)}"
        elif type(identifier) == str:
            hostname_lookup = string_checker(identifier)
        else:
            raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
        
        # Try to get data from catalog
        cached_data = catalog.get_host_data(catalog_name, hostname_lookup)
        
        if cached_data is not None and len(cached_data) > 0:
            param_dict = {}
            for _, row in cached_data.iterrows():
                planet = row['planet_letter']
                param_dict[planet] = {
                    'period': row['period'],
                    'peri_epoch': row['peri_epoch'],
                    'omega_st': row['omega_st'],
                    'transit_epoch': row['transit_epoch'],
                    'comp_peri_epoch': row['comp_peri_epoch'],
                    'a': row['a'],
                    'e': row['e']
                }
            return param_dict
    except (FileNotFoundError, KeyError, Exception) as e:
        print(f"ℹ️  Catalog not found or incomplete, querying astroquery... ({e})")
    
    # Determine input type and get hostname and TIC ID
    hostname = None
    tic_id = None
    
    if type(identifier) == float or type(identifier) == int:
        # Input is TIC ID, query for hostname
        tic_id = int(identifier)
        try:
            host_query = nea.query_criteria(table=table, select="hostname",
                                            where=f"tic_id='TIC {str(tic_id)}'")
            if len(host_query) > 0:
                hostname = str(host_query['hostname'][0])
            else:
                hostname = f"TIC {tic_id}"
        except Exception as e:
            print(f"⚠️  Could not query hostname for TIC ID {tic_id}. Error: {e}")
            hostname = f"TIC {tic_id}"
        
        if compute_epoch == False:
            query = nea.query_criteria(table=table, select="pl_letter,pl_orbper,pl_orbtper,pl_orblper,pl_orbeccen,pl_orbsmax,pl_tranmid",
                                        where=f"tic_id='TIC {str(int(identifier))}'")
        elif compute_epoch == True:
            query = nea.query_criteria(table=table, select="pl_letter,pl_orbper,pl_orbtper,pl_orblper,pl_orbeccen,pl_tranmid,pl_orbsmax,pl_tranmid",
                                        where=f"tic_id='TIC {str(int(identifier))}'")
    elif type(identifier) == str:
        identifier = string_checker(identifier)
        hostname = identifier
        
        # Query for TIC ID
        try:
            tic_query = nea.query_criteria(table=table, select="tic_id",
                                           where=f"hostname like '{str(identifier)}'")
            if len(tic_query) > 0:
                tic_str = str(tic_query['tic_id'][0])
                # Extract numeric part from "TIC XXXXXXX"
                if 'TIC' in tic_str:
                    tic_id = int(tic_str.split()[-1])
                else:
                    tic_id = None
            else:
                tic_id = None
        except Exception as e:
            print(f"⚠️  Could not query TIC ID for hostname {identifier}. Error: {e}")
            tic_id = None
        
        if compute_epoch == False:
            query = nea.query_criteria(table=table, select="pl_letter,pl_orbper,pl_orbtper,pl_orblper,pl_orbeccen,pl_orbsmax,pl_tranmid",
                                        where=f"hostname like '{str(identifier)}'")
        elif compute_epoch == True:
            query = nea.query_criteria(table=table, select="pl_letter,pl_orbper,pl_orbtper,pl_orblper,pl_orbeccen,pl_tranmid,pl_orbeccen,pl_orbsmax,pl_tranmid",
                                        where=f"hostname like '{str(identifier)}'")
    else:
        raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
    
    planets = set(query['pl_letter'])
    param_dict = {}
    for planet in planets:
        mask = (query['pl_letter'] == planet)
        query_planet = query[mask]
        if table == 'pscomppars':
            if compute_epoch == True:
                comp_peri_epoch = transit_to_periastron_epoch(float(query_planet['pl_orbper'].value), float(query_planet['pl_orbeccen'].value), float(query_planet['pl_orblper'].value))
            elif compute_epoch == False:
                comp_peri_epoch = np.nan
            param_dict[str(planet)] = {'period': float(query_planet['pl_orbper'].value), 'peri_epoch': float(query_planet['pl_orbtper'].value), 
                                       'omega_st': float(query_planet['pl_orblper'].value),'transit_epoch': float(np.nanmedian(query_planet['pl_tranmid'].value)), 'comp_peri_epoch': comp_peri_epoch, 'a': float(query_planet['pl_orbsmax'].value), 'e': float(query_planet['pl_orbeccen'].value)}
        elif table == 'ps':
            if compute_epoch == False:
                comp_peri_epoch = np.nan
            elif compute_epoch == True:
                comp_peri_epoch = transit_to_periastron_epoch(float(np.nanmedian(query_planet['pl_orbper'].value)), float(np.nanmedian(query_planet['pl_orbeccen'].value)), float(np.nanmedian(query_planet['pl_orblper'].value)))+ float(np.nanmax(query_planet['pl_tranmid'].value))
            param_dict[str(planet)] = {'period': float(np.nanmedian(query_planet['pl_orbper'].value)), 'peri_epoch': float(np.nanmax(query_planet['pl_orbtper'].value)), 
                                       'omega_st': float(np.nanmedian(query_planet['pl_orblper'].value)),  'transit_epoch': float(np.nanmedian(query_planet['pl_tranmid'].value)),'comp_peri_epoch': comp_peri_epoch, 'a': float(np.nanmedian(query_planet['pl_orbsmax'].value)), 'e': float(np.nanmedian(query_planet['pl_orbeccen'].value))}
    if len(param_dict) ==  0:
        print(f"⚠️  No transit parameters found for {hostname}")
    # Save to catalog
    try:
        catalog = Catalog()
        catalog_name = 'Planet_Params.csv'
        catalog_path = os.path.join(catalog.catalog_dir, catalog_name)
        
        cache_data = []
        for planet, params in param_dict.items():
            row = {'hostname': hostname, 'tic_id': tic_id, 'planet_letter': planet}
            row.update(params)
            cache_data.append(row)
        
        cache_df = pd.DataFrame(cache_data)
        
        # Append or create catalog
        if os.path.exists(catalog_path):
            existing_df = pd.read_csv(catalog_path)
            planets_to_update = list(param_dict.keys())
            # Remove old planet parameter entries for this host's planets
            existing_df = existing_df[~((existing_df['hostname'] == hostname) & 
                                       (existing_df['planet_letter'].isin(planets_to_update)))]
            # Merge with existing data, preserving other columns
            for idx, new_row in cache_df.iterrows():
                matching_rows = existing_df[(existing_df['hostname'] == hostname) & 
                                           (existing_df['planet_letter'] == new_row['planet_letter'])]
                if len(matching_rows) > 0:
                    # Update existing row
                    for col in cache_df.columns:
                        existing_df.loc[matching_rows.index[0], col] = new_row[col]
                else:
                    # Add new row
                    existing_df = pd.concat([existing_df, pd.DataFrame([new_row])], ignore_index=True)
            cache_df = existing_df
        
        cache_df.to_csv(catalog_path, index=False)
    except Exception as e:
        print(f"⚠️  Failed to save to catalog. Error: {e}")
    
    return param_dict

def MAST_data_extractor(input_dir):
    """Function to extract MAST data from nested directories and move them to a higher-level directory."""

    for root, dirs, files in os.walk(input_dir):
        # Get the top-level folder name
        rel_path = os.path.relpath(root, input_dir)
        if rel_path == '.':
            continue  # Skip the root directory itself
        
        top_folder = rel_path.split(os.sep)[0]
        print(top_folder)
        
        # Move all files from this level up to the top-level folder
        for file in files:
            src_file = os.path.join(root, file)
            dst_file = os.path.join(input_dir, top_folder, file)
            
            try:
                shutil.move(src_file, dst_file)
            except (NotADirectoryError, FileExistsError):
                continue
