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
            print(f"✅ Host parameters cached to: {cache_file}")
            print(f"   Hostname: {hostname}, TIC ID: {tic_id}")
        except Exception as e:
            print(f"⚠️  Failed to save cache file {cache_file}. Error: {e}")
    
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
    
    param_dict = {'periods': np.array(periods),'planets': list(planets), 'transit_mids': np.array(transit_mids), 'durations': np.array(durations)}
    return param_dict

def string_checker(identifier):
    """Function to standardize stellar identifiers by inserting spaces where appropriate."""
    if (identifier.startswith('L') and 'LHS' not in identifier and 'LP' not in identifier and 'Lupus' not in identifier) and ' ' not in identifier:
        identifier = identifier[0:2] + ' ' + identifier[2:]
    if (identifier.startswith('HD') or identifier.startswith('CD') or identifier.startswith('Gl') or identifier.startswith('WD') or identifier.startswith('HR') or identifier.startswith('GJ')) and ' ' not in identifier:
        identifier = identifier[0:2] + ' ' + identifier[2:]
    if (identifier.startswith('HIP') or identifier.startswith('LHS') or identifier.startswith('TAP') or identifier.startswith('TIC')) and ' ' not in identifier:
        identifier = identifier[0:3] + ' ' + identifier[3:]
    if (identifier.startswith('EPIC') or identifier.startswith('Ross') or identifier.startswith('Wolf')) and ' ' not in identifier:
        identifier = identifier[0:4] + ' ' + identifier[4:]
    if (identifier.startswith('Gliese')) and ' ' not in identifier:
        identifier = identifier[0:7] + ' ' + identifier[7:]
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
    if identifier.startswith("Barnard'sstar") and ' ' not in identifier:
        identifier = "Barnard's star"
    if identifier.startswith("51Peg") and ' ' not in identifier:
        identifier = "51 Peg"
    if identifier.startswith("61Vir") and ' ' not in identifier:
        identifier = "61 Vir"
    if identifier.startswith("Teegarden") and ' ' not in identifier:
        identifier = "Teegarden's Star"
    if identifier.startswith("V830") and ' ' not in identifier:
        identifier = "V830 Tau"
    if identifier.startswith("V830") and ' ' not in identifier:
        identifier = "V830 Tau"
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
    comp_peri_epochs = []
    param_dict = {}
    for planet in planets:
        mask = (query['pl_letter'] == planet)
        query_planet = query[mask]
        if table == 'pscomppars':
            if compute_epoch == True:
                comp_peri_epochs.append(transit_to_periastron_epoch(float(query_planet['pl_orbper'].value), float(query_planet['pl_orbeccen'].value), float(query_planet['pl_orblper'].value)))
            elif compute_epoch == False:
                comp_peri_epochs.append(np.nan)
            param_dict[str(planet)] = {'period': float(query_planet['pl_orbper'].value), 'peri_epoch': float(query_planet['pl_orbtper'].value), 
                                       'omega_st': float(query_planet['pl_orblper'].value),'transit_epoch': float(np.nanmedian(query_planet['pl_tranmid'].value)), 'comp_peri_epoch': comp_peri_epochs, 'a': float(query_planet['pl_orbsmax'].value), 'e': float(query_planet['pl_orbeccen'].value)}
        elif table == 'ps':
            if compute_epoch == False:
                comp_peri_epochs = np.nan
            elif compute_epoch == True:
                comp_peri_epochs = transit_to_periastron_epoch(float(np.nanmedian(query_planet['pl_orbper'].value)), float(np.nanmedian(query_planet['pl_orbeccen'].value)), float(np.nanmedian(query_planet['pl_orblper'].value)))+ float(np.nanmax(query_planet['pl_tranmid'].value))
            param_dict[str(planet)] = {'period': float(np.nanmedian(query_planet['pl_orbper'].value)), 'peri_epoch': float(np.nanmax(query_planet['pl_orbtper'].value)), 
                                       'omega_st': float(np.nanmedian(query_planet['pl_orblper'].value)),  'transit_epoch': float(np.nanmedian(query_planet['pl_tranmid'].value)),'comp_peri_epoch': comp_peri_epochs, 'a': float(np.nanmedian(query_planet['pl_orbsmax'].value)), 'e': float(np.nanmedian(query_planet['pl_orbeccen'].value))}
    
    # Save to cache file if requested
    if cache_file is not None:
        try:
            # Convert param_dict to DataFrame with hostname and TIC ID
            cache_data = []
            for planet, params in param_dict.items():
                row = {'hostname': hostname, 'tic_id': tic_id, 'planet_letter': planet}
                row.update(params)
                cache_data.append(row)
            
            cache_df = pd.DataFrame(cache_data)
            
            # Create directory if it doesn't exist
            cache_dir = os.path.dirname(cache_file)
            if cache_dir and not os.path.exists(cache_dir):
                os.makedirs(cache_dir)
            
            # Save to CSV
            cache_df.to_csv(cache_file, index=False)
            print(f"✅ Planet parameters cached to: {cache_file}")
            print(f"   Hostname: {hostname}, TIC ID: {tic_id}")
        except Exception as e:
            print(f"⚠️  Failed to save cache file {cache_file}. Error: {e}")
    
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


Query_Host_Params('AU Mic', cache_file='/ugrad/whitsett.n/ardor/test_cache/Host_Params.csv')