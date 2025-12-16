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
from collections import namedtuple
from astroquery.mast import Catalogs
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
    return np.array(periods), np.array(transit_mids), np.array(durations)

def Query_Host_Params(identifier, table = "pscomppars"):
    """Query the NASA Exoplanet Archive for stellar parameters of a host star.
    This function retrieves stellar parameters (effective temperature, stellar radius, and stellar mass)
    for a star associated with a given identifier from the NASA Exoplanet Archive.
    Star : namedtuple
        A named tuple containing the following stellar parameters:
        - Teff : float
            Effective temperature of the star in Kelvin.
        - st_rad : float
            Stellar radius in solar radii.
        - st_mass : float
            Stellar mass in solar masses.
    - For table='ps', the function returns the median effective temperature and mass, and the maximum
      stellar radius from potentially multiple entries per star.
    """

    if type(identifier) == float or type(identifier) == int:
        query = nea.query_criteria(table=table, select="st_teff,st_rad,st_mass,st_tefferr1,st_raderr1,st_masserr1,st_tefferr2,st_raderr2,st_masserr2",
                                    where=f"tic_id='TIC {str(int(identifier))}'")
    elif type(identifier) == str:
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
    Star = namedtuple("Star", ["Teff", "Radius", "Mass"])
    return Star((Teff,Teff_u,Teff_l), (st_rad,st_rad_u,st_rad_l), (st_mass,st_mass_u,st_mass_l))



def Query_Host_Params_TOI(identifier):
    """Query the NASA Exoplanet Archive for stellar parameters of a host star using TOI identifier.

    This function retrieves stellar parameters (effective temperature and stellar radius with uncertainties)
    for a star associated with a given TOI (TESS Object of Interest) identifier from the NASA Exoplanet Archive.

    Parameters
    ----------
    identifier : int or str
        The TOI identifier number. Will be converted to integer format for the query.

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
        query = nea.query_criteria(table=table, select="st_teff,st_rad,st_tefferr1,st_raderr1,st_tefferr2,st_raderr2",
                                    where=f"toipfx='{str(int(identifier))}'")
    except:
        raise ValueError("Identifier must be the TOI ID as an integer.")
    print(query)
    Teff = float(query['st_teff'][0].value)
    Teff_u = float(query['st_tefferr1'][0].value)
    Teff_l = -float(query['st_tefferr2'][0].value)
    st_rad = float(query['st_rad'][0].value)
    st_rad_u = float(query['st_raderr1'][0].value)
    st_rad_l = -float(query['st_raderr2'][0].value)
    Star = namedtuple("Star", ["Teff", "Radius"])
    return Star((Teff,Teff_u,Teff_l), (st_rad,st_rad_u,st_rad_l))

def Query_Transit_Solution_TOI(identifier):
    """
    Query the NASA Exoplanet Archive for transit solution parameters of TOI planets.
    
    This function retrieves transit parameters (orbital period, transit mid-time, and transit duration)
    for planets associated with a given TOI (TESS Object of Interest) identifier from the NASA Exoplanet Archive.
    
    Parameters
    ----------
    identifier : int or float
        The TOI identifier number. Will be converted to integer format for the query.
    
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
        query = nea.query_criteria(table="toi", select="toi,pl_orbper,pl_tranmid,pl_trandurh",
                                    where=f"toipfx='{str(int(identifier))}'")
    else:
        raise ValueError("Identifier must be a string (hostname) or numeric (TIC ID).")
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
    return np.array(periods), np.array(transit_mids), np.array(durations)