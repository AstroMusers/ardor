# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:06:28 2023

@author: Nathan
"""

from astroquery.mast import Observations
import os
import shutil
import pandas as pd
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
