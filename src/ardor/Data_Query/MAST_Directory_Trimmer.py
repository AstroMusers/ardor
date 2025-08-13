# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:30:47 2023

@author: Nate Whitsett
"""

import os
import shutil


def MAST_data_extractor(input_dir):
    ##Takes a list of output folders from MAST
    ##and consolidates the data files into a
    ##single directory
        target_folders = os.listdir(input_dir)
        for folders in target_folders:
            sub_target_folders = os.listdir(input_dir + '/' + folders)
            print(folders)
            for project in sub_target_folders:
                # try:
                #     if project == 'mastDownload':
                #         shutil.rmtree(input_dir + '/' + folders + '/mastDownload')
                # except:
                #     continue
                try:
                    sub_sub_target_folders = os.listdir(input_dir + '/' + folders + '/' + project)
                    for data_folders in sub_sub_target_folders:
                        sub_sub_sub_target_folders = os.listdir(input_dir + '/' + folders + '/' + project + '/' + data_folders)
                        for data in sub_sub_sub_target_folders:
                            sub_sub_sub_sub_target_folders = os.listdir(input_dir + '/' + folders + '/' + project + '/' + data_folders + '/' + data)
                            for data_files in sub_sub_sub_sub_target_folders:
                                shutil.move(input_dir + '/' + folders + '/' + project + '/' + data_folders + '/' + data + '/' + data_files, input_dir + '/' + folders)
                except NotADirectoryError:
                    continue
        
                    
                
MAST_data_extractor("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/TESS Data/All_Exoplanet_Hosts/")
                