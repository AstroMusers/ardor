# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:43:03 2023

@author: Nate Whitsett
"""

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import random as rd

def CVZ(RA, DEC, radius=8.05, resolution=100):
    deg = np.linspace(0, 2*np.pi, num=resolution)
    coords = []
    RA_list = []
    DEC_list = []
    pole = False
    r = 0
    if DEC + radius > 90 or DEC - radius < -90:
        pole = True
    for theta in deg:
        new_RA = RA + np.cos(theta)*radius 
        new_DEC = DEC + np.sin(theta)*radius
        if pole == False:
            if new_RA > 360:
                new_RA -= 360
            elif new_RA < 0:
                new_RA += 360
        elif pole == True:
            new_RA = (theta/(2*np.pi))*360
            r = (90-np.abs(DEC))*np.cos(np.pi*new_RA/180 - np.pi*RA/180) + np.sqrt(radius**2 - (90-np.abs(DEC))**2*(np.sin(np.pi*new_RA/180 - np.pi*RA/180))**2)
            if DEC < 0:
                new_DEC = -90 + r
            if DEC > 0:
                new_DEC = 90 - r
        coords.append((new_RA, new_DEC))
        RA_list.append(new_RA)
        DEC_list.append(new_DEC)
        
    return coords, np.array(RA_list), np.array(DEC_list)

def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = np.remainder(RA+360-org,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection)
    ax.scatter(np.radians(x),np.radians(Dec), s=10)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)

def cost_function(CVZ_RA, CVZ_DEC, CVZ_radius, targets_RAs, targets_DECs, target_weight=1):
    coord = SkyCoord(CVZ_RA, CVZ_DEC, unit='deg')
    CVZ_RA = float(coord.ra.degree)
    CVZ_DEC = float(coord.dec.degree)
    index = 0 
    target_count = 0
    for index in range(len(targets_RAs)):
        target_RA = targets_RAs[index] - CVZ_RA
        target_DEC = targets_DECs[index] - CVZ_DEC
        dist = np.sqrt(target_RA**2 + target_DEC**2)
        if dist < CVZ_radius:
            target_count += 1*target_weight
    return target_count

def cost_function_map(targets_RAs, targets_DECs, output_dir=None, resolution=50):
    CF_map = []
    RA_span = np.linspace(0, 360, num=resolution)
    DEC_span = np.linspace(-90, 90, num=resolution)
    for i in range(len(DEC_span)):
        values_DEC = []
        for j in range(len(RA_span)):
            a = cost_function(RA_span[j], DEC_span[i], 8, targets_RAs, targets_DECs)
            values_DEC.append(a)
        CF_map.append(values_DEC)
    plt.imshow(CF_map)
    plt.xscale()
    if output_dir != None:
        np.savetxt(output_dir, CF_map, delimiter=',')

TOIs = pd.read_csv('/data/whitsett.n/ardor/templates/csv-file-toi-catalog.csv')
TOI_RA = np.array(TOIs['TIC Right Ascension'])
TOI_DEC = np.array(TOIs['TIC Declination'])

cost_function_map(TOI_RA, TOI_DEC)
# values_RA = []
# values_DEC = []
# RA_span = np.linspace(0, 360)
# DEC_span = np.linspace(-90, 90)
# plt.scatter(TOI_RA, TOI_DEC)
# for i in range(len(DEC_span)):
#     values_DEC = []
#     for j in range(len(RA_span)):
#         a = cost_function(RA_span[j], DEC_span[i], 8, TOI_RA, TOI_DEC)
#         values_DEC.append(a)
#     values_RA.append(values_DEC)
# np.savetxt('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/cost_func.csv', values_RA, delimiter=',')