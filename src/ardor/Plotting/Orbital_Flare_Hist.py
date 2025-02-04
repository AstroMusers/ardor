# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 12:20:14 2025

@author: whitsett.n
"""

import numpy as np
import matplotlib.pyplot as plt
import ardor.Utils.Utils as UT
import pandas as pd
from matplotlib import rcParams
from collections import namedtuple
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
import ardor.Statistical_Tests.MLE as MLE
for fontpath in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    if 'lmroman10-regular'.lower() in fontpath.lower():
        path = fontpath
        print(path)
for fontpath in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    if 'lmroman10-italic'.lower() in fontpath.lower():
        italicpath = fontpath
        print(italicpath)

font = FontProperties(fname=path)
italicfont = FontProperties(fname=italicpath)
rcParams["mathtext.fontset"] = "cm"
data_dir = "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares_New.csv"
# data_dir = "C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/src/ardor/SPI_Forward_Models/Sim_Outputs_Kappa/M_Type/Kappa_Sim_M.csv"
data = pd.read_csv(data_dir)
flares = data.loc[data['Host_ID'] == "TOI-833", "Transit_Phase"]
def Polar_Flare_Plot(flare_phases, e = 0.01, a = 0.065, bins = 10, title = "Title", star_rad = 1, color_map = 'Reds', peri_err = None, star_type = "G", transit = False, save_dir = None):
    star_rad = star_rad*0.00465
    # Generate random data for angles and radii
    flares = UT.range_shift(flare_phases, 0, 1, 0, 2*np.pi)
    if peri_err != None:
        if len(peri_err) == 1:
            peri_err[0] = UT.range_shift(peri_err[0], 0, 1, 0, 2*np.pi)
        if len(peri_err) == 2:
            peri_err[0] = UT.range_shift(peri_err[0], 0, 1, 0, 2*np.pi)
            peri_err[1] = UT.range_shift(peri_err[1], 0, 1, 0, 2*np.pi)
    thetas, PDF = MLE.VM_Unbinned_likelihood(flares)
    phases = np.linspace(0, 2*np.pi, num=10000)
    theta = np.linspace(0, 2 * np.pi - 2*np.pi/bins, num= bins)
    orbit = a*(1-e**2)/(1+e*np.cos(phases))
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = (3.5,3.5))
    # Create the histogram
    n, bins, patches = ax.hist(flares, bins=bins,  density=True, alpha=0.75, range=(0,2*np.pi))
    plt.clf()
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize = (3.5,3.5))
    for index, b in enumerate(n):
        n[index] = (n[index]/np.max(n))*a
    cmap = plt.get_cmap(color_map)
    norm = plt.Normalize(min(n), max(n))
    ax.bar(theta, n, width=2*np.pi/len(theta), edgecolor='black',zorder=2, alpha=0.75, color=cmap(norm(n)), label = 'Binned Phases')
    ax.plot(phases, orbit, linestyle='dashdot', color='darkblue', linewidth=1.25, label = 'Orbit')
    ax.set_rlim(0,a*1.64)
    rticks = np.linspace(a*1.3, 0, num=3, endpoint=False)[::-1]
    ax.set_yticks(np.round(rticks,3))
    ax.set_yticklabels(np.round(rticks,3), font=font, fontsize=9)
    ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 5*np.pi/4, 3*np.pi/2, 7*np.pi/4])
    ax.set_xticklabels(['', '', '', '', '', '', '', ''], font=font, zorder=0)
    if transit == False:
        ax.text(np.pi+0.08, a*1.62, 'Periastron\n  Phase', font=font, fontsize=7.5)
        ax.vlines(np.pi, ymin=0, ymax=a, zorder=3,color='black',linestyle='dotted', label = 'Periastron')
    if transit == True:
        ax.text(np.pi+.08, a*1.62, 'Transit\n Phase', font=font, fontsize=9)
        ax.vlines(np.pi, ymin=0, ymax=a, zorder=3,color='black',linestyle='dotted', label = 'Transit')
    # ax.set_title(title, font=font, fontsize=12)
    if peri_err != None and transit == False:
        if len(peri_err) == 1:
            ax.fill_between(np.linspace(np.pi-peri_err[0], np.pi+peri_err[0]), 0 , a, color='black', alpha=0.5,zorder=0, label='$\delta$ Periastron')
        elif len(peri_err) == 2:
            ax.fill_between(np.linspace(np.pi-peri_err[0], np.pi+peri_err[1]), 0 , a, color='black', alpha=0.5,zorder=0, label='$\delta$ Periastron')
    plt.fill(theta, star_rad*np.ones(len(theta)), zorder=4, color='red')
    ax.plot(thetas, (PDF/max(PDF))*(a), linestyle = 'dotted', color='red', label = 'VM-PDF')

    ax.legend(handlelength = 1.5, columnspacing=1, labelspacing = .2, prop=dict(family='Serif', size=7), ncol=2, loc='upper center')
    plt.text(0.13, 0.006, "AU", font = font)
    if save_dir != None:
        plt.savefig(save_dir + '/' + title + '_Polar_Hist.png', dpi=500, bbox_inches='tight')
        plt.show()

Polar_Flare_Plot(flares, 0, 0.017, bins=10, title="TOI-833", peri_err = None, star_rad = 0.6, transit = True, save_dir = 'C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Publication_2024')

