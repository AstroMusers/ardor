<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 12:20:14 2025

@author: Nate Whitsett
"""

import numpy as np
import matplotlib.pyplot as plt
import ardor.Utils.Utils as UT
import pandas as pd
from scipy.stats import vonmises
from matplotlib import rcParams
from collections import namedtuple
from matplotlib import font_manager, patches

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
def orbit_pos(e,a,theta):
    return a*(1-e**2)/(1+e*np.cos(theta+np.pi))
font = FontProperties(fname=path)
italicfont = FontProperties(fname=italicpath)
rcParams["mathtext.fontset"] = "cm"
def a_from_period(R, period):
    return (((period*24*60*60)**2*6.67e-11*(1.989e30)*(R)**(5/4))/(4*np.pi**2))**(1/3)/(1.496e11)
def Polar_Flare_Plot(flare_phases, e = 0.01, a = 0.065, title = "Title", star_rad = 1, color_map = 'Reds', peri_err = None,
                     star_teff = 4000, transit = False, save_dir = None, alfven=[0.1, 0.025, 0.025], 
                     TOI = False, emax = 0, roche = [], dynamic_lim = False, legend = False):
    if alfven != None:
        alfven_guess = alfven[0]
        alfven_upper = alfven[1]
        alfven_lower = alfven[2]
    phases = np.linspace(0, 2*np.pi, num=10000)
    thetas = np.linspace(-np.pi, np.pi, num= 10000)
    ### STAR
    if star_teff > 3000 and star_teff < 3900:
        print(1)
        star_cl = 'orangered'
    if star_teff >= 3900 and star_teff < 5200:
        star_cl = 'goldenrod'
    if star_teff >= 5200 and star_teff < 6000:
        star_cl = 'yellow'
    if star_teff >= 6000 and star_teff < 7600:
        star_cl = 'palegoldenrod'
    elif star_teff > 7600:
        star_cl = 'paleturquoise'
    star_rad = star_rad*0.00465
    
    
    # Generate random data for angles and radii
    if peri_err != None:
        if len(peri_err) == 1:
            peri_err[0] = UT.range_shift(peri_err[0], 0, 1, 0, 2*np.pi)
        if len(peri_err) == 2:
            peri_err[0] = UT.range_shift(peri_err[0], 0, 360, 0, 2*np.pi)
            peri_err[1] = UT.range_shift(peri_err[1], 0, 360, 0, 2*np.pi)
    print(e, a, len(phases))
    orbit = orbit_pos(e,a,phases)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize = (3.5,3.5))
    ax.axis('off')
    print(flare_phases)
    for index, phase in enumerate(flare_phases):
        adj_phase = UT.range_shift(phase, 0, 1, 0, 2*np.pi)
        if index == 0:
            if len(flare_phases) > 30:
                ax.vlines(adj_phase, ymin=0, ymax=orbit_pos(e,a,adj_phase), linewidth=0.75, color='purple', linestyle='--', label='Flares')
            else:
                ax.vlines(adj_phase, ymin=0, ymax=orbit_pos(e,a,adj_phase), linewidth=1.25, color='purple', linestyle='--', label='Flares')
        else:
            if len(flare_phases) > 30:
                ax.vlines(adj_phase, ymin=0, ymax=orbit_pos(e,a,adj_phase), linewidth=0.75, color='purple', linestyle='--')
            else:
                ax.vlines(adj_phase, ymin=0, ymax=orbit_pos(e,a,adj_phase), linewidth=1.25, color='purple', linestyle='--')
    if TOI == False:
        ax.plot(phases, orbit, linestyle='dashdot', color='darkblue', linewidth=1.25, label = 'Orbit')
    ax.set_rlim(0,0.15)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_xticklabels([], font=font, zorder=0)
    if transit == False:
        # ax.text(np.pi+0.08, orbit_pos(e,a,np.pi)+a*0.85, 'Periastron\n  Phase', font=font, fontsize=9.5)
        ax.vlines(np.pi, ymin=0, ymax=orbit_pos(e,a,np.pi), zorder=3,color='black',linestyle='dotted', label = 'Periastron')
    if transit == True:
        # ax.text(np.pi+.08, a*1.62, 'Mid-Transit', font=font, fontsize=10)
        ax.vlines(np.pi, ymin=0, ymax=orbit_pos(e,a,np.pi), zorder=3,color='black',linestyle='dotted', label = 'Mid-Transit')
    ax.set_title(title, font=font, fontsize=14)
    loc, kappa, TS = MLE.VM_Unbinned_likelihood(flare_phases)
    loc = UT.range_shift(loc, 0,1,-np.pi, np.pi)
    PDF = vonmises.pdf(thetas, loc=loc+np.pi, kappa=kappa)
    ax.plot(thetas, (PDF/np.max(PDF))*orbit_pos(e,a,loc+np.pi), linestyle = '-', color='cyan', label = 'VM PDF', linewidth=1)
    # plt.text(0.24, a*1.35, "AU", font = font)
    # plt.text(0.375, 0.0, "a = " + str(round(a,3)) + " AU",font = font, fontsize=11,transform=ax.transAxes)
    # plt.text(0.375, 0.0, "a = " + str(round(a,3)) + " AU",font = font, fontsize=11,transform=ax.transAxes)
    print(loc)
    plt.fill(thetas, star_rad*np.ones(len(thetas)), zorder=4, c=star_cl)
    if alfven != None:
        plt.fill(thetas, (alfven_guess + alfven_upper)*np.ones(len(thetas)), c = 'red',alpha = 0.1)
        plt.fill(thetas, (alfven_guess - np.abs(alfven_lower))*np.ones(len(thetas)), c = 'red', alpha = 0.2)
        plt.fill(thetas, alfven_guess *np.ones(len(thetas)), c = 'red',alpha = 0.2)
        plt.plot(thetas, alfven_guess*np.ones(len(thetas)), c = 'red', linewidth = 0.75, label = 'Alfvén Surface')
    if TOI == True:
        for samples in range(500):
            e_samp = np.random.uniform(0, emax)
            samp_orbit = orbit_pos(e_samp, a, thetas+np.random.uniform(0,360))
            if samples == 0:
                plt.plot(thetas, samp_orbit, color='black', alpha=0.02, zorder=-1)
                plt.plot(0, 0, color='black', label = 'Sampled Orbits', zorder=-1, linewidth=0.75)
            else:
                plt.plot(thetas, samp_orbit, color='black', alpha=0.02,zorder=-1)
        if dynamic_lim == True:
            plt.plot(thetas, roche*np.ones(len(thetas)), c = 'blue', linewidth = 0.75,linestyle='--', label = 'Dyn. Lim.')
        else:
            plt.plot(thetas, roche*np.ones(len(thetas)), c = 'green', linewidth = 0.75,linestyle='--', label = 'Roche Lim.')
    if peri_err != None and transit == False:
        if len(peri_err) == 1:
            ax.fill_between(np.linspace(np.pi-peri_err[0], np.pi+peri_err[0]), 0 , orbit_pos(e,a,np.pi), color='black', alpha=0.5,zorder=0, label='$\delta$ Periastron')
        elif len(peri_err) == 2:
            ax.fill_between(np.linspace(np.pi-np.abs(peri_err[0]), np.pi+np.abs(peri_err[1])), 0 , orbit_pos(e,a,np.pi), color='black', alpha=0.5,zorder=0, label='$\delta$ Periastron')
    if legend == True:
        ax.legend(handlelength = 1.5, columnspacing=1, labelspacing = .2, prop=dict(family='Serif', size=9), ncol=2, loc = 'lower center')
    if save_dir != None:
        plt.savefig(save_dir + '/' + title + '_Polar_Hist.png', dpi=500, bbox_inches='tight')
        plt.show()



flare_dir = "C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_News.csv"
planet_params = "C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOI_List.csv"
# peri_params_dir = "C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Arg_Peri_List_New.csv"
flares = pd.read_csv(flare_dir)
params = pd.read_csv(planet_params)
# peri_params = pd.read_csv(peri_params_dir)
# (flare_phases, e = 0.01, a = 0.065, bins = 10, title = "Title", star_rad = 1, color_map = 'Reds', peri_err = None, star_type = "G", transit = False, save_dir = None):

############### HOSTS###################### 
### TRANSIT
# target_hosts = ['TOI-1807','TOI-1518','TOI-1452','TOI-1288','KOI-7368', 'TOI-540']
# title_hosts = ['TOI-1807','TOI-1518','TOI-1452','TOI-1288','KOI-7368', 'TOI-540']
### PERI
# target_hosts = ['TOI-1062','AUMic', 'WASP-34', 'Kepler-1972', 'HD260655', 'Gl49', 'TOI-2145']
# title_hosts = ['TOI-1062','AU Mic', 'WASP-34', 'Kepler-1972', 'HD 260655', 'Gl 49', 'TOI-2145']
letters = ['b','b','b','b','b','b', 'b']

### TOI
target_hosts = [6276.01, 2329.01,1232.01, 723.01, 1520.01,1887.01, 5051.01, 1187.01]
title_hosts = ['6276.01', '2329.01','1232.01', '723.01', '1520.01','1887.01', '5051.01', '1187.01']
# letters = ['b','b','b','b','b','b', 'b','b']
emax_list = [0.2458, 0.6459,0.82, 0.5, 0.71,0.699, 0.3, 0.626]
a_list = [0.048, 0.0658,0.117, 0.0229, 0.0693, 0.0287, 0.0418, 0.041]
roche_list = [0.0598, 0.023,0.0214, 0.0114, 0.01968,0.01, 0.029, 0.01533]
for index, targets in enumerate(target_hosts[0:2]):
    
    ###SYSTEM PARAMS Alfven_Rad,0.13,-0.135
    # alfven_guess = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'Alfven_Rad'].item()
    # alfven_upper = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), '0.13'].item()
    # alfven_lower = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), '-0.135'].item()
    # alfven = [alfven_guess, alfven_upper, alfven_lower]
    period = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'pl_orbper'].item()
    # e = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'pl_orbeccen'].item()
    # a = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'pl_orbsmax'].item()
    # peri_err_l = peri_params.loc[(peri_params['Host_ID'] == targets) & (peri_params['pl_letter'] == letters[index]), 'pl_orblpererr1'].item()
    # peri_err_u = peri_params.loc[(peri_params['Host_ID'] == targets) & (peri_params['pl_letter'] == letters[index]), 'pl_orblpererr2'].item()
    ###STELLAR PARAMS
    star_rad = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'st_rad']
    st_teff = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'st_teff']
    ###FLARE PHASES
    epochs = np.array(flares.loc[flares['Host_ID'] == targets, 'Flare_Epoch'])
    tran_mid = params.loc[(params['Host_ID'] == targets) & (params['pl_letter'] == letters[index]), 'pl_tranmid'].item()
    # if np.isnan(tran_mid) == True:
    #     tran_mid = peri_params.loc[(peri_params['Host_ID'] == str(targets)) & (peri_params['pl_letter'] == letters[index]), 'Comp_Epoch_Peri'].item()
    phases = ((epochs - (float(tran_mid) + float(period))) % float(period))/float(period)
    if index == -1:
        continue
        # Polar_Flare_Plot(np.array(phases),title=title_hosts[index] + ' ' + letters[index], e=e,a=a, peri_err = [peri_err_l,peri_err_u],
        #                   star_rad=float(star_rad), star_teff=float(st_teff), transit=False, 
        #                   save_dir="C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Publication_2024/Polar_Hists"
        #                   , alfven=alfven, legend=True)
    else:
        Polar_Flare_Plot(np.array(phases),title='TOI-' + title_hosts[index], e=0,a=a_list[index],peri_err = None,
                          star_rad=float(star_rad), star_teff=float(st_teff), transit=True, roche=roche_list[index], TOI = True, emax = emax_list[index], 
                          save_dir="C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Publication_2024/Polar_Hists"
                          , legend=False, alfven=None)
    plt.show()
    plt.clf()
=======
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

>>>>>>> 843bcfd59f29e0297b65b05360ddcf3ed4f7ad77
