# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 12:20:14 2025

@author: Nate Whitsett
"""

import numpy as np
import matplotlib.pyplot as plt
import ardor.Utils.Utils as UT
from scipy.stats import vonmises
from matplotlib import rcParams
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
import ardor.Statistical_Tests.MLE as MLE
from matplotlib.animation import FuncAnimation, PillowWriter
for fontpath in font_manager.findSystemFonts():
    if 'lmroman10-regular' in fontpath.lower():
        path = fontpath
    if 'lmroman10-italic' in fontpath.lower():
        italicpath = fontpath

# Register and get name
font_manager.fontManager.addfont(path)
font_manager.fontManager.addfont(italicpath)

font = FontProperties(fname=path)
font_name = font.get_name()

rcParams["font.family"] = font_name
rcParams["mathtext.fontset"] = "cm"

def St_Polar_Flare_Plot(
    flare_phases,
    a=0.02,
    title="Title",
    star_rad=1,
    color_map='Reds',
    star_teff=4000,
    save_dir=None,
    plot_stats=False,
    bottom_text_list=None
):
    # Define thetas for PDF
    thetas = np.linspace(-np.pi, np.pi, num=10000)

    # STAR color
    if 3000 < star_teff < 3900:
        star_cl = 'orangered'
    elif 3900 <= star_teff < 5200:
        star_cl = 'goldenrod'
    elif 5200 <= star_teff < 6000:
        star_cl = 'yellow'
    elif 6000 <= star_teff < 7600:
        star_cl = 'palegoldenrod'
    elif star_teff >= 7600:
        star_cl = 'paleturquoise'
    else:
        star_cl = 'goldenrod'

    star_rad = star_rad * 0.00465

    # Set up polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(4, 4))
    ax.axis('off')

    # Orbit
    phases = np.linspace(0, 2 * np.pi, num=10000)
    orbit = orbit_pos(0, a, phases)
    ax.plot(phases, orbit, linestyle='dashdot', color='darkblue', linewidth=1.25, label='Phase')

    # Plot flare phases
    for index, phase in enumerate(flare_phases):
        adj_phase = UT.range_shift(phase, 0, 1, 0, 2 * np.pi)
        lw = 0.75 if len(flare_phases) > 30 else 1.25
        if index == 0:
            ax.vlines(adj_phase, ymin=0, ymax=orbit_pos(0, a, adj_phase),
                      linewidth=lw, color='purple', linestyle='--', label='Flares')
        else:
            ax.vlines(adj_phase, ymin=0, ymax=orbit_pos(0, a, adj_phase),
                      linewidth=lw, color='purple', linestyle='--')

    ax.set_rlim(0, 0.15)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title(title, fontsize=14,y=0.975)

    # VM PDF example
    loc, kappa, TS = MLE.VM_Unbinned_likelihood(flare_phases)
    loc = UT.range_shift(loc, 0, 1, -np.pi, np.pi)
    PDF = vonmises.pdf(thetas, loc=loc + np.pi, kappa=kappa)
    ax.plot(thetas, (PDF / np.max(PDF)) * orbit_pos(0, a, loc + np.pi),
            linestyle='-', color='cyan', label='VM PDF', linewidth=1)
    ax.vlines(np.pi, ymin=0, ymax=orbit_pos(0, a, adj_phase),
              color='black', label='Spot')
    
    # 2. True anomaly at transit
    f_transit = np.pi
    
    # 3. Angle in polar plot:
    theta_center = f_transit
    
    # 4. Edges
    theta1 = theta_center - np.pi/2
    theta2 = theta_center + np.pi/2
    theta_fill = np.linspace(theta1, theta2, 500)
    r_fill = orbit_pos(0, a, theta_fill)
    ax.fill_between(theta_fill, 0, r_fill, color='gray', alpha=0.7, label='Spot LoS Cone')
    ax.vlines(theta_center, ymin=0, ymax=orbit_pos(0,a,theta_center), color='black', linestyle='-', linewidth=1.0)
    plt.legend(loc = 'upper left', handlelength = 2, ncol=2)
    # Fill star
    ax.fill(thetas, star_rad * np.ones(len(thetas)), zorder=4, color=star_cl)
    # ======= DYNAMIC BOTTOM TEXT =======
    if bottom_text_list:
        n = len(bottom_text_list)
        for i, text in enumerate(bottom_text_list):
            x = 0.3 + i * 0.5 / max(n - 1, 1)
            y = 0.15  # adjust this for spacing below the plot
            fig.text(x, y, text, ha='center', va='center', fontsize=11)

    if save_dir:
        plt.savefig(f"{save_dir}/{title}_Polar_Hist.png", dpi=500, bbox_inches='tight')
    plt.show()

def orbit_pos(e,a,theta):
    return a*(1-e**2)/(1+e*np.cos(theta + np.pi))
def a_from_period(R, period):
    return (((period*24*60*60)**2*6.67e-11*(1.989e30)*(R)**(5/4))/(4*np.pi**2))**(1/3)/(1.496e11)
def Polar_Flare_Plot(flare_phases, e = 0.01, a = 0.065, omega_p = 0, title = "Title", star_rad = 1, color_map = 'Reds', peri_err = None,
                     star_teff = 4000, transit = False, save_dir = None, alfven=[0.1, 0.025, 0.025], 
                     TOI = False, emax = 0, roche = [], dynamic_lim = False, legend = False,bottom_text_list=None, alfven_bool = None):
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
    orbit = orbit_pos(e,a,phases)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize = (3.5,3.5))
    ax.axis('off')
    for index, phase in enumerate(flare_phases):
        adj_phase = UT.range_shift(phase, 0, 1, -np.pi, np.pi)
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
    PDF = vonmises.pdf(thetas, loc=loc, kappa=kappa)
    ax.plot(thetas, (PDF/np.max(PDF))*orbit_pos(e,a,loc), linestyle = '-', color='cyan', label = 'VM PDF', linewidth=1)
    plt.fill(thetas, star_rad*np.ones(len(thetas)), zorder=4, c=star_cl)
    omega_p = np.deg2rad(omega_p)
    
    # 2. True anomaly at transit
    f_transit = np.pi/2 - omega_p
    
    # 3. Angle in polar plot:
    theta_center = f_transit
    
    # 4. Edges
    theta1 = theta_center - np.pi/2
    theta2 = theta_center + np.pi/2
    
    # 5. Grid for fill
    theta_fill = np.linspace(theta1, theta2, 500)
    r_fill = orbit_pos(e, a, theta_fill)
    # 6. Fill it
    ax.fill_between(theta_fill, 0, r_fill, color='gray', alpha=0.7, label='Line-of-Sight Cone')
    ax.vlines(theta_center, ymin=0, ymax=orbit_pos(e,a,theta_center), color='black', linestyle='-', linewidth=1.0, label = r"$t_{{mid}}$")
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
    # ======= DYNAMIC BOTTOM TEXT =======
    if bottom_text_list:
        n = len(bottom_text_list)
        for i, text in enumerate(bottom_text_list):
            x = 0.3 + i * 0.5 / max(n - 1, 1)
            y = 0.15  # adjust this for spacing below the plot
            fig.text(x, y, text, ha='center', va='center', fontsize=11)
    if peri_err != None and transit == False:
        if len(peri_err) == 1:
            ax.fill_between(np.linspace(np.pi-peri_err[0], np.pi+peri_err[0]), 0 , orbit_pos(e,a,np.pi), color='black', alpha=0.5,zorder=-1, label='$\delta$ Periastron')
        elif len(peri_err) == 2:
            ax.fill_between(np.linspace(np.pi-np.abs(peri_err[0]), np.pi+np.abs(peri_err[1])), 0 , orbit_pos(e,a,np.pi), color='black', alpha=0.5,zorder=-1, label='$\delta$ Periastron')
    if legend == True:
        ax.legend(handlelength = 1.5, columnspacing=1, labelspacing = .2, prop=dict(family='Serif', size=9), ncol=2, loc = 'upper center')
    if save_dir != None:
        plt.savefig(save_dir + '/' + title + '_Polar_Hist.png', dpi=500, bbox_inches='tight')
        plt.show()