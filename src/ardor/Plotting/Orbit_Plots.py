# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:56:10 2024

@author: natha
"""
import Orbit_Model_Library as OML
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import patches as pt
import matplotlib.animation as animation


###############################################################################
############ PLOTS FOR INVERSE CUBIC PROBABILITY DENSITY ######################
###############################################################################
font = {'family': 'serif',
        'weight': 'normal',
        'size': 11,
        }
fig = plt.figure(figsize = (7,7), dpi=100)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222, projection = 'polar')
ax3 = plt.subplot(223)
ax4 = plt.subplot(224, projection = 'polar')
plt.subplots_adjust(wspace=0, hspace=0)
P_B = [0, 2.5, 5, 10, 50]
P_e = [0, 0.1, 0.25, 0.5]
colors = ['blue', 'orange', 'green', 'red']
# P_a = 
# P_T = 
mpl.rc('font',family='serif', size=9)
phase_lists = []
density_lists = []
orbit_lists = []
theta_lists = []
for ratios in P_B:
    star = OML.Star(1, 1, 1, radius=1, age=1e9, B = 5)
    planet = OML.Planet(1, 7, 0.03, 0.3, B=ratios, Star=star)
    phase, density = OML.probability_density(star, planet, 1)
    if ratios < 5 and ratios > 0:
        ax1.plot(phase, density, label = round(ratios/5, 2), linestyle = '--')
    elif ratios == 0:
        ax1.plot(phase, density, label = round(ratios/5, 2))
    else:
        ax1.plot(phase, density, label = int(ratios/5), linestyle = '--')
    ax2.plot(planet.phase, planet.orbit, linestyle = '--', color = 'black')
    ax2.add_patch(pt.Circle((0,0), 0.004, transform=ax2.transData._b, alpha = 0.6, color = 'orange'))
    # ax2.add_patch(pt.Circle((min(planet.position),0), 0.002, transform=ax2.transData._b, alpha = 0.6, color='purple'))
    phase_list = []
    density_list = []
    for i in range(50):
        phase_list.append(phase[int(len(planet.phase)/50*i)])
        density_list.append(density[int(len(planet.phase)/50*i)])
    phase_lists.append(phase_list)
    density_lists.append(density_list)
    
for index, ecc in enumerate(P_e):
    star = OML.Star(1, 1, 1, radius=1, age=1e9, B = 5)
    planet = OML.Planet(1, 7, 0.03, ecc, B=10, Star=star)
    phase, density = OML.probability_density(star, planet, 1)
    if ecc == 0:
        ax3.plot(phase, density, label = ecc)
    else:
        ax3.plot(phase, density, label = ecc, linestyle = '--')
    phase_list = []
    density_list = []
    orbit_list = []
    theta_list= []
    for i in range(50):
        if i < 10:            
            phase_list.append(phase[int(len(planet.phase)/70*i)])
            density_list.append(density[int(len(planet.phase)/70*i)])
            orbit_list.append(planet.orbit[int(len(planet.phase)/70*(i-25))])
            theta_list.append(planet.phase[int(len(planet.phase)/70*(i-25))])
        if i > 10 and i <= 20:            
            phase_list.append(phase[int(len(planet.phase)/50*i)])
            density_list.append(density[int(len(planet.phase)/50*i)])
            orbit_list.append(planet.orbit[int(len(planet.phase)/50*(i-25))])
            theta_list.append(planet.phase[int(len(planet.phase)/50*(i-25))])
        if i > 20 and i <= 30:            
            phase_list.append(phase[int(len(planet.phase)/50*i)])
            density_list.append(density[int(len(planet.phase)/50*i)])
            orbit_list.append(planet.orbit[int(len(planet.phase)/50*(i-25))])
            theta_list.append(planet.phase[int(len(planet.phase)/50*(i-25))])
        if i > 40:            
            phase_list.append(phase[int(len(planet.phase)/50*i)])
            density_list.append(density[int(len(planet.phase)/50*i)])
            orbit_list.append(planet.orbit[int(len(planet.phase)/50*(i-25))])
            theta_list.append(planet.phase[int(len(planet.phase)/50*(i-25))])
    phase_lists.append(phase_list)
    density_lists.append(density_list)
    orbit_lists.append(orbit_list)
    theta_lists.append(theta_list)
    ax4.plot(planet.phase, planet.orbit, linestyle = '--')
    ax4.add_patch(pt.Circle((0,0), 0.004, transform=ax4.transData._b, alpha = 0.6, color = 'orange'))
    # ax4.add_patch(pt.Circle((min(planet.position),0), 0.002, transform=ax4.transData._b, alpha = 0.6, color=colors[index]))



### ANIMATE THE PLOTS
artists = []
star = OML.Star(1, 1, 1, radius=1, age=1e9, B = 5)
planet = OML.Planet(1, 7, 0.03, 0.3, B=ratios, Star=star)
for i in range(50):
    
    ## TOP PLOT
    theta = planet.phase[int(len(planet.phase)/50*(i-25))]
    radius = planet.orbit[int(len(planet.phase)/50*(i-25))]
    data, = ax2.plot(theta, radius, 'o', color='black')
    data2, = ax1.plot(phase_lists[0][i], density_lists[0][i], 'o', color='blue')
    data3, = ax1.plot(phase_lists[1][i], density_lists[1][i], 'o', color='orange')
    data4, = ax1.plot(phase_lists[2][i], density_lists[2][i], 'o', color='green')
    data5, = ax1.plot(phase_lists[3][i], density_lists[3][i], 'o', color='red')
    data6, = ax1.plot(phase_lists[4][i], density_lists[4][i], 'o', color='purple')
    
    ## BOTTOM PLOT
    data7, = ax3.plot(phase_lists[5][i], density_lists[5][i], 'o', color='blue')
    data8, = ax3.plot(phase_lists[6][i], density_lists[6][i], 'o', color='orange')
    data9, = ax3.plot(phase_lists[7][i], density_lists[7][i], 'o', color='green')
    data10, = ax3.plot(phase_lists[8][i], density_lists[8][i], 'o', color='red')
    data11, = ax4.plot(theta_lists[0][i], orbit_lists[0][i], 'o', color='blue')
    data12, = ax4.plot(theta_lists[1][i], orbit_lists[1][i], 'o', color='orange')
    data13, = ax4.plot(theta_lists[2][i], orbit_lists[2][i], 'o', color='green')
    data14, = ax4.plot(theta_lists[3][i], orbit_lists[3][i], 'o', color='red')
    artists.append([data,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14])
    
    
### Labels and legends for the regular plots
ax3.set_xlabel('Phase')
fig.text(0.05, 0.5, 'Probability Density', va='center', rotation='vertical', 
          fontdict=font)
ax1.legend(borderpad = 0.5, labelspacing=0.25, title = '$\mathrm{B_{pl}/B_{star}}$',
            handlelength = 0.75)
ax3.legend(borderpad = 0.5, labelspacing=0.25, title = 'e',
            handlelength = 0.75)
ax1.text(.02, .97, '(a)', ha='left', va='top', transform=ax1.transAxes,
            fontdict = font)
ax2.text(.02, .97, '(b)', ha='left', va='top', transform=ax2.transAxes,
            fontdict = font)
ax3.text(.02, .97, '(c)', ha='left', va='top', transform=ax3.transAxes,
            fontdict = font)
ax4.text(.02, .97, '(d)', ha='left', va='top', transform=ax4.transAxes,
            fontdict = font)


### Set r labels and maximums for polar plots ###
ax2.set_rmax(0.05)
ax2.set_rticks([0.01, 0.02, 0.03, 0.04, 0.05])
ax4.set_rmax(0.05)
ax4.set_rticks([0.01, 0.02, 0.03, 0.04, 0.05])
ax2.tick_params(axis='y', labelsize=9)
ax4.tick_params(axis='y', labelsize=9)
ax2.set_xticks([])
ax4.set_xticks([])

###Plot Lines denoting Alfven surface
ax2.add_patch(pt.Circle((0,0), star.Alfven, transform=ax2.transData._b, color='firebrick', linestyle = ':', fill=False, label = 'Alfven Surface', linewidth=2))
ax4.add_patch(pt.Circle((0,0), star.Alfven, transform=ax4.transData._b, color='firebrick', linestyle = ':', fill=False, label = 'Alfven Surface', linewidth=2))
ax2.text(0.6, 0.051, 'AU')
ax4.text(0.6, 0.051, 'AU')
ax2.text(1.95, 0.045, 'Alfvén Surface', color = 'firebrick', fontdict={'family': 'serif'})
ax4.text(1.95, 0.045, 'Alfvén Surface', color = 'firebrick', fontdict={'family': 'serif'})
# plt.savefig('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Publication_2024/Inverse_Cubic.png',
#             dpi=800, bbox_inches='tight')
ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=100)
ani.save(filename="C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/Repos/ardor/test.gif", writer="pillow")  
