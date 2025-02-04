# -*- coding: utf-8 -*-
"""
Created on Sun May 14 18:28:00 2023

@author: Nate Whitsett
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import simpson

##This defines the stellar parameters relevant in a magnetic star-planet induced flare.
##The parameters are mass (in solar masses), distance (in parsecs), luminosity
##(in log(L_{solar})), spectral class (O,B,A,F,G,K,M), and optional parameters are
##radius (solar radii), age (seconds, 1e9 -> 1 bly), rotational period (days),
##and surface magnetic field strength (Gauss)

##For the most part, all parameters are converted to CGS units if they are not
##already, except for luminosity, which is kept in Watts. This is because the
##interaction energy is quoted in Joules/Watts. It is converted later. 
##Typical terminal stellar windspeed is estimated to be approximately solar-like
##at 400 km/s.

##The mass-loss rate is assumed to be primarily a function of spectral type,
##with more massive, hotter stars undergoing more mass loss, and cooler, 
##less massive and convective stars, losing less. This is based on empirical
##trends from various sources.

##Most stellar magnetic fields are not known; thus, it is approximated by 
##stellar age, following an empirical power law of -1.32. Additionally,
##magnetic fields can be approximated by another empirical power law of
##power -0.655. 

##A crucial assumption of this code is that the stellar magnetic fields are
##approximately dipolar; consequently, the Alfven surface can be estimated
##using only the 'magnetic confinement parameter', eta. This is a function of
##surface stellar magnetic field strength, the star's radius, the terminal
##windspeed, and the mass loss rate. The Alfven surface is then determined by
##a 1/4 power law with respect to eta, which is scaled to the Sun's Alfven
##surface of 18.8 solar radii. The Alfven surface is given in AU.
class Star:
    def __init__(self, mass, dist, lumin, radius=None, age=None, p_rot=None, B=None, err = False, alfven = None):
        """
        The `Star` class initializes with basic physical parameters found on, e.g., the Exoplanet Archive
        and generates the predicted polar magnetic field strength as well as the median Alfven surface
        radius based on various empirical trends found in literature. The results generally carry large 
        uncertainties, and should be taken with caution, particularly the alfven surface estimates.

        Parameters
        ----------
        mass : float
            Stellar mass. In solar masses.
        dist : float
            Stellar distance. In pc.
        lumin : float
            Stellar luminosity. In log(solar).
        radius : float, optional
            Radius of the star. In solar radii. If not provided, approximated
            by mass/radius relation.
        age : float, optional
            Stellar age, in years. 1 Gy = 1e9
        p_rot : float, optional
            Stellar rotation speed, in km/s
        B : float, optional
            Stellar polar field strength, in Gauss
        err : bool, optional
            Randomly generates values from the distributions in the approximations
            used in this code. Can be used to generate, e.g., confidence
            intervals in a 'for' loop.
        alfven : float, optional
            The median Alfven radius, in AU. Can be approximated by passing
            the parameters in the Star() class. NOTE: Results generally have
            high uncertainty.

        Returns
        -------
        Star class object.

        """
        if err == False:
            self.mass = mass
            self.stlumin = 10**(lumin)
            self.lumin = 10**(lumin)*3.8e26
            if radius == None:
                radius = mass**0.8
                self.radius = (mass**0.8)*6.957e10  
            self.radius = radius*6.957e10
            self.age = age
            self.dist = dist*3.086e18
            self.p_rot = p_rot
            self.windspeed = np.sqrt(6.67e-8*self.mass*1.989e33/self.radius)
            self.brightness = self.lumin/(4*np.pi*(self.dist)**2)
            ##Uncertainty of 0.5
            if self.age == None:
                self.massloss = (10**8.33)/(4.5e9**2.33)
            elif self.age != None:
                self.massloss = (10**8.33)/(self.age**2.33)
            if B != None:
                self.B = B
            if age == None and B == None and p_rot != None:
                self.B = 10**1.98/(p_rot**(1.32))
            if p_rot == None and B == None and self.age != None:
                self.B = 10**6.63/(self.age**0.655)
            elif age != None and p_rot != None and B == None:
                self.B = (10**1.98/(0.25*p_rot**1.32)+ 0.75*10**6.63/(age**0.655))
            
            if B != None and self.windspeed != None and self.massloss != None:
                self.eta = ((self.B)**2 * (self.radius)**2)/((self.windspeed)*(self.massloss*6.306*10**25))
            if alfven != None:
                self.Alfven = alfven
            else:
                self.Alfven = self.radius*6.68459e-14*(0.3+(self.eta+0.25)**(1/4))*2.1
        elif err == True:
            self.mass = mass
            self.stlumin = 10**(lumin)
            self.lumin = 10**(lumin)*3.8e26
            if radius == None:
                radius = mass**0.8
                self.radius = (mass**0.8)*6.957e10  
            self.radius = radius*6.957e10
            self.age = age
            self.dist = dist*3.086e18
            self.p_rot = p_rot
            self.windspeed = np.sqrt(6.67e-8*self.mass*1.989e33/self.radius)
            self.brightness = self.lumin/(4*np.pi*(self.dist)**2)
            ##Uncertainty of 0.5
            if self.age == None:
                self.massloss = (10**8.33)/(4.5e9**2.33)
            elif self.age != None:
                self.massloss = (10**8.33)/(self.age**2.33)
            if B != None:
                self.B = B
            if age == None and B == None:
                ##Error of 0.07
                self.B = 10**1.98/(p_rot**(np.random.normal(loc=1.32, scale = 0.14)))
            if p_rot == None and B == None:
                ##Error of 0.0225
                self.B = 10**6.63/(age**(np.random.normal(loc = 0.655, scale = 0.045)))
            if p_rot == None and age == None and B == None:
                return 'Input age, rotational period, or B field'
            elif age != None and p_rot != None and B == None:
                ##Error 0.655
                self.B = (10**1.98/(0.25*p_rot**(np.random.normal(loc=1.32, scale = 0.07))))+ 0.75*10**6.63/(age**(np.random.normal(loc = 0.655, scale = 0.0225)))
            self.eta = ((self.B)**2 * (self.radius)**2)/((self.windspeed)*(self.massloss*6.306*10**25))
            self.Alfven = self.radius*6.68459e-14*(0.29+(self.eta+0.25)**(1/4))*(np.random.normal(loc = 11, scale = 1)/5.3)
        
        
##This defines the planet. Since little is known about the parameters which
##dictate exoplanetary magnetic fields, this is manually inputted. The radius,
##eccentricity and semi-major axis, as well as the occultation phase, are the
##only other planetary parameters. Since this class keeps track of the orbital
##distance to the barycenter (assumed to be the host star), an 'orbital resolution'
##is defined which allows the user to simulate the instrument's phase
##sensitivity. 

##The radius is given in Jupiter radii, the semi-major axis in AU, and the
##magnetic field strength in Gauss. The phase related parameters are in radians.
class Planet:
    def __init__(self, radius, period, a, e, B, arg_periastron=0, orbit_length=100,inclination=90, star = None, compute_period = True):
        '''
        The `Planet` class initializes with basic physical parameters found on, e.g., the Exoplanet Archive
        and generates the orbital geometry of a system, including orbital distances as a function of time 
        and phase which is generated using the Newton-Raphson method.
        Parameters
        ----------
        radius : float
            The radius of the planet, in Jupiter radii.
        period : float
            The period of the planet, in days.
        a : float
            The semi-major axis of the planet, in AU.
        e : float
            The eccentriciy of the planet.
        B : float
            The polar magnetic field strength, in Gauss
        arg_periastron : float, optional
            The argument of periastron, in degrees. The default is 0.
        orbit_resolution : float, optional
            The step size of the planetary orbit. The default is 0.01.
        inclination : float, optional
            The inclination of the orbit, in degrees. The default is 90.
        Star : Star class, optional
            The Star class object of the host of the planet. The default is None.

        Returns
        -------
        Planet class

        '''
        if star == None:
            star = Star(1, 1, 1, B = 1)
        if (star != None and compute_period == True):
            self.period = np.sqrt((4*np.pi**2*(a*1.496e11)**3)/(6.67*10**(-11)*star.mass*1.989e30))/(24*60*60)
        elif period != None:
            self.period = period
        self.radius = radius*7.149e9
        self.e = e
        self.a = a
        self.B = B
        self.orbit_length = orbit_length
        self.arg_periastron = arg_periastron*0.0174533
        if arg_periastron > 90:
            self.true_anomaly = np.pi*2 + np.pi/2 - arg_periastron*0.0174533
        elif arg_periastron <= 90:
            self.true_anomaly = np.pi/2 - arg_periastron*0.0174533
        self.inclination = inclination*0.0174533
        self.periastron = None
        self.periastron_time = None
        self.time, self.position = orbit_pos_v_time(self.period, self.e, self.a, orbit_length = self.orbit_length, phase=False, arg_periastron=arg_periastron)
        self.phase = (self.time/self.period)*np.pi*2
        self.orbit = a*(1-e**2)/(1+e*np.cos(self.phase))
        self.magnetosphere = []
        if Star != None:
            for dist in self.position:
                self.magnetosphere.append(2*(2*1.16)**(1/3)*(self.B/(star.B*(3.4*star.radius/(dist*1.496e+13))**2))*(self.radius/7.149e9))
            self.v = np.sqrt(6.67e-8*star.mass*1.989e+33*(2/(self.position*1.496e+13)-1/(self.a*1.496e+13)))/(1e5)
            self.magnetosphere = np.array(self.magnetosphere)
        

        



##This interaction function takes as input a star and planet class, and will
##generate an array corresponding to a phase-dependent luminosity increase in 
##the host star. The luminosity increase follows the analysis by Antonio
##Lanza of weak-moderate stellar field strength. The energy is proportional
##to (B_{star}^2)(R_{star}^3) * (B_{planet}/B_{star})*(1/(star_planet distance/R_{star})^2)
##That is, the larger the star, the stronger the B field of the star, the
##stronger the B field of the planet to the star, and the closeness of approach
##of the planet all dictate the energy output. 

##The flare itself is given some energy based on the above relation. Then the
##total flare energy is released over some quasi-random timescale
##(corresponding to minutes-hours given a few day orbital period).

##The assumption is that the flare interaction cannot occur if the planet
##is not within its host's Alfven radius. The flare occurence can be modified
##by any scaling factor one wishes. As it stands, the probability is just a 
##linear increse towards some arbitrary probability as the orbit gets closer 
##to the star, though there is no real basis behind it other than intuition.
##The probability is normalized to each step, so increasing the step size
##will not increase the total probability of a flare.
def M_func(e,E,M):
    '''

    Parameters
    ----------
    e : float
        The eccentricity.
    E : float
        The eccentric anomaly.
    M : float
        The mean anomaly.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return E - np.sin(E)*e - M
def M_prime(e,E):
    return 1- np.cos(E)*e

def M_newton(e,M):
    E = M
    for steps in range(100):
        E = E - M_func(e,E,M)/M_prime(e,E)
    return E

def SPI_Metric(distance, B_star, B_planet=10):
    return np.log10(B_star*B_planet/distance**3)

def orbit_pos_v_time(period, e, a, orbit_length = 10, phase=False, arg_periastron = 0):
    time = 0
    time_list = []
    position = []
    while time < period:
        n = np.pi*2/period
        M = n*time
        E = (M_newton(e, M) + 1) + 2.14281
        dist = a*(1-e*np.cos(E))
        time += 1/orbit_length
        time_list.append(time)
        position.append(dist)
    if phase == True:
        return np.array(time_list)/period, np.array(position)
    else:
        return np.array(time_list), np.array(position)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def asymmetric_error(value, upper, lower, lumin_check = False):
    if lumin_check == False:
        while value < 0:
            if np.isnan(value) == True:
                return value
            if np.isnan(lower) == True:
                lower = 0.1*value
            if np.isnan(lower) == True:
                upper = 0.1*value
            check = np.random.random()
            if check > 0.5:
                err = np.abs(value - np.random.normal(loc = value, scale = np.abs(upper)))
                value = value + err
            elif check <= 0.5:
                err = np.abs(value - np.random.normal(loc = value, scale = np.abs(lower)))
                value = value - err
    if lumin_check == True:
        if np.isnan(value) == True:
            return value
        if np.isnan(lower) == True:
            lower = 0.1*value
        if np.isnan(lower) == True:
            upper = 0.1*value
        check = np.random.random()
        if check > 0.5:
            err = np.abs(value - np.random.normal(loc = value, scale = np.abs(upper)))
            value = value + err
        elif check <= 0.5:
            err = np.abs(value - np.random.normal(loc = value, scale = np.abs(lower)))
            value = value - err
    return value
def interaction(star, planet):
    period = planet.period
    flare_time = 24*60*(period)*(2/np.pi)* ((2*1.16)**(1/3) * (((planet.a*1.496e13))/(star.radius))**(-(1/3)) * ((1-planet.e)**(7/6))/((1+planet.e)**(1/2)) * (planet.B/star.B)**(1/3) * (planet.radius)/(star.radius))
    total_energy = (star.B**2*(star.radius)**3)*((((planet.B/star.B)*0.04)/((planet.a*(1-planet.e)*1.496e13)/star.radius)**2))*(1+(1/3)*((planet.a*(1-planet.e))/(2*(2*1.16)**(1/3)*(planet.B/(star.B*(planet.a*(1-planet.e)/star.radius)**(-2)))**(1/3)*planet.radius)))
    return flare_time, total_energy


def phase_curve(star, planet, interaction=None, linear_parameter=0.1):
    base = []
    subtraction = []  
    transit_loc = int(len(planet.orbit)*((planet.true_anomaly)/(2*np.pi)))
    ingress_phase = int(np.arccos(1-(planet.radius*2)**2/(2*(planet.orbit[transit_loc]*1.496e13)**2))/(2*np.pi)*len(planet.orbit))
    transit_phase = int(np.arccos(1-(star.radius*2)**2/(2*(planet.orbit[transit_loc]*1.496e13)**2))/(2*np.pi)*len(planet.orbit))
    transit_displacement = int(len(planet.orbit)*0.5 - (ingress_phase + transit_phase/2))
    limb_darkening = []
    index = 0
    true_anomaly = planet.true_anomaly
    if true_anomaly < np.pi:
        periastron_index = -true_anomaly/np.pi
    if true_anomaly >= np.pi:
        periastron_index = (2*np.pi-true_anomaly)/(np.pi)
    impact_parameter = (planet.a*1.496e13/star.radius)*np.cos(planet.inclination)*((1-planet.e**2)/(1+planet.e*np.sin(planet.arg_periastron)))
    for thetas in np.arange(-np.pi/2, np.pi/2, np.pi/(2*ingress_phase+transit_phase)):
        limb_darkening.append(1-linear_parameter*(1-np.cos(thetas)))
    limb_impact = 1-linear_parameter*(1-np.cos(impact_parameter*(np.pi/2)))
    limb_darkening.append(limb_darkening[0])
    for steps in range(len(planet.phase)):
        base.append(np.random.normal(1, 0.0001))
    for phases in range(len(planet.orbit)):
        if phases < transit_displacement:
            subtraction.append(0)
        elif phases >= transit_displacement and phases < (transit_displacement + ingress_phase):
            subtraction.append((-((phases-transit_displacement)/ingress_phase)*(planet.radius**2/star.radius**2))*limb_darkening[index]*limb_impact)
            index += 1
        elif phases < (ingress_phase + transit_displacement) + transit_phase:
            subtraction.append(-(planet.radius**2/star.radius**2)*limb_darkening[index]*limb_impact)
            index += 1
        elif phases >= (ingress_phase + transit_displacement) + transit_phase and phases <= (2*ingress_phase + transit_displacement + transit_phase):
            subtraction.append(-(((transit_displacement+transit_phase+2*ingress_phase)-phases)/ingress_phase)*(planet.radius**2/star.radius**2)*limb_darkening[index]*limb_impact)
            index += 1
        elif phases > (2*ingress_phase + transit_displacement) + transit_phase:
            subtraction.append(0)
    while len(subtraction) < len(base):
        subtraction.append(0)
    base= np.array(base)
    subtraction = np.array(subtraction)
    curve = np.add(base,subtraction)
    if interaction is not None:
        interaction = np.array(interaction)
        return np.add(interaction,curve), periastron_index+1
    else:
        return curve, periastron_index+1

def probability_density(star, planet, periastron_index, length = 10):
    probability_density = []
    bool_list = []
    for distance in planet.position:
        if distance > star.Alfven:
            bool_list.append(0)
        if distance <= star.Alfven:
            bool_list.append(1)
        probability = 1/(distance)**3
        probability_density.append(probability)
    percent = np.sum(bool_list)/len(bool_list)
    bool_list = np.array(bool_list)
    phase = planet.time/planet.period
    planet.periastron = periastron_index
    index = int(len(planet.orbit)*periastron_index/2)
    rotated_probability = probability_density[-index:] + probability_density[:-index]
    rotated_bool = np.concatenate((bool_list[int(len(bool_list)/2):],  bool_list[:int(len(bool_list)/2)]))
    uniform = np.ones(len(phase))
    rotated_probability = np.array(rotated_probability)
    rotated_probability = (rotated_probability/simpson(rotated_probability))
    probability_dist = rotated_probability*(star.Alfven/min(planet.position))*rotated_bool
    if probability_dist.sum() == 0:
        probability_dist = np.ones(len(probability_dist))
    integral = simpson(probability_dist, x= phase)
    probability_dist = probability_dist/integral
    probability_dist = (probability_dist*(planet.B/star.B)*percent + uniform)/simpson(probability_dist*(planet.B/star.B)*percent+ uniform, x= phase)
    return probability_dist, planet.phase/(2*np.pi)


def interaction_checker(star, planet):
    check = False
    for instances in planet.orbit:
        if instances < star.Alfven:
            check = True
    return check


def orbit_periastron_recenter(time):
    time = np.roll(time, int(len(time)/2))
    return time

def plot_orbit_flares(star, planet, flare_phase_list, fig, ax, alfven_radius = False, overlapping_plots = False):
    phase = planet.phase/(2*np.pi)
    time, position = orbit_pos_v_time(planet.period, planet.e, planet.a, phase=True, time_step=planet.period/len(planet.orbit))
    time = orbit_periastron_recenter(time)
    ax.plot(phase*2*np.pi, planet.orbit2/star.Alfven, linestyle='--', color = 'black')
    if overlapping_plots == False:
        ax.scatter(1000, 1000, marker='x', s=150, color='blue', label = 'Super-Alfvenic')
        ax.scatter(1000, 1000, marker='x', s=150, color='red', label = 'Sub-Alfvenic')
    # if alfven_radius == True:
    #     alfven = plt.Circle((0, 0), star.Alfven, transform=ax.transData._b, fill=True, linewidth=2, zorder=10, color='green', alpha=0.2, label='$\mathrm{R_{A}}$')
    #     plt.gca().add_artist(alfven)
    for phases in flare_phase_list:
        distance = position[find_nearest(time, phases)]
        if phases <= 0.5:
            phase_rad = phase[find_nearest(planet.orbit2[:int(len(planet.orbit2)/2)], distance)]*2*np.pi
        if phases > 0.5:
            phase_rad = phase[find_nearest(planet.orbit2[int(len(planet.orbit2)/2):], distance) ]*2*np.pi
            phase_rad = phase_rad + np.pi
        if alfven_radius == False:
            plt.scatter(phase_rad*2*np.pi, distance, marker='x', s=80)
        else:
            if distance > star.Alfven:
                plt.scatter(phase_rad, distance/star.Alfven, marker='x', s=150, color='blue')
            elif distance <= star.Alfven:
                plt.scatter(phase_rad, distance/star.Alfven, marker='x', s=150, color='red')
    tick_range = np.round(np.linspace(0, np.max(position)*1.1/(star.Alfven), num = 5), 2)
    ax.set_rticks(tick_range)
    ax.set_rlabel_position(-22.5)
    angle = np.deg2rad(67.5)
    L = ax.legend(loc="lower left",
          bbox_to_anchor=(.4 + np.cos(angle)/2, .4 + np.sin(angle)/2), prop={'size': 16})
    plt.setp(L.texts, family='serif')
    ax.set_ylim([0,max(position)*1.1/(star.Alfven)])

