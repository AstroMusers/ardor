# -*- coding: utf-8 -*-
"""
Created on Sun May 14 18:28:00 2023

@author: Nate Whitsett
"""

import numpy as np
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
            if mass == None:
                self.mass = radius**(1/0.8)*1.989e33
            elif mass != None:
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
            self.density = (self.mass*1.989e33)/((4/3)*np.pi*(6.957e10*self.radius)**3)
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
        if radius > 0.142744 and radius < 1.11:
            self.mass = radius**(3/2)
        if radius < 0.142744:
            self.mass = radius**(100/27)
        if radius > 1.142744:
            self.mass = radius**(0.06)
        self.density = (self.mass*1.89e30)/((4/3)*np.pi*(self.radius)**3)
        self.inclination = inclination*0.0174533
        self.periastron = None
        self.periastron_time = None
        self.time, self.position, self.rot_time = orbit_pos_v_time(self.period, self.e, self.a, orbit_length = self.orbit_length, phase=False, arg_periastron=arg_periastron)
        self.phase = (self.time/self.period)*np.pi*2
        self.orbit = a*(1-e**2)/(1+e*np.cos(self.phase))
        self.magnetosphere = []
        if star != None:
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
    Function relating mean anomaly and eccentric anomaly.
    Parameters
    ----------
    e : float
        Orbital eccentricity.
    E : float
        The eccentric anomaly, in degrees.
    M : float
        The mean anomaly, in degrees.

    Returns
    -------
    function
        Mean anomaly/eccentriciy equation for NR method.

    '''
    return E - np.sin(E)*e - M
def M_prime(e,E):
    '''
    Derivative of mean anomaly equation to be used in NR method.

    Parameters
    ----------
    e : float
        Orbital eccentricity.
    E : float
        Eccentric anomaly, in degrees.

    Returns
    -------
    function
        Derivative for mean anomaly equation.

    '''
    return 1- np.cos(E)*e

def M_newton(e,M):
    '''
    Newton-Raphson method to solve for eccentric anomaly given mean anomaly.

    Parameters
    ----------
    e : float
        Orbital eccentricity.
    M : float
        Mean anomaly, in degrees.

    Returns
    -------
    E : float
        Eccentric anomaly, in degrees.

    '''
    if e < 0.8:
        E = M
    else:
        E = np.pi
    for steps in range(100):
        E = E - M_func(e,E,M)/M_prime(e,E)
    return E
def true_anomaly(E, e):
    '''
    Computes the true anomaly from the eccentric anomaly.

    Parameters
    ----------
    E : float
        Eccentric anomaly, in degrees.
    e : float
        Orbital eccentricity.

    Returns
    -------
    true anomaly : float
        The true anomaly, in degrees.

    '''
    beta = e/(1+np.sqrt(1-e**2))
    return E + 2*np.arctan((beta*np.sin(E))/(1-beta*np.cos(E)))
def elliptical_dist(a, e, theta):
    '''
    Generates instantenous separation of a planetary orbit as function of 
    true anomaly.

    Parameters
    ----------
    a : float
        Semi-major axis, in AU.
    e : float
        Orbital eccentricity.
    theta : float
        True anomaly, in degrees

    Returns
    -------
    instantanous separation : float
        The instantanous separation at the specified true anomaly, in AU.

    '''
    return a*(1-e**2)/(1+np.cos(theta)*e)
def sol(A, K):
    return A[K % len(A):] + A[:K % len(A)]
def orbit_pos_v_time(period, e, a, orbit_length , phase=False, arg_periastron = 0):
    '''
    Uses Newton-Raphson method to generate a time-dependent instantanous separation
    equation.

    Parameters
    ----------
    period : float
        Orbital period, in days.
    e : float
        Orbital eccentricity.
    a : float
        Semi-major axis, in AU.
    orbit_length : float
        How fine of a grid is returned.
    phase : bool, optional
        Whether or not to return as a function of normalized
        orbital phase (0 to 1), with 0.5 marking periastron. The default is False.
    arg_periastron : float, optional
        The argument of periastron, omega, in degrees. The default is 0.

    Returns
    -------
    time : numpy array
        The time, in either phase or days, assuming periastron
        is mid-phase.
    separation : numpy array
        The instantanous separation, in AU.
    rotated time : numpy array
        The time array rotated to the appropriate argument of periastron.

    '''
    time = 0
    time_list = []
    position = []
    n = np.pi*2/period
    while time < (period):
        M = n*(time)
        E = (M_newton(e, M))
        nu = true_anomaly(E,e)
        time += 1/orbit_length
        time_list.append(time)
        position.append(elliptical_dist(a,e,nu))
    rot = int((arg_periastron/(2*np.pi))*len(time_list))
    rot_time = sol(time_list, rot)
    if phase == True:
        return np.array(time_list)/period, np.array(position), np.array(rot_time)/period
    else:
        return np.array(time_list), np.array(position),np.array(rot_time)


def find_nearest(array, value):
    '''
    Takes in an array and value and returns the index of the closest value
    in the array.

    Parameters
    ----------
    array : np array
        Array to search.
    value : float
        Value to find the closest in the array.

    Returns
    -------
    idx : int
        Index in the array that most closely matches the value passed.

    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def interaction(star, planet):
    '''
    Determines the flare energy predicted by an induced flare following the prescription
    from Lanza 2018 for main sequence host stars.

    Parameters
    ----------
    star : Star class object.
        The host star of the interaction.
    planet : Planet class objected
        The planet of the interaction.

    Returns
    -------
    flare_time : float
        The time scale of SPI induced flare, in seconds.
    total_energy : float
        Total energy of the SPI induced flare, in ergs.

    '''
    period = planet.period
    flare_time = 24*60*(period)*(2/np.pi)* ((2*1.16)**(1/3) * (((planet.a*1.496e13))/(star.radius))**(-(1/3)) * ((1-planet.e)**(7/6))/((1+planet.e)**(1/2)) * (planet.B/star.B)**(1/3) * (planet.radius)/(star.radius))
    total_energy = (star.B**2*(star.radius)**3)*((((planet.B/star.B)*0.04)/((planet.a*(1-planet.e)*1.496e13)/star.radius)**2))*(1+(1/3)*((planet.a*(1-planet.e))/(2*(2*1.16)**(1/3)*(planet.B/(star.B*(planet.a*(1-planet.e)/star.radius)**(-2)))**(1/3)*planet.radius)))
    return flare_time, total_energy


def arg_peri_to_epoch(arg, e, a, period):
    '''
    Computes the delay in time between the transit epoch and the epoch of peri-
    astron given the argument of periastron, eccentricity, semi-major axis,
    and period.

    Parameters
    ----------
    arg : float
        Argument of periastron, in degrees. The arguement of ascending node is 
        when the planet is at quadrature, with the planet approaching the observer.
    transit_epoch : float
        Transit epoch, in BJD.
    e : float
        Orbital eccentricity.
    a : float
        Semi-major axis, in AU.
    period : float
        Orbital period, in days.

    Returns
    -------
    delta_epoch : float
        The time difference between the transit epoch and the periastron epoch.
        Returns in days. 

    '''
    if e == 0 or np.isnan(e) == True:
        print('No argument of periastron!')
    else:
        arg_theta = (arg/360)*2*np.pi
        time, position = orbit_pos_v_time(period, e, a, orbit_length=200, phase=True)
        omegas = np.linspace(0, 2*np.pi, num=len(time))
        orbit = a*(1-e**2)/(1+e*np.cos(omegas-arg_theta))
        transit_theta = np.pi/2
        transit_distance = orbit[find_nearest(omegas, transit_theta)]
        if arg >= 90 and arg < 270:
            transit_time = time[find_nearest(position[:int(len(position)/2)], transit_distance)]
        else:
            transit_time = time[find_nearest(position[int(len(position)/2):], transit_distance) + int(len(position)/2)]
        delta_epoch = transit_time
        return delta_epoch
    
    
def transit_phase_to_peri_shift(arg, e, a, period,num=20, larg = 0, uarg = 0):
    '''
    Takes in the argument of periastron, eccentricity, semi-major axis, and
    period, and returns the shift in normalized phase (0,1) from the transit 
    phase. The transit phase is assumed to be at 0.5. It can also take in
    lower and upper uncertainties of the argument of periastron and return
    the lower and upper uncertainty in the phase shift.
    
    This is computed by comparing two orbital equations:
        - The first computes the instantaneous separation of the planet
        with respect to true anomaly. 
        - The second computes the instantanous separation of the plnaet
        with respect to time.
    To find the phase shift given the argument of periastron (true anomaly),
    we compute the two equations both temporaly and angularly. Based on the
    solution to the angular equation given the true anomaly, we match the 
    temporal solution at the same instantanous separation and find the phase
    shift in time that leads to the same separation solution.
    Parameters
    ----------
    arg : float
        Argument of periastron, in degrees. The arguement of ascending node is 
        when the planet is at quadrature, with the planet approaching the observer.
    transit_epoch : float
        Transit epoch, in BJD.
    e : float
        Orbital eccentricity.
    a : float
        Semi-major axis, in AU.
    period : float
        Orbital period, in days.
        DESCRIPTION.
    num : int, optional
        The number of true anomalies to use to compute the temporal
        grid. The default is 20, evenly spaced between 0 and 360.
    larg : float, optional
        Lower uncertainty of omega, in degrees. The default is 0.
    uarg : float, optional
        Upper uncertainty of omega, in degrees. The default is 0.

    Returns
    -------
    cent_shift : float
        Shift in phase from transit phase to periastron phase.
    lshift : float
        Lower error shift in phase.
    ushift : float
        Upper error shift in phase.

    '''
    lower_arg = arg - larg
    upper_arg = arg+ uarg
    if lower_arg < 0:
        lower_arg += 1
    if lower_arg > 1:
        lower_arg -= 1
    if upper_arg < 0:
        upper_arg += 1
    if upper_arg > 1:
        upper_arg -= 1
    if lower_arg > upper_arg:
        lower_arg, upper_arg = upper_arg, lower_arg
    omegas = np.linspace(lower_arg, upper_arg, num =num)
    shift = []
    for omega in omegas:
        a = arg_peri_to_epoch(omega, e, a, period, larg = larg, uarg = uarg)
        shift.append(a)
    cent_shift = shift[find_nearest(omegas, arg)]
    lower_shift = shift[find_nearest(omegas, lower_arg)]
    upper_shift = shift[find_nearest(omegas, upper_arg)]
    if cent_shift > 1:
        cent_shift -= 1
    lshift = np.abs(lower_shift-cent_shift)
    ushift = np.abs(cent_shift - upper_shift)
    return cent_shift, lshift, ushift
