# Star-Planet-Interactions

The scientific goal of this program is to give an expectation of what the phase curve of a transiting exoplanet, which could be inducing stellar flares in its host star, might look like. This is accomplished by consolidating different empirical and theoretical trends which allow common system parameters (e.g. semi-major axis, stellar mass, eccentricity) to be used in generating relevant interaction parameters (e.g. the stellar Alfven Surface, periastron location relative to transit, stellar magnetic field) which are needed in predicting the relative probability and power genarted of an induced flare. This allows for swift identification of salient systems.

The current needs of the science goal only require one file and only basic python packages, and builds from simple Keplerian mechanics and accepted exoplanet parameters. The code makes the star and planet as their own classes with set parameters, with the intention that the input are values from databases such as the Exoplanet Archive. The classes are broken down as such:
******************************************************************************************************************************************************************************************************************************
Star(mass, dist, lumin, radius=None, age=None, p_rot=None, B=None)

Inputs

-mass: Mass of the host star, in solar masses

-dist: Disance of the system from Earth, in parsecs

-lumin: luminosity of the host star, given as log10(Solar Luminosity)

-radius: Optional parameter in solar radii, otherwise estimated by the mass

-age: Optional parameter in seconds, used to estimate stellar magnetic field strength

-p_rot: Optrional parameter in days, used to estimate stellar magnetic field strength

-B: Optional parameter in Gauss, if known for the particular system, otherwise estimated by age and/or the stellar rotation
******************************************************************************************************************************************************************************************************************************
Outputs/Estimaed Parameters

-windspeed: This is the estimated terminal windspeed of the stellar wind, which uses a relation between stellar mass and radius of the star. It is required to estimate the magnetic confinement parameter

-brightness: The observed power of the star in Watts, which uses the distance and luminosity with the inverse square law. This is used to generate the transit phase curve

-massloss: The estimated mass loss of the star, given in solar masses per year. This is again required to estimate the magnetic confinement parameter

-eta: The magnetic confinment parameter. This is a unitless parameter which describes the ratio between the magnetic energy and the bulk kinetic energy of the outflowing stellar wind. The larger the ratio, the more that the star's magnetic field domiantes over the kinetic energy of its stellar wind. The numerator defines the strength of the magnetic confinement and depends on the squre of the radius and surface magnetic field strength, and the denominator defines the kinetic energy of the stellar wind.

-Alfven: This is the estimated Alfven radius in AU, the boundary where perturbations can return back to the stellar surface. This computation assumes the star's magnetic field is dipolar, and that the surface is spherical. It is simply calculated by a power law involving eta and the star's radius.
******************************************************************************************************************************************************************************************************************************
******************************************************************************************************************************************************************************************************************************

Planet(radius, a, e, B, arg_periastron=0, orbit_resolution=0.1,inclination=90)

Inputs:
-radius: The radius of the planet, in Jupiter radii. Used in generating the phase curve

-a: The semi-major axis, in AU. Used in determining the orbit

-e: The eccentricity of the orbit, again used in determining the orbit

-B: Polar magnetic field of the planet. The interaction of the star and planet relies on the presence of a planetary magnetosphere. This value is a free parameter, and can be set to any realistic value.

-arg_periastron: An optional parameter, the argument of periastron, in degrees. Used to determine the orbit orientation relative to the phase curve. While the code can still function without this input, its usefulness is limited in trying to predict where flares should occur in the phase curve of the observed planet.

-orbit_resoltuion: An optional parameter used in determining how many points are used in generating the orbit of the planet, and thus the resolution of the phase curve and any other generated functions

-inclination: An optional parameter, in degrees, used in determining the impact parameter of the planet, which consequently determines the transit depth
******************************************************************************************************************************************************************************************************************************
Outputs/Estimated Parameters:
-phase: an array which goes from 0 to 2pi based on the orbit resolution specified. Used in forming the orbit

-orbit: Takes the Keplerian elliptical orbit which takes semi major axis, eccentricity and angle. This also serves as the planet-star distance

-true_anomaly: This is the angle between when transit occurs and the periastron, in radians. Since the inclination of transiting systems is approximately 90 degrees, this is estimated by pi/2 - argument of periastron. Since the argument of periastron has a range of 0 to 2pi, the code also wraps the angle arround so that the true anomaly is also between 0 and 2pi.

-orbit2: This array rotates the orbit such that the transit occurs at 180 degrees. This allows a polar plot of the system to be generated with 0 degrees corresponding to -1 in the phase curve, and 180 degrees corresponds to transit at 0.

******************************************************************************************************************************************************************************************************************************
******************************************************************************************************************************************************************************************************************************

These two classes must now interact with each other to generate the desired information. This is done through 3 functions: interaction(star, planet), phase_curve(star, planet, interaction, linear_parameter=0.1), and probability(star, planet, periastron_index).These are explained below:

******************************************************************************************************************************************************************************************************************************

interaction(star, planet)

The goal of this function is to generate the characteristics of a theoretical flare of a given planetary system. This includes the total energy of the flare, the timescale of the flare event, the power generated over the timescale, and the brightness increase as seen from the phase curve. 

-Flare Time: The flare time is calculated as a fraction of the planet's orbit, and involves the planet's semi-major axis, the star's radius, the planet's magnetic field,  the star's magnetic field, the planet's radius, and the orbital eccentricity. It is given in seconds.

-Flare energy: The flare energy is a function of the stellar magnetic field, the stellar radius, the ratio of the planet's magnetic field and star's magnetic field, and the star-planet distance. It is given in Joules.

-Flare Power: This is just the flare energy over the flare time. It is given in Watts.

-Brightness Increase: This is an array which acts as an overlay to the transit curve. Most of the time, it is just an array of 0's, but when a flare occurs, a number of arrays are filled in by the expected flare power. The number of blocks the flare occurs is the same ratio as the flare time to the period. 

-Probabalistic model: the rate at which flares occur is quasi-random, since it is unlikely that in a given orbit a close-by planet would induce an event every time. The code itself has a somewhat arbitrary probability at each step of the orbit on whether or not a flare occurs, but it does increase as the planet approaches the star below the Alfven radius. The most important feature is that the probability of a flare occuring is 0 when the planet exists outside of the Alfven radius. This is a realistic assumption, since this surface represents the boundary at which a perturbation can reach the stellar surface. Ultimately, this part of the code is relatively unimportant, since the probability density is addressed more thoroughly in the probability function.

******************************************************************************************************************************************************************************************************************************
******************************************************************************************************************************************************************************************************************************

phase_curve

The goal of this function is to generate a semi realistic phase curve of the transiting planet. It incorporates the geometry of the planet over the stellar disk, taking into account distance at time of transit, ingress and egress, the impact parameter, as well as a linear limb darkening term.
