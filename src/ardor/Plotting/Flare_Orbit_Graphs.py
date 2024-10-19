# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:27:30 2024

@author: whitsett.n
"""

# flares_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_MCMC_Flares.csv'
# exo_system_dir = 'C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Flare_Catalog_Files/All_Exoplanets/All_Exoplanet_Parameter_Reference.csv'

# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(subplot_kw={'projection': 'polar'})

# flares = pd.read_csv(flares_dir)
# exo_systems = pd.read_csv(exo_system_dir)
# HD76700 = Star(0.99, 1, .257188, 1.37, age=11e9)
# HD76700b = Planet(0.856, 3.97097, 0.05, 0.090, B=100, stellar_mass=0.99, orbit_resolution=0.005, arg_periastron=180)
# HD76700_flare_epochs = flares.loc[flares['Host_ID'] == 'HD76700']
# epochs = np.array(HD76700_flare_epochs['Periastron_Phase'])
# plot_orbit_flares(HD76700, HD76700b, epochs, fig = fig, ax = ax, alfven_radius=True, overlapping_plots=True)


# HD215497 = Star(0.67, 1, -.41, 1.37, age=7e9)
# HD215497b = Planet(0.214,3.93404, 0.05, 0.16, B=100, stellar_mass=0.67, orbit_resolution=0.005, arg_periastron=180)
# HD215497_flare_epochs = flares.loc[flares['Host_ID'] == 'HD215497']
# epochs2 = np.array(HD215497_flare_epochs['Periastron_Phase'])
# plot_orbit_flares(HD215497, HD215497b, epochs2, fig = fig, ax = ax, alfven_radius=True, overlapping_plots = True)

# AUMic = Star(0.5, 1, -0.99, radius=0.82, age=0.085, alfven = 0.133)
# AUMicb = Planet(1.1297182611,5.6335, 0.06814, 0.186, B=100, stellar_mass=0.5, orbit_resolution=0.005, arg_periastron=180)
# AUMic_flare_epochs = flares.loc[flares['Host_ID'] == 'AUMic']
# epochs4 = np.array(AUMic_flare_epochs['Periastron_Phase'])
# plot_orbit_flares(AUMic, AUMicb, epochs4, fig = fig, ax = ax, alfven_radius=True, overlapping_plots = True)


# HATP2 = Star(1.33, 1, 0.53, radius=1.39, age=1.44e9)
# HATP2b = Planet(1.1297182611,5.6335, 0.06814, 0.5172, B=100, stellar_mass=1.33, orbit_resolution=0.005, arg_periastron=180)
# HATP2_flare_epochs = flares.loc[flares['Host_ID'] == 'HAT-P-2']
# epochs3 = np.array(HATP2_flare_epochs['Periastron_Phase'])
# plot_orbit_flares(HATP2, HATP2b, epochs3, fig = fig, ax = ax, alfven_radius=True, overlapping_plots = False)