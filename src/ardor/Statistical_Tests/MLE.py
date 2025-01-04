# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:29:39 2024

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import ardor.Utils.Utils as U
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.stats import vonmises
from scipy.stats import uniform
from matplotlib import pyplot as plt
import os
import numpy as np
import pandas as pd
import math
def round_to_sf(x, sf):
    """
    Rounds a number to a specified number of significant figures.

    Args:
        x: The number to round.
        sf: The number of significant figures.

    Returns:
        The rounded number.
    """
    if x == 0:
        return 0.0
    
    sign = 1 if x > 0 else -1
    x = abs(x)
    
    rounding_position = sf - int(math.floor(math.log10(x))) - 1
    
    rounded_x = round(x, rounding_position)
    
    return rounded_x * sign
def Inverse_Cubic_Function(ratio, loc, e, star, planet):
    model, phase = SPI.SPI_Cubic(ratio, loc,e, star, planet, length=10)
    f = interp1d(np.linspace(0, 1, num=len(model)) ,model, kind='linear')
    return f

def VM_pdf(x, loc, kappa):
    return vonmises.pdf(x, kappa= kappa, loc=loc)
def cubic_pdf(x, ratio, loc, e, star, planet):
    return Inverse_Cubic_Function(ratio, loc,e, star, planet)(x)

def VM_likelihood(params, data):
    loc, kappa = params

    return -np.sum(np.log(VM_pdf(data, loc, kappa)))

def cubic_likelihood(params, data, star, planet):
    ratio, loc,e = params
    return -np.sum(np.log(cubic_pdf(data, ratio, loc, e, star, planet)))
def null_likelihood(data):
    return -np.sum(np.log(uniform.pdf(data, loc = -np.pi, scale=2*np.pi)))

def VM_Unbinned_likelihood(flares):
    results = minimize(VM_likelihood, ([0, 1]), args = (flares), bounds = [(-np.pi, np.pi), (1e-10, None)])
    thetas = np.linspace(0, np.pi*2, num=100)
    print(results)
    return thetas, vonmises.pdf(thetas, loc = results.x[0], kappa = results.x[1])

# data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares_New.csv")
# flares = np.array(data.loc[data['Host_ID'] == 'AUMic', ['Transit_Phase']])
# flares = (((flares - 0) * (2*np.pi)) / 1) + -np.pi

# star = OML.Star(1, 1, 1, radius=1, age=1e9, alfven=0.1)
# planet = OML.Planet(1, 6.7, a=0.0649, e=0.01, B=10)
# hosts = set(np.array(data["Host_ID"]))
# kappa = []
# loc = []
# TS = []
# host_list = []
# SNR_list = []
# e = []
# ratio = []
# loc = []
# fit = "VM"
# runs = os.listdir("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/G_Type_Inverse_Cubic")
# for index, sims in enumerate(runs):
#     if sims[0] == "G":
#         data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/G_Type_Inverse_Cubic/" + sims)
#         print(data)
#         hosts = set(np.array(data['Host_ID']))
    
#         if fit == "VM":
#             for host in hosts:
#                 flares = np.array(data.loc[data['Host_ID'] == host, ['Periastron_Phase']])
#                 flares = U.range_shift(flares, 0, 1, -np.pi, np.pi)
#                 if np.isnan(np.mean(flares)) == False and len(flares) > 3:
#                     results = minimize(VM_likelihood, ([0, 1]), args = (flares), bounds = [(-np.pi, np.pi), (1e-10, None)])
#                     null = null_likelihood(flares)
#                     host_list.append(host)
#                     loc.append(round_to_sf(U.range_shift(results.x[0], -np.pi, np.pi, 0, 1), 3))
#                     kappa.append(round_to_sf(results.x[1], 3))
#                     TS.append(round_to_sf(2*(null - results.fun), 2))
#             new_data = pd.DataFrame({"Host_ID": host_list, "Kappa": kappa, "Center": loc, "TS_{VM}":TS})
#             new_data.to_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/src/ardor/Statistical_Tests/VM_Unbinned_LK_Sims" + str(index) + ".csv", index=False)
    
# if fit == "Cubic":
#     for host in hosts:
#         flares = np.array(data.loc[data['Host_ID'] == host, ['Periastron_Phase']])
#         if np.isnan(np.mean(flares)) == False and len(flares) > 3:
#             star = OML.Star(1, 1, 1, radius=1, age=1e9, alfven=0.1)
#             planet = OML.Planet(1, 6.7, a=0.0649, e=0.01, B=10)
#             results = minimize(cubic_likelihood, ([1, 0.5, 0.01]), args = (flares, star,planet), bounds = [(1e-8,None), (1e-8, 0.999), (0,0.999)])
#             null = null_likelihood(flares)
#             host_list.append(host)
#             e.append(round_to_sf(results.x[2],3))
#             ratio.append(round_to_sf(results.x[0],3))
#             loc.append(round_to_sf(results.x[1],3))
#             TS.append(round_to_sf(2*(null - results.fun), 2))
#     new_data = pd.DataFrame({"Host_ID": host_list, "B_{P}/B_{S}": ratio, "e": e, "Center": loc, "TS_{Cubic}":TS})
#     new_data.to_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Python Scripts/ardor/src/ardor/Statistical_Tests/Cubic_Unbinned_LK.csv", index=False)