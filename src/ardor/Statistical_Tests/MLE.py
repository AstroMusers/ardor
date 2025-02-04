# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:29:39 2024

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import ardor.Utils.Utils as U
import ardor.Statistical_Tests.K_Tests as K
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
    return thetas, vonmises.pdf(thetas, loc = results.x[0], kappa = results.x[1])



def VM_Unbinned_likelihood_List(DataFrame, output_dir):
    hosts = list(set(DataFrame["Host_ID"]))
    host_list = []
    kappa = []
    TS = []
    loc = []
    new_ID = []
    for host in hosts:
        flares = DataFrame.loc[DataFrame["Host_ID"] == host, "Transit_Phase"]
        flares = U.range_shift(flares, 0, 1, -np.pi, np.pi)
        if np.isnan(np.mean(flares)) == False and len(flares) > 3:
            results = minimize(VM_likelihood, ([0, 2]), args = (flares), bounds = [(-np.pi, np.pi), (0, None)])
            null = null_likelihood(flares)
            if null == results.fun:
                print(1)
            if host == 'AUMic':
                print(host, results)
                print(flares)
            host_list.append(host)
            loc.append(round_to_sf(U.range_shift(results.x[0], -np.pi, np.pi, 0, 1), 3))
            kappa.append(round_to_sf(results.x[1], 3))
            TS.append(round_to_sf(2*(null - results.fun), 2))
            new_ID.append(host)
    new_data = pd.DataFrame({"Host_ID": new_ID, "Kappa": kappa, "Center": loc, "TS_{VM}":TS})
    new_data.to_csv(output_dir, index=False)
def Cubic_Unbinned_likelihood_List(DataFrame, output_dir, star, planet):
    hosts = list(set(DataFrame["Host_ID"]))
    host_list = []
    e = []
    ratio = []
    loc = []
    new_ID = []
    TS= []
    for host in hosts:
        flares = DataFrame.loc[DataFrame["Host_ID"] == host, "Periastron_Phase"]
        if np.isnan(np.mean(flares)) == False and len(flares) > 3:
            results = minimize(cubic_likelihood, ([10, 0.5,0.1]), args = (flares, star, planet), bounds = [(1e-10, 1e8), (0, 1), (0,0.999)], tol=1e-8)
            null = null_likelihood(flares)
            host_list.append(host)
            loc.append(results.x[1])
            ratio.append(results.x[0])
            e.append(results.x[2])
            TS.append(round_to_sf(2*(null - results.fun), 2))
            new_ID.append(host[-3:])
    new_data = pd.DataFrame({"Host_ID": new_ID, "ratio": ratio, "Center": loc, "e":e, "TS_{VM}":TS})
    new_data.to_csv(output_dir, index=False)



data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_New.csv")
VM_Unbinned_likelihood_List(data, "C:/Users/Nate Whitsett/Desktop/VM_Likelihood_All_TOI.csv")

