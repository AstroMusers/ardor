# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:29:39 2024

@author: whitsett.n
"""

import ardor.SPI_Forward_Models.SPI_Simulation as SPI
import ardor.SPI_Forward_Models.Orbit_Model_Library as OML
import ardor.Flares.Flare as Flare
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
def Inverse_Cubic_Function(ratio, loc, e, a, period):
    model, phase = SPI.SPI_Cubic(ratio, loc, e, a, period, length=10)
    f = interp1d(np.linspace(0, 2*np.pi, num=len(model)) ,model, kind='linear')
    return f
def VM_pdf(x, loc, kappa):
    return vonmises.pdf(x, kappa= kappa, loc=loc)
def cubic_pdf(x, ratio, loc, e, a, period):
    return Inverse_Cubic_Function(ratio, loc, e, a, period)(x)
def VM_likelihood(params, data):
    loc, kappa = params

    return -np.sum(np.log(VM_pdf(data, loc, kappa)))
def cubic_likelihood(params,data, e, a,period):
    loc, ratio = params
    return -np.sum(np.log(cubic_pdf(data, ratio, loc, e, a, period)))
def null_likelihood(data):
    return -np.sum(np.log(uniform.pdf(data, loc = -np.pi, scale=2*np.pi)))
def null_likelihood_cubic(data):
    return -np.sum(np.log(uniform.pdf(data, loc = 0, scale=2*np.pi)))
def VM_Unbinned_likelihood(flares):
    flares = U.range_shift(flares, 0, 1, -np.pi, np.pi)
    if np.isnan(np.mean(flares)) == False and len(flares) > 3:
        results = minimize(VM_likelihood, ([0, 2]), args = (flares), bounds = [(-np.pi, np.pi), (0, None)])
        null = null_likelihood(flares)
        if null == results.fun:
            print(1)
        loc = round_to_sf(U.range_shift(results.x[0], -np.pi, np.pi, 0, 1), 3)
        kappa = round_to_sf(results.x[1], 3)
        TS = 2*(null - results.fun)
    return loc, kappa, np.sqrt(TS)
def Cubic_Unbinned_likelihood(flares,e,a,period):
    flares = U.range_shift(flares, 0, 1, 0, 2*np.pi)
    results = minimize(cubic_likelihood, ([0, 0]), args = (flares,e, a, period), bounds = [(0, 2*np.pi), (0,1e4)], tol=1e-8)
    null = null_likelihood_cubic(flares)
    loc = results.x[0]
    ratio = results.x[1]
    TS =2*(null - results.fun)
    return loc, ratio, TS
def VM_Unbinned_likelihood_List(DataFrame, output_dir, phase_header="Transit_Phase"):
    hosts = list(set(DataFrame["Host_ID"]))
    host_list = []
    kappa = []
    TS = []
    loc = []
    new_ID = []
    for host in hosts:
        flares = DataFrame.loc[DataFrame["Host_ID"] == host, phase_header]
        flares = U.range_shift(flares, 0, 1, -np.pi, np.pi)
        if np.isnan(np.mean(flares)) == False and len(flares) > 3:
            results = minimize(VM_likelihood, ([0, 2]), args = (flares), bounds = [(-np.pi, np.pi), (0, None)])
            null = null_likelihood(flares)
            if null == results.fun:
                print(1)
            host_list.append(host)
            loc.append(round_to_sf(U.range_shift(results.x[0], -np.pi, np.pi, 0, 1), 3))
            kappa.append(round_to_sf(results.x[1], 3))
            TS.append(round_to_sf(2*(null - results.fun), 2))
            new_ID.append(host)
    new_data = pd.DataFrame({"Host_ID": new_ID, "Kappa": kappa, "Center": loc, "TS_{VM}":TS})
    new_data.to_csv(output_dir, index=False)
def Cubic_Unbinned_likelihood_List(DataFrame, output_dir, e, a,period, TOI = False):
    hosts = list(set(DataFrame["Host_ID"]))
    host_list = []
    ratio = []
    loc = []
    new_ID = []
    TS= []
    for host in hosts:
        flares = DataFrame.loc[DataFrame["Host_ID"] == host, "Periastron_Phase"]
        results = minimize(cubic_likelihood, ([10, 5]), args = (flares,e, a, period), bounds = [(1e-10, 1e8), (0,2*np.pi)], tol=1e-8)
        null = null_likelihood_cubic(flares)
        host_list.append(host)
        loc.append(results.x[1])
        ratio.append(results.x[0])
        TS.append(round_to_sf(2*(null - results.fun), 2))
        new_ID.append(host[-3:])
    new_data = pd.DataFrame({"Host_ID": new_ID, "ratio": ratio, "Center": loc, "TS_{VM}":TS})
    new_data.to_csv(output_dir, index=False)
# print(Cubic_Unbinned_likelihood(np.array([0.14396843 ,0.1899146  ,0.34667978 ,0.35896735 ,0.3778688 , 0.43336178,
#  0.92344073 ,0.96847297]), 0.5, 0.3, 2))
def UBL_Sim_Kappa(DataFrame, output_dir, star=''):
    # interest_hosts = [['p_KU', 'p_KS', 'p_AD', 'p_KS_Samp', 'p_AD_Samp', 'N']]
    data_frame = pd.DataFrame()
    set_hosts = set(DataFrame['Host_ID'])
    kappa0_kappa = []
    kappa25_kappa = []
    kappa5_kappa = []
    kappa1_kappa = []
    kappa2_kappa = []
    kappa4_kappa = []
    kappa8_kappa = []
    loc0_kappa = []
    loc25_kappa = []
    loc5_kappa = []
    loc1_kappa = []
    loc2_kappa = []
    loc4_kappa = []
    loc8_kappa = []
    TS0_kappa = []
    TS25_kappa = []
    TS5_kappa = []
    TS1_kappa = []
    TS2_kappa = []
    TS4_kappa = []
    TS8_kappa = []
    kappa_list = [kappa0_kappa,
        kappa25_kappa,
        kappa5_kappa,
        kappa1_kappa,
        kappa2_kappa,
        kappa4_kappa,
        kappa8_kappa]
    loc_list = [loc0_kappa,
        loc25_kappa,
        loc5_kappa,
        loc1_kappa,
        loc2_kappa,
        loc4_kappa,
        loc8_kappa]
    TS_list = [TS0_kappa,
        TS25_kappa,
        TS5_kappa,
        TS1_kappa,
        TS2_kappa,
        TS4_kappa,
        TS8_kappa]
    for index, hosts in enumerate(set_hosts):

        if float(hosts[-1]) == 8 or float(hosts[-1]) == 4 or float(hosts[-1]) == 0 or float(hosts[-1]) == 2 or float(hosts[-1]) == 1:
            kappa = float(hosts[-1])
        elif float(hosts[-3:]) == 0.5:
            kappa = 0.5
        elif float(hosts[-4:]) == 0.25:
            kappa = 0.25
        peri_phases = np.array(DataFrame.loc[DataFrame['Host_ID'] == hosts, 'Phase'])
        peri_phases = Utils.range_shift(peri_phases, 0, 1, -np.pi, np.pi)
        results = MLE.VM_Unbinned_likelihood(peri_phases)
        loc = Utils.range_shift(results.x[0],-np.pi, np.pi, 0, 1)
        kappas= results.x[1]
        null = MLE.null_likelihood(peri_phases)
        TS = 2*(null - results.fun)
        if kappa == 0:
            kappa0_kappa.append(kappas)
            loc0_kappa.append(loc)
            TS0_kappa.append(TS)
        elif kappa == 0.25:
            kappa25_kappa.append(kappas)
            loc25_kappa.append(loc)
            TS25_kappa.append(TS)
        elif kappa == 0.5:
            kappa5_kappa.append(kappas)
            loc5_kappa.append(loc)
            TS5_kappa.append(TS)
        elif kappa == 1:
            kappa1_kappa.append(kappas)
            loc1_kappa.append(loc)
            TS1_kappa.append(TS)
        elif kappa == 2:
            kappa2_kappa.append(kappas)
            loc2_kappa.append(loc)
            TS2_kappa.append(TS)
        elif kappa == 4:
            kappa4_kappa.append(kappas)
            loc4_kappa.append(loc)
            TS4_kappa.append(TS)
        elif kappa == 8:
            kappa8_kappa.append(kappas)
            loc8_kappa.append(loc)
            TS8_kappa.append(TS)
    print(kappa_list)
    for index, kappas in enumerate([0,0.25,0.5,1,2,4,8]):
        data_frame[str(star) + "Kappa_" + str(kappas)] = kappa_list[index]
        data_frame[str(star) + "Center_" + str(kappas)] = loc_list[index]
        data_frame[str(star) + "TS_{VM}_" + str(kappas)] = TS_list[index]
    data_frame.to_csv(output_dir, index = False)
