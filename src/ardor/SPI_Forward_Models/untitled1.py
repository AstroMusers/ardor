# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 14:17:04 2025

@author: Nate Whitsett
"""

from ardor.Statistical_Tests.K_Tests import K_Tests
from ardor.Flares.Flare import phase_folder
import pandas as pd
from ardor.Utils.Utils import df_return, add_list_to_csv
from ardor.Statistical_Tests.MLE import SPI_Metric
import re
import numpy as np
from matplotlib import pyplot as plt
G_data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/Output/VM/G_Type/Kappa_Sim_G.csv")
M_data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/Output/VM/M_Type/Kappa_Sim_M.csv")
set_hosts = set(M_data['Host_ID'])
es = [0.2, 0.5, 0.8]
for e in es:
    KUs = []
    betas = []
    for host in set_hosts:
        nums = re.findall(r"\d+\.?\d*", host)
        values = [float(n) if "." in n else int(n) for n in nums]
        name, count, str_e = values
        flares = np.array(df_return(M_data, host, ['Flare_Epoch']))
        if float(str_e) == e and len(flares) > 3:
            results = K_Tests(flares, [5], [2457000], KS=False, AD=False, sampling=False, output_message=False)
            beta = SPI_Metric(5, e, 0.05, phase_folder(flares,5,2457001.75), 180, 1,transit=False, los_factor=20, plot=True, s)
            betas.append(beta)
            print(host, e)
            if results[0] is not None:
                KUs.append(results[0])
            elif results[0] is None:
                KUs.append(1)
    print(KUs)
    add_list_to_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_Sim.csv", 
                    f'beta_e_{e}', betas, output_path="C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_Sim.csv")
    add_list_to_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_Sim.csv", 
                    f'KU_e_{e}', KUs, output_path="C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Simulations/M_Sim.csv")