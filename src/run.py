# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 17:36:39 2025

@author: whitsett.n
"""
import pandas as pd
import numpy as np
import ardor.Utils.Utils as U
import ardor.Flares.Flare as Flare
import ardor.Statistical_Tests.MLE as MLE
import ardor.Plotting.Orbital_Flare_Hist as Plot
import ardor.Statistical_Tests.K_Tests as K
from matplotlib import pyplot as plt
# data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_News.csv")
# catalog = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOI_List.csv")
data = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")
catalog = pd.read_csv("C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Alfven_Catalog.csv")
host_list = []            
count = 0
host_name = []
host_KU = []
host_TS = []
lengths = []
dlogz_list = []
period_list = []
beta_list = []
spot_epoch= 2457000 + 1579.598146 + 1.3
for host in set(data['Host_ID']):
    # if host not in host_list:
    if host == "TOI-1062":
        flares = data.loc[(data['Host_ID'] == host), 'Flare_Epoch']
        planets = catalog.loc[(catalog['Host_ID'] == host), 'pl_letter']
        dlogz = np.mean(data.loc[(data['Host_ID'] == host), 'dlogZ'])
        for idx, planet in enumerate(planets):     
            period, t_peri,t_tran, a, e, omega, omega_u,omega_l, Rst, st_per, selection, temp, alfven, alfven_u, alfven_l = U.df_return(catalog, host, 
                                                      ['pl_orbper', 't_peri', 'pl_tranmid', 'pl_orbsmax', 'pl_orbeccen', 'omega', 'omega_u', 
                                                       'omega_l','st_rad', 'der_st_per', 'selection', 'st_teff', 'Alfven_Rad', '0.13', '-0.135'], planet = planet)
            # period, t_peri, Rst, st_per,selection, temp, st_rad = U.df_return(catalog, host, 
            #                                           ['pl_orbper', 'pl_tranmid','st_rad', 'der_st_per','selection','st_teff','st_rad'], planet = planet)

            phases = Flare.phase_folder(flares, period, t_peri)
            phases_st = Flare.phase_folder(flares, st_per, 2457000)
            result = MLE.VM_Unbinned_likelihood(phases)
            result_st = MLE.VM_Unbinned_likelihood(phases_st)
            # periods = [period[idx], st_per[0]]
            if result is not None and idx == 0 and e[0] > 0 and len(flares) > 3:
                if result[2] > 0:
                    # if np.isnan(t_tran[0]) == True:
                    #     beta = MLE.SPI_Metric(period[0], e[0], a[0], phases, omega[0], Rst[0], los_factor=2, transit = False)
                    # if np.isnan(t_tran[0]) == False:
                    # for shifts in np.linspace(0,1, num=10):
                    beta = MLE.SPI_Metric(period[0], e[0], a[0], phases, omega[0], Rst[0], t_tran = t_tran[0], los_factor=5, transit = True)
                    print(host, beta, result[1], result[2])
                    beta_list.append(beta)
                if result[2] > 3:
                #     print(host)
                    # Plot.Polar_Flare_Plot(phases, e=e[0], a=a[0], omega_p=omega[0], star_rad = Rst[0], peri_err = [omega_l[0],omega_u[0]], title=host,
                    #                       alfven=[alfven, alfven_u, alfven_l])
                    test = K.K_Tests(flares, period, [t_peri])
                    test_st = K.K_Tests(flares, st_per, [2457000])
                    TS = U.round_to_sig_figs(result[2], 2)
                    pku = U.round_to_sig_figs(test[1], 2)
                    TS_st = U.round_to_sig_figs(result_st[2], 2)
                    pku_st = U.round_to_sig_figs(test_st[1], 2)
                    period_disp = U.round_to_sig_figs(period[0], 3)
                    print(host, planet, pku, TS)
                    # st_per[0] = U.round_to_sig_figs(st_per[0], 4)
                    # bottom_text_st = [fr"$T_{{p}} =${st_per[0]} d", fr"$\sqrt{{TS}}=${TS_st}", fr"$p_{{KU}}=${pku_st}"]
                    # bottom_text = [fr"$T_{{p}} =${period_disp} d", fr"$\sqrt{{TS}}=${TS}", fr"$p_{{KU}}=${pku}"]
                    # Plot.St_Polar_Flare_Plot(phases_st, title=host, a = a[0], star_rad=Rst[0], star_teff = temp[0], bottom_text_list=bottom_text_st)
                    # plt.show()
                    # Plot.Polar_Flare_Plot(phases, e=e[0], a=a[0], omega_p=omega[0], title=host, star_rad=Rst[0], peri_err=[omega_u[0],omega_l[0]], 
                    #                       star_teff=temp[0], bottom_text_list=bottom_text, alfven=None, legend=True)
                    # plt.show()
                    host_TS.append(result[2])
                    host_name.append(host)
                    lengths.append(len(phases))
                    dlogz_list.append(dlogz)
                    period_list.append(st_per)
                    count += 1
        host_list.append(host)
# fig, ax = plt.subplots(figsize=(4,4))
# ax.hist(beta_list)
# ax.set_xlabel(r"$\beta_{{SPI}}$", fontsize=14)
# ax.set_ylabel("Count", fontsize=14)
# plt.title(r'$\beta_{{SPI}}$ for Hosts with $e>0$')
# plt.show()
# X = np.column_stack((host_name,period_list, lengths, host_KU, host_TS, dlogz_list))
# np.savetxt('C:/Users/whitsett.n/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/KS_KU_Tests/K_Tests/St_Rot/TOI_Rot_K_Cat.csv', X, delimiter=',', fmt='%s')

