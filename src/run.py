# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 17:36:39 2025

@author: Nate Whitsett
"""
import pandas as pd
import numpy as np
import ardor.Utils.Utils as U
import ardor.Flares.Flare as Flare
import ardor.Statistical_Tests.MLE as MLE
import ardor.Plotting.Orbital_Flare_Hist as Plot
import ardor.Statistical_Tests.K_Tests as K
import os
from matplotlib import pyplot as plt
# data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_News.csv")
# catalog = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/TOI_List.csv")
data = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares.csv")
catalog = pd.read_csv("C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/Reference_Lists/Alfven_Catalog.csv")
host_list = []            
count = 0
host_name = []
host_KU = []
host_TS = []
lengths = []
dlogz_list = []
period_list = []
beta_list = []
e_list = []
spot_epoch= 2457000 + 1579.598146 + 1.3
## Stat Significant
# candidates = ['AUMic','WASP-34', 'Kepler-1972', 'HD260655', 'Gl49', 'TOI-2145']
# label_host = ['AU Mic','WASP-34', 'Kepler-1972', 'HD 260655', 'Gl 49', 'TOI-2145']
## Beta Significant
candidates = ['HD163607']
label_host = ['HD 163607']
## TOI Hosts
# candidates = [2329.01, 1232.01,723.01,1520.01,1887.01,5051.01,1187.01]

# a_list = [0.066,0.118,0.023,0.069,0.029,0.042,0.041]
# emax = [0.65,.82,0.5,0.71,0.7,0.3,0.63]
# rmax = [0.023,0.021,0.011,0.019,0.01,0.029,0.015]
letter = ['b','b','b','b','b','b']
for indx, host in enumerate(candidates):
    # if host == 'Gl49':
    flares = data.loc[(data['Host_ID'] == host), 'Flare_Epoch'].to_numpy()
    planets = catalog.loc[(catalog['Host_ID'] == host), 'pl_letter']
    # planets = ['b']
    dlogz = np.mean(data.loc[(data['Host_ID'] == host), 'dlogZ'])
    obs_time = np.mean(data.loc[(data['Host_ID'] == host), 'Obs_Time'])
    for idx, planet in enumerate(planets):

        ## Collect astrophysical parameters to use to analyze
        t_peri1, period, t_peri2,t_tran, a, e, omega, omega_u,omega_l, Rst, der_st_per, selection, temp, alfven, alfven_u, alfven_l, cat_strot = U.df_return(catalog, host, 
                                                  ['t_peri', 'pl_orbper', 'pl_orbtper', 'pl_tranmid', 'pl_orbsmax', 'pl_orbeccen', 'pl_orblper', 'omega_u', 
                                                    'omega_l','st_rad', 'der_st_per', 'selection', 'st_teff', 'Alfven_Rad', '0.13', '-0.135', 'st_rotp'], planet = planet)
        # period, t_tran, Rst, der_st_per,selection, temp, st_rad = U.df_return(catalog, host, 
        #                                           ['pl_orbper', 'pl_tranmid','st_rad', 'der_st_per','selection','st_teff','st_rad'], planet = planet)
        if np.isnan(der_st_per) == False:
            st_per = der_st_per
        elif np.isnan(der_st_per) == False:
            st_per = der_st_per
        else:
            st_per = np.nan  
        if period > st_per:
            synodic_period = (np.abs(1/st_per - 1/period))**-1
        elif period < st_per:
            synodic_period = (np.abs(1/period - 1/st_per))**-1
        try:
            if np.isnan(t_peri1[0]) == False and np.isnan(t_peri2[0]) == False:
                t_peri = t_peri2[0]
            elif np.isnan(t_peri1[0]) == False and np.isnan(t_peri2[0]) == True:
                t_peri = t_peri1[0]
            elif np.isnan(t_peri1[0]) == True and np.isnan(t_peri2[0]) == False:
                t_peri = t_peri2[0]
            else:
                t_peri = t_peri1[0]
        except:
            continue
        phases = Flare.phase_folder(flares, period[0], t_peri)
        # st_phasesses = Flare.phase_folder(flares, st_per, 2457000)
        # syn_phases = Flare.phase_folder(np.array(flares), synodic_period, 2457000)        
        result = MLE.VM_Unbinned_likelihood(phases)
        # result_st = MLE.VM_Unbinned_likelihood(st_phases)
        # result_syn = MLE.VM_Unbinned_likelihood(syn_phases)
        tests = K.K_Tests(flares, period,[t_peri], output_message=False)
        # test_st = K.K_Tests(flares, st_per,[2457000], output_message=False)
        # test_syn = K.K_Tests(flares, synodic_period,[2457000], output_message=False)
        if len(flares) >= 4 and planet == 'b':
            # beta = MLE.SPI_Metric(period[0], e[0], a[0], np.array(phases), omega[0], Rst[0], t_tran=t_tran[0], transit = False, 
            #                       los_factor = 0, plot=False)
            print(t_tran)
            if np.isnan(t_tran) == True:
                transit = False
            elif np.isnan(t_tran) == False:
                transit = True
            plt.hist(phases)
            # Plot.Polar_Flare_Plot(phases + 0.5, e=e[0], a=a[0], omega_p=omega[0]+180, star_rad = Rst[0], 
            #                           peri_err = [omega_l[0],omega_u[0]], title=rf'{label_host[indx]} {letter[indx]}', star_teff=temp[0],
            #                           alfven=[alfven[0], alfven_u[0],alfven_l[0]], color='inferno', transit=transit, save_dir=os.getcwd()
            #                         , alfven_bool=True, scale = True, legend = False)
                # beta = MLE.SPI_Metric(period[0], e[0], a[0], np.array(phases), omega[0], Rst[0], t_tran=t_tran[0], transit = False, 
                #                       los_factor = 0, plot=True)
                # print(host, result_st[2], len(flares), st_per)
                # tests = K.K_Tests(flares, period,[t_peri], output_message=False)
                # beta_list.append(beta)
                # e_list.append(e[0])
                # host_name.append(host)
                # plt.show()
                # if beta > 0.2:
                #     MLE.SPI_Metric(period[0], e[0], a[0], np.array(phases), omega[0], Rst[0], t_tran=t_tran[0], transit = False, 
                #                           los_factor = 0, plot=True)
                #     print(host)
                #     plt.show()
                #     Plot.Polar_Flare_Plot(phases + 0.5, e=e[0], a=a[0], omega_p=omega[0]+180, star_rad = Rst[0], 
                #                              peri_err = [omega_l[0],omega_u[0]], title=rf'{host} {planet}', star_teff=temp[0],
                #                              alfven=[alfven[0], alfven_u[0],alfven_l[0]], color='inferno', transit=True, save_dir=os.getcwd()
                #                             , alfven_bool=True)
                #     plt.show()
                # if synodic_period < 30 and len(flares) > 3 and np.isnan(synodic_period) == False:
                
                # tests = K.K_Tests(flares, period, [t_tran[0]], output_message=False)
                # Plot.Polar_Flare_Plot(phases + 0.5, e=0, a=a_list[indx], omega_p=0, star_rad = Rst[0], title=f'TOI-{candidates[indx]}', star_teff=temp[0], color='inferno', transit=True, 
                #                       save_dir=os.getcwd(), TOI=True, dynamic_lim=True, emax=emax[indx], alfven_bool=False, roche=[a_list[indx]*(1-emax[indx])], legend=True, scale=True)
                # Plot.Polar_Flare_Plot(phases + 0.5, e=e[0], a=a[0], omega_p=omega[0]+180, star_rad = Rst[0], 
                #                          peri_err = [omega_l[0],omega_u[0]], title=f'{label_host[indx]} {planet}', star_teff=temp[0],
                #                          alfven=[alfven[0], alfven_u[0],alfven_l[0]], color='inferno', transit=True, save_dir=None
                #                         )
                # plt.show()
                # Plot.St_Polar_Flare_Plot(syn_phases, a=0.05, star_rad=Rst[0], star_teff=temp[0], title='TOI-' + str(host)[:-3])
                # host_TS.append(result_syn[2])
                # host_name.append(host)
                # lengths.append(len(syn_phases))
                # dlogz_list.append(dlogz)
                # period_list.append(synodic_period)
                # host_KU.append(tests[1])
                # count += 1
                # plt.show()
# fig, ax = plt.subplots(figsize=(4,4))
# ax.hist(beta_list, histtype='step', color='black', bins=15)
# ax.set_xlabel(r"$\beta_{{SPI}}$", fontsize=14)
# ax.set_ylabel("Count", fontsize=14)
# # plt.title(r'$\beta_{{SPI}}$ for Hosts with $e>0$')
# plt.savefig(os.path.join(os.getcwd(), 'beta_samp.png'), dpi=400)
# plt.show()
# X = np.column_stack((host_name,beta_list, e_list))
# X = np.column_stack((host_name,period_list, lengths, host_KU, host_TS, dlogz_list))
# np.savetxt('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/KS_KU_Tests/K_Tests/Syn_Rot/Beta_SPI.csv', X, delimiter=',', fmt='%s')

