# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 13:26:00 2024

@author: natha
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.colors as colors
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from matplotlib import rcParams
from collections import namedtuple
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
import math
for fontpath in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    print(fontpath.lower())
    if 'lmroman10-regular'.lower() in fontpath.lower():
        path = fontpath
for fontpath in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    if 'lmroman10-italic'.lower() in fontpath.lower():
        italicpath = fontpath
        print(italicpath)
font = FontProperties(fname='/ugrad/whitsett.n/fonts/latin-modern-roman/lmroman10-regular.otf')
italicfont = FontProperties(fname='/ugrad/whitsett.n/fonts/latin-modern-roman/lmroman10-italic.otf')
rcParams["mathtext.fontset"] = "cm"
def linear(x, alpha, beta):
    return alpha*x + beta
def round_to_sig_figs(x, sig_figs):
    """Round a number to a specified number of significant figures."""
    if x == 0:
        return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (sig_figs - 1))
def FFD_Generator(data_dir, dict_params = None, third_variables=['name', 'N', 'logZ', 'Teff', 'Age', 'Met', 'Obs_Time']):
    data = pd.read_csv(data_dir)
    hosts = data['Host_ID']
    energy_lists = []
    x_lists = []
    error_lists = []
    params_list = []
    param_err_list_a = []
    param_err_list_b = []
    meta_data = []
    hosts = set(hosts)
    FFD_Tuple = namedtuple('FFD', third_variables)
    for host in hosts:
        energy_count_list = []
        err_list = []
        obs_time = np.array(data.loc[data['Host_ID'] == host, 'Obs_Time'])[0]/(24*60)
        if obs_time == 0 or np.isnan(obs_time) == True:
            continue
        energies = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'Energy']).dropna())
        dlogz = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'log(Z)']))
        Teff = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'Teff']))[0]
        age = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'st_age']))[0]
        met = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'st_met']))[0]
        count = 0
        for index, values in enumerate(dlogz):
            if values < 5:
                if len(energies) == 0:
                    continue
                energies = np.delete(energies, count)
                dlogz = np.delete(dlogz, count)
                count -= 1
            count += 1
            
        energies[::-1].sort()
        energies = energies[~np.isnan(energies)]
        energies = np.log10(energies)
        if len(energies) > 10:
            a = np.linspace(energies[0], energies[-1], num = 6)
            for bins in a:
                b = sum(x >= bins for x in energies)
                energy_count_list.append(np.log10(b/(obs_time)))
                err_list.append(1/b)
            coeff, pcov = curve_fit(linear, a, energy_count_list, sigma = err_list, absolute_sigma=True)
            perr = np.sqrt(np.diagonal(pcov))
            a_err, b_err = perr
            param_err_list_a.append(a_err)
            param_err_list_b.append(b_err)
            energy_lists.append(energy_count_list)
            x_lists.append(a)
            error_lists.append(err_list)
            params_list.append((coeff[0], coeff[1]))
            meta_data.append(FFD_Tuple(host,len(energies), dlogz, Teff,age, met,obs_time))
    return energy_lists, x_lists, error_lists, params_list, param_err_list_a,param_err_list_b, meta_data












############### EXOPLANET HOST FFD GENERATOR###################################
# rcParams["mathtext.default"] = 'rm'
# data = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares_New.csv')
# a, b,c,d,a_err, b_err,e= FFD_Generator('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/Exoplanet_Hosts/All_Exoplanet_MCMC_Flares_New.csv', dict_params=True)
# fig, ax = plt.subplots(2, 3, layout='constrained')
# fig.set_size_inches(7, 3.5)
# x = np.linspace(30,37)
# teff = []
# met = []
# age = []
# for index, target in enumerate(e):
#     teff.append(e[index].Teff)
#     met.append(e[index].Met)
#     age.append(e[index].Age)
# met = np.array(met)
# age = np.array(age)
    
# norm = matplotlib.colors.Normalize(vmin=min(teff), vmax=max(teff), clip=True)
# mapper = cm.ScalarMappable(norm=norm, cmap='autumn')
# teff_color = np.array([(mapper.to_rgba(v)) for v in teff])
# for index, data in enumerate(range(len(a))):
#     if teff[index] < 4000:
#         ax[0][0].errorbar(b[index], a[index], yerr=c[index], fmt='o', alpha = 0.75, markersize = 3, capsize = 3,  c=teff_color[index])
#         ax[0][0].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, color = teff_color[index], linewidth=0.5)
#     else:
#         ax[1][0].errorbar(b[index], a[index], yerr=c[index], fmt='o', c=teff_color[index],alpha = 0.5, markersize = 3, capsize=3)
#         ax[1][0].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, c=teff_color[index],linewidth=0.5)
# cb = plt.colorbar(mapper, ax=ax[0][0], orientation='horizontal', location='top', aspect = 18, pad = 0.03, shrink=0.75, label = r'$\mathit{T}_{eff} \,(\mathrm{K})$')
# cb.ax.tick_params(labelsize=8, labelfontfamily='serif')
# norm = matplotlib.colors.Normalize(vmin=min(met[~np.isnan(met)]), vmax=max(met[~np.isnan(met)]))
# mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
# met_color = np.array([(mapper.to_rgba(v)) for v in met])
# for index, data in enumerate(range(len(met))):
#     if np.isnan(met[index]) == True:
#         continue
#     if met[index] < 0:
#         ax[0][1].errorbar(b[index], a[index], yerr=c[index], fmt='o', alpha = 0.75, markersize = 3, capsize = 3,  c=met_color[index])
#         ax[0][1].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, color = met_color[index], linewidth=0.5)
#     else:
#         ax[1][1].errorbar(b[index], a[index], yerr=c[index], fmt='o', alpha = 0.75, markersize = 3, capsize = 3,  c=met_color[index])
#         ax[1][1].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, color = met_color[index], linewidth=0.5)
# cb = plt.colorbar(mapper, ax=ax[0][1], orientation='horizontal', location='top', aspect = 18, pad = 0.03, shrink=0.75, label = r'$\mathrm{Metallicity\,[Fe/H]}$')
# cb.ax.tick_params(labelsize=8, labelfontfamily='serif')
# norm = matplotlib.colors.Normalize(vmin=0, vmax=15, clip=True)
# mapper = cm.ScalarMappable(norm=norm, cmap='plasma')
# teff_color = np.array([(mapper.to_rgba(v)) for v in age])
# for index, data in enumerate(range(len(age))):
#     if np.isnan(age[index]) == True:
#         continue
#     if age[index] < 2:
#         ax[0][2].errorbar(b[index], a[index], yerr=c[index], fmt='o', alpha = 0.75, markersize = 3, capsize = 3,  c=teff_color[index])
#         ax[0][2].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, color = teff_color[index], linewidth=0.5)
#     else:
#         ax[1][2].errorbar(b[index], a[index], yerr=c[index], fmt='o', alpha = 0.75, markersize = 3, capsize = 3,  c=teff_color[index])
#         ax[1][2].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, color = teff_color[index], linewidth=0.5)
# cb = plt.colorbar(mapper, ax=ax[0][2], orientation='horizontal', location='top', aspect = 18, pad = 0.03, shrink=0.75, label = r'$\mathit{t}_{age} (\mathrm{Gyr})$')
# cb.ax.tick_params(labelsize=8, labelfontfamily='serif')
# ax[0][0].set_xlim(30.5, 37)
# ax[0][1].set_xlim(30.5, 37)
# ax[0][2].set_xlim(30.5, 37)
# ax[1][0].set_xlim(30.5, 37)
# ax[1][1].set_xlim(30.5, 37)
# ax[1][2].set_xlim(30.5, 37)
# ax[0][0].set_ylim(-1.8, 5.3)
# ax[0][1].set_ylim(-1.8, 5.3)
# ax[0][2].set_ylim(-1.8, 5.3)
# ax[1][0].set_ylim(-1.8, 5.3)
# ax[1][1].set_ylim(-1.8, 5.3)
# ax[1][2].set_ylim(-1.8, 5.3)

# ax[0][0].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# ax[0][1].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# ax[1][0].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# ax[1][1].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# ax[1][2].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# fig.text(-.02, 0.5, r'$\log_{10} (cum. \, \mathit{N}_{flare} \geq \mathit{E}_{\mathrm{bol}}) (\mathrm{d}^{-1}) $', ha='center', va='center', rotation='vertical', font = font, fontsize = 12)
# fig.text(0.5, -0.03, r'$log_{10} (\mathit{E}_{\mathrm{bol}}) [erg]$', ha='center', va='center', font = font, fontsize = 12)
# # ax[0][1].get_yaxis().set_ticks([])
# # ax[1][1].get_yaxis().set_ticks([])
# # ax[1][2].get_yaxis().set_ticks([])
# # ax[0][2].get_yaxis().set_ticks([])
# # ax[0][2].get_xaxis().set_ticks([])
# # ax[0][0].get_xaxis().set_ticks([])
# # ax[0][1].get_xaxis().set_ticks([])
# fig.get_layout_engine().set(w_pad=0 / 72, h_pad=0 / 72, hspace=0,
#                             wspace=0)

# panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
# for i, ax in enumerate(ax.flat):
#     ax.annotate(panel_labels[i], xy=(0.87, 0.95), xycoords='axes fraction', 
#                 fontsize=12, va='top', font=font)
#     if i == 0 or i== 1 or i == 2:
#         for tick in ax.xaxis.get_major_ticks():
#             tick.tick1line.set_visible(False)
#             tick.tick2line.set_visible(False)
#             tick.label1.set_visible(False)
#             tick.label2.set_visible(False)
#     if i == 1 or i== 2 or i == 4 or i == 5:
#         for tick in ax.yaxis.get_major_ticks():
#             tick.tick1line.set_visible(False)
#             tick.tick2line.set_visible(False)
#             tick.label1.set_visible(False)
#             tick.label2.set_visible(False)
#     ax.grid(True)
# plt.savefig('Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Publication_2024/FFD_Hosts.png', dpi=500, bbox_inches='tight')
# plt.show()





######################## TOI FFD GENERATOR ####################################
# rcParams["mathtext.default"] = 'rm'
# data = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_New.csv')
# a, b,c,d, a_err, b_err, e = FFD_Generator('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Induced_Flares/Flare_Catalogs/All TOIs/All_TOI_MCMC_Flares_New.csv', third_variables=['name','N','logZ', 'Teff', 'Obs_Time'])
# names = []
# a_param = []
# b_param = []
# N = []
# dlogz_list = []
# T_list = []
# age = []
# for index, name in enumerate(e):
#     names.append(e[index].name)
#     a_param.append("$" + str(round_to_sig_figs(d[index][0],3)) + "\pm" + str(round_to_sig_figs(a_err[index], 1)) + "$")
#     b_param.append("$" + str(round_to_sig_figs(d[index][1],3)) + "\pm" + str(round_to_sig_figs(b_err[index], 1)) + "$")
#     N.append(e[index].N)
#     T_list.append(e[index].Teff)
#     dlogz_list.append(round_to_sig_figs(np.median(e[index].logZ),2))
#     age.append(e[index].Age)
# output = pd.DataFrame({"Name": names, "N": N, "a": a_param,"b": b_param, "Teff": T_list, "dlogZ": dlogz_list})
# output.to_csv("C:/Users/natha/Desktop/FFD_Host_Data.csv", index = False)
# fig, ax = plt.subplots(2, 1, layout='constrained')
# fig.set_size_inches(3.5, 3.5)
# x = np.linspace(30,37)
# teff = []
# met = []
# age = []
# for index, target in enumerate(e):
#     teff.append(e[index].Teff)
    
# norm = matplotlib.colors.Normalize(vmin=min(teff), vmax=max(teff), clip=True)
# mapper = cm.ScalarMappable(norm=norm, cmap='autumn')
# teff_color = np.array([(mapper.to_rgba(v)) for v in teff])
# for index, data in enumerate(range(len(a))):
#     if teff[index] < 4000:
#         ax[0].errorbar(b[index], a[index], yerr=c[index], fmt='o', alpha = 0.75, markersize = 3, capsize = 3,  c=teff_color[index])
#         ax[0].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, color = teff_color[index], linewidth=0.5)
#     else:
#         ax[1].errorbar(b[index], a[index], yerr=c[index], fmt='o', c=teff_color[index],alpha = 0.5, markersize = 3, capsize=3)
#         ax[1].plot(x, linear(x, d[index][0], d[index][1]), linestyle='-', alpha=0.75, c=teff_color[index],linewidth=0.5)
# cb = plt.colorbar(mapper, ax=ax[0], orientation='horizontal', location='top', aspect = 18, pad = 0.03, shrink=0.75, label = '$\mathit{T}_{eff} \,(\mathrm{K})$')
# cb.ax.tick_params(labelsize=8, labelfontfamily='serif')
# ax[0].set_xlim(32, 36)
# ax[1].set_xlim(32, 36)
# ax[0].set_ylim(-1.8, 3.9)
# ax[1].set_ylim(-1.8, 3.9)

# ax[0].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# ax[1].tick_params(axis='both', which='major', labelfontfamily='serif', labelsize =9)
# fig.text(-.03, 0.5, r'$\log_{10} (cum.\,\mathit{N}_{flare} \geq \mathit{E}_{bol}) (\mathrm{d}^{-1}) $', ha='center', va='center', rotation='vertical', font = italicfont, fontsize = 12)
# fig.text(0.5, -0.04, r'$\log_{10} (\mathit{E}_{\mathrm{bol}}) [erg]$', ha='center', va='center', font = italicfont, fontsize = 12)
# fig.get_layout_engine().set(w_pad=0 / 72, h_pad=0 / 72, hspace=0,
#                             wspace=0)

# panel_labels = ['(a)', '(b)']
# for i, ax in enumerate(ax.flat):
#     ax.annotate(panel_labels[i], xy=(0.9, 0.95), xycoords='axes fraction', 
#                 fontsize=12, va='top', font=font)
#     if i == 0:
#         for tick in ax.xaxis.get_major_ticks():
#             tick.tick1line.set_visible(False)
#             tick.label1.set_visible(False)
#     ax.grid(True)

# plt.savefig('Users/natha/OneDrive - Washington University in St. Louis/Desktop/Research/Publication Documents/Publication_2024/FFD_TOI.png', dpi=500, bbox_inches='tight')
# plt.show()