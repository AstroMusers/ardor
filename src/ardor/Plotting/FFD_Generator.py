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
from pathlib import Path
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
def FFD_Generator(data_dir, output_csv=None):
    data = pd.read_csv(data_dir)
    hosts = data['Host_ID']
    energy_lists = []
    x_lists = []
    error_lists = []
    params_list = []
    param_err_list_a = []
    param_err_list_b = []
    host_list = []
    n_flares_list = []
    hosts = set(hosts)
    for host in hosts:
        energy_count_list = []
        err_list = []
        obs_time = np.array(data.loc[data['Host_ID'] == host, 'Obs_Time'])[0]/(24*60)
        if obs_time == 0 or np.isnan(obs_time) == True:
            continue
        energies = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'Energy']).dropna())
        dlogz = np.array(pd.Series(data.loc[(data['Host_ID'] == host), 'dlogZ']))
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
            host_list.append(host)
            n_flares_list.append(len(energies))

    coeff_catalog = pd.DataFrame({
        "Host_ID": host_list,
        "N_flares": n_flares_list,
        "alpha": [params[0] for params in params_list],
        "alpha_err": param_err_list_a,
        "beta": [params[1] for params in params_list],
        "beta_err": param_err_list_b,
    })
    if not coeff_catalog.empty:
        coeff_catalog = coeff_catalog.sort_values("Host_ID").reset_index(drop=True)
    if output_csv is not None:
        Path(output_csv).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)
        coeff_catalog.to_csv(output_csv, index=False)

    return energy_lists, x_lists, error_lists, params_list, param_err_list_a, param_err_list_b, coeff_catalog
