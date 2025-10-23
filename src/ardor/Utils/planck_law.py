# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:59:11 2023

@author: Nathan
"""
import numpy as np
from scipy.integrate import simpson
from matplotlib import pyplot as plt
import pandas as pd
# import pyckles as pk
from scipy.integrate import simpson
def planck_law(lambda0, T):
    '''
    Returns the spectral flux density at a specified wavelength 
    and blackbody temperature.

    Parameters
    ----------
    lambda0 : float
        Wavelength (cm)
    T : float
        Blackbody temperature (K)

    Returns
    -------
    array
        Returns the spectral density at a specified blackbody temperature, T,
        at a specified wavelength, lambda0

    '''
    h = 6.626e-27
    k = 1.38e-16
    c = 2.997e10
    return ((h*c**2)/(lambda0**5))/(np.exp((h*c)/(lambda0*k*T))-1)

def planck_integrator(lambda_min, lambda_max, Teff):
    lambda_min = lambda_min
    lambda_max = lambda_max
    lambda_range = np.linspace(lambda_min, lambda_max, num=50)
    R_A = []
    for wavelengths in lambda_range:
       R_A.append(planck_law(wavelengths, Teff))
    integral2 = simpson(np.array(R_A), lambda_range)
    return integral2

def TESS_to_AB(TESS_mag, spectral_type, TESS_min=970, TESS_Max=1770, ULTRA_min=230, ULTRA_max=350):
    for subspec in spectra:
        ULT_flux = spec_lib[subspec].data["flux"][ULTRA_min:ULTRA_max]
        ULT_wave = spec_lib[subspec].data["wavelength"][ULTRA_min:ULTRA_max]
        TESS_flux = spec_lib[subspec].data["flux"][TESS_min:TESS_Max]
        TESS_wave = spec_lib[subspec].data["wavelength"][TESS_min:TESS_Max]
        spectra_dict[str(subspec)] = simpson(ULT_flux, ULT_wave)/simpson(TESS_flux,TESS_wave)

# TOI_List = pd.read_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/ULTRASAT/Exoplanets.csv')
# TOI_Teff = TOI_List['st_teff']
# TOI_TESS_Mag = TOI_List['sy_tmag']
# TOI_AB_Mag = []
# TOI_spec = TOI_List['st_spectype']
# spec_lib = pk.SpectralLibrary("pickles")
# spec_lib.available_spectra
# spectra_dict = dict()
# spectra = spec_lib.available_spectra.tolist()
# for subspec in spectra:
#     ULT_flux = spec_lib[subspec].data["flux"][230:350]
#     ULT_wave = spec_lib[subspec].data["wavelength"][230:350]
#     TESS_flux = spec_lib[subspec].data["flux"][900:1770]
#     TESS_wave = spec_lib[subspec].data["wavelength"][900:1770]
#     spectra_dict[str(subspec)] = simpson(ULT_flux, ULT_wave)/simpson(TESS_flux,TESS_wave)

# ratio = []
# count = 0
# for index in range(len(TOI_Teff)):
#     if pd.isna(TOI_spec[index]) == False:
#         ratio = TOI_spec[index].replace(' ','')
#         ratio = ratio[0] + '0V'
#         if len(ratio) < 3:
#             ratio = ratio + str('II')
#         try:
#             TOI_AB_Mag.append(-2.5*np.log10((2461*10**(-.4*TOI_TESS_Mag[index]))*spectra_dict[str(ratio)]) + 8.9)
#             count += 1
#             print(count)
#         except:
#             TOI_AB_Mag.append(-2.5*np.log10((2461*10**(-.4*TOI_TESS_Mag[index]))*(planck_integrator(0.23e-4, 0.29e-4, TOI_Teff[index])/planck_integrator(0.6e-4, 1e-4, TOI_Teff[index]))) + 8.9)
#     elif np.isnan(TOI_spec[index]) == True: 
#         TOI_AB_Mag.append(-2.5*np.log10((2461*10**(-.4*TOI_TESS_Mag[index]))*(planck_integrator(0.23e-4, 0.29e-4, TOI_Teff[index])/planck_integrator(0.6e-4, 1e-4, TOI_Teff[index]))) + 8.9)

# TOI_List['AB Mag'] = np.array(TOI_AB_Mag)
# TOI_List.to_csv('C:/Users/natha/OneDrive - Washington University in St. Louis/Desktop/ULTRASAT/Exoplanets.csv', index=False)
# TOI_AB_Mag = TOI_List['AB Mag']
# TOI_Dist = TOI_List['sy_dist']
# bins = np.linspace(5, 25, num=15)
# plt.title('Kepler Field TESS Mag vs ULTRASAT AB Mag')
# plt.hist(TOI_TESS_Mag, bins, alpha=0.5, label='TESS Magnitude')
# plt.hist(TOI_AB_Mag, bins, alpha=0.5, label='AB Magnitude')
# # plt.hist(ratio, bins)
# plt.legend()
# plt.xlabel('Magnitude')
# plt.ylabel('Count')
# plt.axvline(21.8, linestyle='dashed', linewidth=1)
# plt.show()

