# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:59:11 2023

@author: Nathan
"""
import numpy as np
from scipy.integrate import simpson


def planck_law(lambda0, T):
    '''
    Returns the spectral flux density at a specified wavelength
    and blackbody temperature.

    Parameters
    ----------
    lambda0 : float or array-like
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
    lambda0 = np.asarray(lambda0, dtype=float)
    return ((h * c**2) / (lambda0**5)) / np.expm1((h * c) / (lambda0 * k * T))


def planck_integrator(lambda_min, lambda_max, Teff, num=50):
    """Integrate Planck spectral flux density over a wavelength range."""
    lambda_range = np.linspace(lambda_min, lambda_max, num=num)
    spectral_density = planck_law(lambda_range, Teff)
    return simpson(spectral_density, x=lambda_range)
