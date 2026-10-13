'''
Hold the analytic flare model

Could also put other simple flare models in here, e.g. polynomial from Balona (2015)
'''

import numpy as np

try:
    from tdpy.exoplanet import davenport_flare_area, evaluate_davenport_flare
except ImportError:
    def evaluate_davenport_flare(time, peak_time, fwhm, amplitude):
        t = np.asarray(time, dtype=float)
        phase = (t - peak_time) / fwhm
        rise = (
            1.00000
            + 1.94053 * phase
            - 0.175084 * phase**2
            - 2.24588 * phase**3
            - 1.12498 * phase**4
        )
        decay = 0.689008 * np.exp(-1.60053 * phase) + 0.302963 * np.exp(-0.278318 * phase)
        shape = np.where(phase < -1.0, 0.0, np.where(phase <= 0.0, rise, decay))
        shape = np.clip(shape, 0.0, None)
        return shape * np.abs(amplitude)

    def davenport_flare_area(fwhm, amplitude):
        return np.abs(amplitude) * np.abs(fwhm) * 1.8268848560866364


def aflare(t, p):
    """
    This is the Analytic Flare Model from the flare-morphology paper.
    Reference Davenport et al. (2014) http://arxiv.org/abs/1411.3723

    Note: this model assumes the flux before the flare is zero centered

    Parameters
    ----------
    t : 1-d array
        The time array to evaluate the flare over
    p : 1-d array
        p == [tpeak, fwhm (units of time), amplitude (units of flux)] x N

    Returns
    -------
    flare : 1-d array
        The flux of the flare model evaluated at each time
    """
    Nflare = int(np.floor(len(p) / 3.0))
    flare = np.zeros_like(t, dtype=float)
    for i in range(Nflare):
        flare = flare + evaluate_davenport_flare(t, p[i * 3], p[i * 3 + 1], p[i * 3 + 2])
    return flare


def aflare1(t, tpeak, fwhm, ampl):
    '''
    The Analytic Flare Model evaluated for a single-peak (classical).
    Reference Davenport et al. (2014) http://arxiv.org/abs/1411.3723

    Parameters
    ----------
    t : 1-d array
        The time array to evaluate the flare over
    tpeak : float
        The time of the flare peak
    fwhm : float
        The "Full Width at Half Maximum", timescale of the flare
    ampl : float
        The amplitude of the flare

    Returns
    -------
    flare : 1-d array
        The flux of the flare model evaluated at each time
    '''
    return evaluate_davenport_flare(t, tpeak, fwhm, ampl)


def flare_total_area(fwhm, ampl):
    """Analytical integrated area for Davenport et al. (2014) flare model."""
    return davenport_flare_area(fwhm, ampl)


def exp_decay(x, a, b, c):
    """Exponential decay function for flare fitting: a * exp(-b * x) + c."""
    return a * np.exp(-b * x) + c
