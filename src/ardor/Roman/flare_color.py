import numpy as np
import pandas as pd
import csv
from astropy import units as u
import synphot as syn
import stsynphot as stsyn
import stpsf
from expecto import get_spectrum
from ardor.Flares.aflare import aflare1, aflare
from ardor.Pipeline_Completeness.Injection_Recovery import amp_log_normal, FWHM_uniform
from scipy.interpolate import make_interp_spline
from pathlib import Path
import os

DEFAULT_ROMAN_RATIO_DIR = Path("/data2/whitsett.n/dustmaps/ratio_models")
DEFAULT_VEGA_PATH = Path("/data2/whitsett.n/Spectra/grp/redcat/trds/calspec/alpha_lyr_stis_011.fits")

ROMAN_WFI = stpsf.WFI()
ROMAN_FILTERS = {
    "F087": ROMAN_WFI._get_synphot_bandpass("F087"),
    "F146": ROMAN_WFI._get_synphot_bandpass("F146"),
    "F213": ROMAN_WFI._get_synphot_bandpass("F213"),
}

def get_data_file(filename: str) -> Path:
    """Return the full path to a ratio-grid data file.

    Resolution order:
    1) ``ROMAN_RATIO`` environment variable.
    2) ``DEFAULT_ROMAN_RATIO_DIR`` fallback.
    """
    env_dir = os.environ.get("ROMAN_RATIO")
    search_dirs = []
    if env_dir:
        search_dirs.append(Path(env_dir).expanduser())
    search_dirs.append(DEFAULT_ROMAN_RATIO_DIR)

    for data_dir in search_dirs:
        candidate = data_dir / filename
        if candidate.exists():
            return candidate

    searched = ", ".join(str(path) for path in search_dirs)
    raise FileNotFoundError(
        f"Could not find '{filename}'. Searched directories: {searched}"
    )

def ensure_synphot_source_spectrum(spectrum):
    if isinstance(spectrum, syn.SourceSpectrum):
        return spectrum

    if hasattr(spectrum, "spectral_axis") and hasattr(spectrum, "flux"):
        return syn.SourceSpectrum(
            syn.models.Empirical1D,
            points=spectrum.spectral_axis.to(u.AA),
            lookup_table=spectrum.flux,
        )

    raise TypeError("Unsupported spectrum type returned by expecto.get_spectrum")


def build_blackbody_spectrum(teff_k):
    return syn.SourceSpectrum(syn.models.BlackBody1D, temperature=teff_k * u.K)


def build_phoenix_spectrum(teff_k, logg, mh, cache=True):
    phoenix_spec = get_spectrum(teff_k, logg, mh, cache=cache)
    return ensure_synphot_source_spectrum(phoenix_spec)


def _get_vega_path() -> Path:
    """Return the Vega spectrum path.

    Resolution order:
    1) ``ROMAN_VEGA`` environment variable.
    2) ``DEFAULT_VEGA_PATH`` fallback.
    """
    env_path = os.environ.get("ROMAN_VEGA")
    if env_path:
        p = Path(env_path).expanduser()
        if p.exists():
            return p
    if DEFAULT_VEGA_PATH.exists():
        return DEFAULT_VEGA_PATH
    raise FileNotFoundError(
        f"Could not find Vega spectrum. Set the ROMAN_VEGA environment variable or "
        f"place the file at {DEFAULT_VEGA_PATH}"
    )


def observe_spectrum_flux(
    spectrum,
    filter_names=["F146","F087", "F213"],
    telescope_diameter=2.4 * u.m,
    exposure_time=66 * u.s,
    normalize_mag=19.0 * u.ABmag,
    normalize_band_name="gaia, grp",
    vega_path=None,
    ebv=0.0,
    extinction_model="mwdense",
):
    observed_fluxes = dict()
    observed_flux_errors = dict()
    for filter_name in filter_names:
        if filter_name not in ROMAN_FILTERS:
            raise ValueError(f"Unsupported filter '{filter_name}'. Valid options: {list(ROMAN_FILTERS)}")

    normalize_band = stsyn.band(normalize_band_name)
    ebv = float(ebv)
    if ebv < 0.0:
        raise ValueError("ebv must be >= 0")

    extinction_curve = None
    if ebv > 0.0:
        reddening_law = syn.ReddeningLaw.from_extinction_model(extinction_model)
        extinction_curve = reddening_law.extinction_curve(ebv)

    # Keep the normalization band as a native obs-mode bandpass to avoid
    # metadata issues (e.g., Missing OBSMODE) in synphot normalization.
    # Extinction is still applied when computing each observed filter flux.

    resolved_vega_path = Path(vega_path) if vega_path is not None else _get_vega_path()
    vegaspec = syn.SourceSpectrum.from_file(str(resolved_vega_path))
    normalized_spectrum = spectrum.normalize(normalize_mag, band=normalize_band, vegaspec=vegaspec)

    telescope_area = np.pi * (telescope_diameter / 2) ** 2
    for filter_name in filter_names:
        bandpass = ROMAN_FILTERS[filter_name]
        if extinction_curve is not None:
            bandpass = bandpass * extinction_curve
        observation = syn.Observation(normalized_spectrum, bandpass)
        observed_flux = observation.countrate(area=telescope_area).to(u.ct / u.s)

        flux_error = np.sqrt((observed_flux * exposure_time).to_value(u.ct)) / exposure_time.to_value(u.s)
        flux_error = flux_error * (u.ct / u.s)
        observed_fluxes[filter_name] = observed_flux
        observed_flux_errors[filter_name] = flux_error

    return observed_fluxes, observed_flux_errors



def create_syn_lc(
    flux_arrays,
    flux_errs,
    cadence=12,
    duration=3,
    spot_amplitude=0.0,
    spot_period_days=None,
    spot_phase=None, normalize = False
):
    """
    Parameters
    ----------
    spot_amplitude : float
        Semi-amplitude of the sinusoidal star-spot modulation as a fraction
        of the baseline flux (e.g. 0.01 = 1% peak-to-peak half-amplitude).
        Set to 0 to disable.
    spot_period_days : float or None
        Rotation period of the star in days. If None, a random period between
        1 and 30 days is drawn from the RNG.
    spot_phase : float or None
        Initial phase offset in radians. If None, a random phase in [0, 2π)
        is drawn from the RNG.
    """
    lcs = []
    errs = []
    time = np.arange(0, duration, cadence / (60 * 24))
    for idx, flux in enumerate(flux_arrays):
        flux_val = flux.value if hasattr(flux, "value") else float(flux)
        flux_err_val = flux_errs[idx].value if hasattr(flux_errs[idx], "value") else float(flux_errs[idx])
        flux_array = np.random.normal(loc=flux_val, scale=flux_err_val, size=len(time))

        if spot_amplitude > 0.0:
            period_min = (np.random.uniform(1.0, 30.0) if spot_period_days is None else float(spot_period_days)) * 24.0 * 60.0
            phase = np.random.uniform(0.0, 2.0 * np.pi) if spot_phase is None else float(spot_phase)
            spot_signal = spot_amplitude * flux_val * np.sin(2.0 * np.pi * time / period_min + phase)
            flux_array += spot_signal
        lcs.append(flux_array)
        errs.append(flux_err_val * np.ones_like(flux_array))
    if normalize:
        lcs_norm = [lc / np.median(lc) for lc in lcs]
        errs_norm = [err / np.median(lc) for lc, err in zip(lcs, errs)]
        return time, lcs_norm, errs_norm
    return time, lcs, errs


def inject_random_flares(
    time,
    flux_arrays,
    n_flares=3,
    upsample_factor=20,
    ratio_params={"Teff": 9000, "EBV": 0.1},
    base_fluxes=[2000, 200, 200],
    behavior="random",
    t_peaks_days=None,
    fwhms_days=None,
    amplitudes_frac=None,
):
    """
    Inject one or more synthetic flare events into one or multiple light curves.

    Parameters
    ----------
    time : array-like
        Time array in days. All input flux arrays must match this length.
    flux_arrays : array-like or sequence of array-like
        One flux array or a list/tuple of flux arrays sampled on ``time``.
        Each channel is injected independently using the same flare event times
        and widths, then scaled by ``flux_ratios``.
    n_flares : int, optional
        Number of flares to inject when ``behavior='random'``. Ignored when
        ``behavior='set'`` (in that case it is inferred from provided arrays).
    upsample_factor : int, optional
        Temporal upsampling factor used to build a smooth flare model.
    ratio_params : list of [temperature, fraction], optional
        Parameters defining the flare amplitude ratios for each channel.
        Each element should be a list or tuple with two values:
        [temperature, E(B-V)]. If a single element is provided, it is
        used for all channels.
    behavior : {'random', 'set'}
        'random' – draw flare parameters stochastically (default).
        'set'    – use the explicitly provided t_peaks_days, fwhms_days, and
                   amplitudes_frac; n_flares is inferred from their length.
    t_peaks_days : array-like, optional
        Peak times in days for each flare. Required when behavior='set'.
    fwhms_days : array-like, optional
        Full-width at half-maximum in days for each flare.
        Required when behavior='set'.
    amplitudes_frac : array-like, optional
        Flare amplitude as a fraction of the per-channel baseline flux for
        each flare. Required when behavior='set'.

    Returns
    -------
    injected_flux_arrays : list of ndarray
        Input flux arrays with flare signal added on the native cadence.
    flare_events : list of dict
        Event definitions with keys ``t_peak_day``, ``fwhm_day``, and
        ``amplitude_frac``.
    flare_models_on_cadence : list of ndarray
        Flare-only models sampled on the original ``time`` grid, one per
        channel after applying ``flux_ratios``.
    smooth_time : ndarray
        Upsampled time grid used for smooth flare models.
    smooth_flare_models : list of ndarray
        Flare-only models evaluated on ``smooth_time``, one per channel after
        applying ``flux_ratios``.
    """
    ratios = {}
    for band in ("F087", "F213"):
        table = pd.read_csv(get_data_file(f"{band}.csv"), index_col=0)
        ebv_vals = np.array(table.columns, dtype=float)
        teff_vals = np.array(table.index, dtype=float)
        teff_idx = (np.abs(teff_vals - ratio_params["Teff"])).argmin()
        ebv_idx = (np.abs(ebv_vals - ratio_params["EBV"])).argmin()
        ratios[band] = table.iloc[teff_idx, ebv_idx]
    print(ratios)
    F087_ratio = ratios["F087"]
    F213_ratio = ratios["F213"]
    behavior = behavior.lower()
    if behavior not in ("random", "set"):
        raise ValueError("behavior must be 'random' or 'set'")

    if behavior == "set":
        if t_peaks_days is None or fwhms_days is None or amplitudes_frac is None:
            raise ValueError(
                "behavior='set' requires t_peaks_days, fwhms_days, and amplitudes_frac"
            )
        t_peaks_days = np.asarray(t_peaks_days, dtype=float)
        fwhms_days = np.asarray(fwhms_days, dtype=float)
        amplitudes_frac = np.asarray(amplitudes_frac, dtype=float)
        lengths = {len(t_peaks_days), len(fwhms_days), len(amplitudes_frac)}
        if len(lengths) != 1:
            raise ValueError(
                "t_peaks_days, fwhms_days, and amplitudes_frac must all have the same length"
            )
        n_flares = len(t_peaks_days)

    if n_flares < 1:
        time = np.asarray(time, dtype=float)
        if not isinstance(flux_arrays, (list, tuple)):
            flux_arrays = [flux_arrays]
        empty_model = np.zeros_like(time)
        empty_models = [empty_model.copy() for _ in flux_arrays]
        return list(flux_arrays), [], empty_models, time, empty_models
    time = np.asarray(time, dtype=float)

    if not isinstance(flux_arrays, (list, tuple)):
        flux_arrays = [flux_arrays]

    flare_events = []
    if behavior == "random":
        for _ in range(n_flares):
            flare_events.append(
                {
                    "t_peak_day": np.random.uniform(np.min(time), np.max(time)),
                    "fwhm_day": FWHM_uniform(),
                    "amplitude_frac": amp_log_normal(),
                }
            )
    else:  # behavior == "set"
        for i in range(n_flares):
            flare_events.append(
                {
                    "t_peak_day": float(t_peaks_days[i]),
                    "fwhm_day": float(fwhms_days[i]),
                    "amplitude_frac": float(amplitudes_frac[i]),
                }
            )

    injected_flux_arrays = []
    flare_models_on_cadence = []
    smooth_flare_models = []
    smooth_time = np.linspace(np.min(time), np.max(time), len(time) * int(upsample_factor))

    for idx, flux_array in enumerate(flux_arrays):
        if hasattr(flux_array, "unit"):
            flux_array = np.full(len(time), flux_array.to_value(u.ct / u.s), dtype=float)
        else:
            flux_array = np.asarray(flux_array, dtype=float)
            if flux_array.ndim == 0:
                flux_array = np.full(len(time), float(flux_array), dtype=float)
        if flux_array.size != len(time):
            raise ValueError("Each flux array must have the same length as time")

        baseline = float(np.median(flux_array))
        flare_params = []
        for event in flare_events:
            if idx == 0:  # F087
                flux_ratio = 1
                amplitude = baseline * event["amplitude_frac"]
            if idx == 1:  # F087
                flux_ratio = F087_ratio
                amplitude = compute_new_flare_amp(amplitude, flux_ratio, base_fluxes[0], base_fluxes[1])
            if idx == 2:
                flux_ratio = F213_ratio
                amplitude = compute_new_flare_amp(amplitude, flux_ratio, base_fluxes[0], base_fluxes[2])
            flare_params.append(
                {
                    "t_peak_day": event["t_peak_day"],
                    "fwhm_day": event["fwhm_day"],
                    "amplitude_ct_per_s": amplitude,
                }
            )

        smooth_time, smooth_flare_model, flare_model_on_cadence = build_smooth_flare_model(
            time,
            flare_params,
            upsample_factor=upsample_factor,
        )


        injected_flux = flux_array + flare_model_on_cadence
        injected_flux_arrays.append(injected_flux)
        flare_models_on_cadence.append(flare_model_on_cadence)
        smooth_flare_models.append(smooth_flare_model)

    return injected_flux_arrays, flare_events, flare_models_on_cadence, smooth_time, smooth_flare_models


def build_smooth_flare_model(time, flare_params, upsample_factor=20):
    time = np.asarray(time, dtype=float)
    if upsample_factor < 1:
        raise ValueError("upsample_factor must be >= 1")

    smooth_time = np.linspace(np.min(time), np.max(time), len(time) * int(upsample_factor))
    smooth_flare_model = np.zeros_like(smooth_time)
    flare_model_on_cadence = np.zeros_like(time)

    for params in flare_params:
        smooth_flare_model += aflare1(
            smooth_time,
            params["t_peak_day"],
            params["fwhm_day"],
            params["amplitude_ct_per_s"],
        )
        flare_model_on_cadence += aflare1(
            time,
            params["t_peak_day"],
            params["fwhm_day"],
            params["amplitude_ct_per_s"],
        )

    return smooth_time, smooth_flare_model, flare_model_on_cadence


def export_blackbody_ratio_grid_csv(
    csv_path="blackbody_wfi_ratio_grid_3000_15000K_step10.csv",
    temp_min=3000,
    temp_max=13660,
    temp_step=10,
    extinct_step=0.1, band = "F087"
):
    temperatures = np.arange(temp_min, temp_max + temp_step, temp_step, dtype=float)
    extinct = np.arange(0, 10 + extinct_step, extinct_step, dtype=float)
    reference_band = stsyn.band("roman, wfi, f146")
    reference_mag = 19.0 * u.ABmag
    telescope_area = np.pi * (2.4 * u.m / 2) ** 2
    reddening_law = syn.ReddeningLaw.from_extinction_model("mwdense")

    # Create 2D grid: rows are temperatures, columns are E(B-V) values
    ratio_grid = np.zeros((len(temperatures), len(extinct)))
    
    for i, teff in enumerate(temperatures):
        blackbody_spec = build_blackbody_spectrum(teff)
        normalized_spec = blackbody_spec.normalize(reference_mag, band=reference_band, force='extrap')
        
        for j, ebv in enumerate(extinct):
            # Apply extinction to the bandpass rather than the source spectrum.
            # Multiplying the analytical BlackBody1D source by an empirical extinction
            # curve converts it to an empirical spectrum with a bounded wavelength range,
            # which can become disjoint from filter bandpasses at high temperatures.
            # Keeping the source analytical and folding extinction into the bandpass avoids
            # the DisjointError entirely.
            extinction_curve = reddening_law.extinction_curve(ebv)
            extincted_band = ROMAN_FILTERS[band] * extinction_curve
            extincted_f146 = ROMAN_FILTERS["F146"] * extinction_curve

            band_rate_ext = syn.Observation(
                normalized_spec, extincted_band, force='extrap'
            ).countrate(area=telescope_area).to_value(u.ct / u.s)
            f146_rate_ext = syn.Observation(
                normalized_spec, extincted_f146, force='extrap'
            ).countrate(area=telescope_area).to_value(u.ct / u.s)
            
            ratio_grid[i, j] = band_rate_ext / f146_rate_ext
        
        print(teff, band_rate_ext / f146_rate_ext)
    
    # Write 2D grid to CSV
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        # Header row: Teff_K followed by E(B-V) values
        header = ["Teff_K"] + [f"EBV_{ebv:.1f}" for ebv in extinct]
        writer.writerow(header)
        # Data rows: temperature followed by F087/F146 ratios for each E(B-V)
        for i, teff in enumerate(temperatures):
            row = [int(teff)] + ratio_grid[i, :].tolist()
            writer.writerow(row)

def create_color_obs(time, flux_f146, flux_f087, flux_f213, phase):
    """Slice three input flux arrays onto Roman-like alternating color cadence.

    At a fixed index cadence (``step=30``), samples are assigned to:
    - F087 at ``phase + n*step``
    - F213 at ``(phase+15 mod 30) + n*step``
    - F146 at all remaining indices

    Parameters
    ----------
    time : array-like
        Time array.
    flux_f146 : array-like
        Flux array for F146 sampled on ``time``.
    flux_f087 : array-like
        Flux array for F087 sampled on ``time``.
    flux_f213 : array-like
        Flux array for F213 sampled on ``time``.
    phase : int
        Starting index for F087 in [0, 29] (or < len(time) for short arrays).

    Returns
    -------
    time_out : dict
        ``{"F146": ..., "F087": ..., "F213": ...}`` sliced time arrays.
    fluxes_out : dict
        ``{"F146": ..., "F087": ..., "F213": ...}`` sliced flux arrays.
    idx_out : dict
        ``{"F146": ..., "F087": ..., "F213": ...}`` integer indices into
        the original full-cadence arrays.
    """
    step = 30
    time = np.asarray(time)
    flux_f146 = np.asarray(flux_f146)
    flux_f087 = np.asarray(flux_f087)
    flux_f213 = np.asarray(flux_f213)

    n = time.shape[0]
    if flux_f146.shape[0] != n or flux_f087.shape[0] != n or flux_f213.shape[0] != n:
        raise ValueError("time and all flux arrays must have the same length")

    start_idx_f087 = int(phase)
    if start_idx_f087 < 0 or start_idx_f087 >= min(step, n):
        raise ValueError("phase must be an integer index in [0, 29] (or < len(time) if shorter)")

    start_idx_f213 = (start_idx_f087 + 15) % step

    idx_f087 = np.arange(start_idx_f087, n, step, dtype=int)
    idx_f213 = np.arange(start_idx_f213, n, step, dtype=int)

    idx_f146 = np.setdiff1d(np.arange(n, dtype=int), np.union1d(idx_f087, idx_f213), assume_unique=False)

    time_out = {
        "F146": time[idx_f146],
        "F087": time[idx_f087],
        "F213": time[idx_f213],
    }
    fluxes_out = {
        "F146": flux_f146[idx_f146],
        "F087": flux_f087[idx_f087],
        "F213": flux_f213[idx_f213],
    }

    idx_out = {
        "F146": idx_f146,
        "F087": idx_f087,
        "F213": idx_f213,
    }

    return time_out, fluxes_out, idx_out


def compute_new_flare_amp(amp, ratio, F146_base, band_base):
    F146_excess = F146_base * (amp)
    band_amp = (F146_excess * ratio) / band_base
    return band_amp