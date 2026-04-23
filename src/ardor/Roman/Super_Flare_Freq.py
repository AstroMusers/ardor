#%%
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ardor.Flares import Flare
from ardor.Flares.aflare import aflare1
from astropy import units as u
from collections import namedtuple
from pathlib import Path
from time import perf_counter
import csv
try:
    from Roman.Nested_Sampling.flare_color import create_syn_lc, inject_random_flares
except ModuleNotFoundError:
    from Roman.Nested_Sampling.flare_color import create_syn_lc, inject_random_flares


def _mark_recovered(injection_dict, detected_indices, tolerance=50):
    for detected in detected_indices:
        for key in injection_dict.keys():
            if abs(key - detected) <= tolerance:
                injection_dict[key][4] = True


def _dict_to_rows(injection_dict):
    rows = []
    for _, value in injection_dict.items():
        rows.append(
            {
                "Amplitude": value[0],
                "FWHM": value[1],
                "Error": value[2],
                "Integral": value[5],
                "Injected?": value[3],
                "Accepted?": value[4],
            }
        )
    return rows


def Injection_Recovery(
    output_dir,
    catalog_path='/ugrad/whitsett.n/Roman/gaia_results.csv',
    n_iterations=5,
    n_stars=None,
    n_flares=15,
    sigma=3,
    chi_square_cutoff=10,
    cadence_min=12,
    duration_days=72,
    old=False,
    seed=42,
    teff_range=None,
):
    """Run injection-recovery using simulated Roman light curves.

    The baseline flux and uncertainty are read from a catalog with `Flux` and
    `Flux_err` columns (ct/s). Synthetic light curves are generated with
    `flare_color.create_syn_lc`, flare templates are injected with
    `flare_color.inject_random_flares`. The generated light curves are then
    converted to relative flux (baseline near 1.0), and recovery is evaluated with
    `Flare.flare_ID` (Tier 1) plus `Flare.tier2` (Tier 2).
    
    Parameters
    ----------
    teff_range : tuple of (float, float) or None
        If provided, filter the catalog to only include stars with effective temperature
        in the range [teff_min, teff_max]. Example: teff_range=(3000, 5000) for K/M dwarfs.
    """

    data = pd.read_csv(catalog_path)
    data = data[np.isfinite(data['Flux']) & np.isfinite(data['Flux_err'])].copy()
    data = data[data['Flux'] > 0]
    
    # Filter by teff range if provided
    if teff_range is not None:
        teff_min, teff_max = teff_range
        data = data[(data['teff_gspphot'] >= teff_min) & (data['teff_gspphot'] <= teff_max)]
        print(f'Filtered to {len(data)} stars in teff range {teff_range}')
    
    # Randomly sample n_stars from the filtered catalog
    if n_stars is not None:
        n_stars = int(min(n_stars, len(data)))
        data = data.sample(n=n_stars, random_state=seed)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if old:
        output_file_t1 = output_path / 'Injection_Recovery_T1_Old.csv'
        output_file_t2 = output_path / 'Injection_Recovery_T2_Old.csv'
    else:
        output_file_t1 = output_path / 'Injection_Recovery_T1.csv'
        output_file_t2 = output_path / 'Injection_Recovery_T2.csv'
    output_time = output_path / 'Time_Stat.csv'

    LightCurve = namedtuple(
        'LightCurve',
        ['time', 'flux', 'detrended_flux', 'error', 'fast_bool', 'obs_time', 'trend', 'sector'],
    )

    rng = np.random.default_rng(seed)
    tier_timing_rows = []
    t1_rows = []
    t2_rows = []
    lc_num = 0

    for _ in range(n_iterations):
        for _, row in data.iterrows():
            t0 = perf_counter()
            flux_q = float(row['Flux']) * (u.ct / u.s)
            flux_err_q = float(row['Flux_err']) * (u.ct / u.s)

            time_min, baseline_flux = create_syn_lc(
                flux_q,
                flux_err_q,
                cadence=cadence_min,
                duration=duration_days,
                rng=rng,
            )
            baseline_median = float(np.nanmedian(baseline_flux))
            if not np.isfinite(baseline_median) or baseline_median <= 0:
                continue
            baseline_flux = np.asarray(baseline_flux, dtype=float) / baseline_median
            relative_err = float(flux_err_q.to_value(u.ct / u.s) / baseline_median)

            injected_flux_arrays, injected_events, _, _, _ = inject_random_flares(
                time_min,
                [baseline_flux],
                n_flares=n_flares,
                rng=rng,
            )
            injected_flux = np.asarray(injected_flux_arrays[0], dtype=float)

            flare_inject_dict_t1 = {}
            flare_inject_dict_t2 = {}
            for event in injected_events:
                event_index = int(np.argmin(np.abs(time_min - event['t_peak_min'])))
                local_slice = slice(max(0, event_index - 50), min(len(baseline_flux), event_index + 50))
                local_err = float(np.std(baseline_flux[local_slice]))
                baseline = float(np.median(baseline_flux))
                amp = baseline * event['amplitude_frac']
                fwhm = float(event['fwhm_min']) / (24 * 60)
                integral = np.trapezoid(
                    aflare1(time_min, event['t_peak_min'], event['fwhm_min'], amp),
                    x=time_min,
                )
                flare_inject_dict_t1[event_index] = [amp, fwhm, local_err, True, False, integral]
                flare_inject_dict_t2[event_index] = [amp, fwhm, local_err, True, False, integral]
            t1 = perf_counter()

            flare_events_t1, lengths = Flare.flare_ID(
                injected_flux,
                sigma,
                fast=False,
                injection=True,
                old=old,
            )
            _mark_recovered(flare_inject_dict_t1, flare_events_t1)
            t2 = perf_counter()

            time_days = 2457000 + (time_min / (24 * 60))
            lc_error = np.full_like(injected_flux, relative_err, dtype=float)
            lc = LightCurve(
                time=time_days,
                flux=injected_flux,
                detrended_flux=injected_flux,
                error=lc_error,
                fast_bool=False,
                obs_time=float(duration_days),
                trend=np.ones_like(injected_flux),
                sector=-1,
            )

            _, tier2_results = Flare.tier2(
                lc,
                flare_events_t1,
                lengths,
                chi_square_cutoff=chi_square_cutoff,
                output_dir=None,
                host_name='Roman_Sim',
                csv=False,
                Sim=True,
                injection=True,
            )
            tier2_indices = tier2_results.get('event_index', [])
            _mark_recovered(flare_inject_dict_t2, tier2_indices)

            for key in list(flare_inject_dict_t2.keys()):
                if not flare_inject_dict_t1[key][4]:
                    del flare_inject_dict_t2[key]
            t3 = perf_counter()

            tier_timing_rows.append(
                {
                    'Tier0_tau_s': t1 - t0,
                    'Tier1_tau_s': t2 - t1,
                    'Tier2_tau_s': t3 - t2,
                }
            )
            t1_rows.extend(_dict_to_rows(flare_inject_dict_t1))
            t2_rows.extend(_dict_to_rows(flare_inject_dict_t2))

            lc_num += 1
            print(f'LC #: {lc_num}')

    pd.DataFrame(t1_rows).to_csv(output_file_t1, index=False)
    pd.DataFrame(t2_rows).to_csv(output_file_t2, index=False)
    pd.DataFrame(tier_timing_rows).to_csv(output_time, index=False)

    return {
        'n_lightcurves': lc_num,
        'n_t1_injections': len(t1_rows),
        'n_t2_survivors': len(t2_rows),
        't1_file': str(output_file_t1),
        't2_file': str(output_file_t2),
        'timing_file': str(output_time),
    }
def Injection_Recovery_Grid(data_dir, grid_dir, label = 'T1', energy = False):
    
    data = pd.read_csv(data_dir)
    amp = list(data['Amplitude'])
    FWHM = list(data['FWHM'])
    injection = list(data['Injected?'])
    Bool = list(data['Accepted?'])
    integral = list(data['Integral'])
    amp_bins = np.logspace(-3, 0, num=17)
    FWHM_bins = np.linspace(0.001388888, 0.041, num=17)
    int_bins = np.logspace(-5, -1, num=17)
    x = []
    y = []
    if energy == False:
        for i in range(len(amp_bins)-1):
            tmp = []
            for index in range(len(amp)):
                if amp[index] > amp_bins[i] and amp[index] < amp_bins[i+1]:
                    tmp.append([amp[index], FWHM[index], Bool[index]])
            x.append(tmp)
        total = 0
        for i in range(len(FWHM_bins)):
            tmp = []
            for cells in x:
                count = 0
                pos = 0
                for flares in cells:
                    if flares[1] > FWHM_bins[i] and flares[1] < FWHM_bins[i+1]:
                        count += 1
                        total += 1
                        if flares[2] == 1:
                            pos += 1
                if i == 16:
                    continue
                if count == 0:
                    tmp.append(np.nan)
                elif count != 0:
                    tmp.append(pos/count * 100)
            y.append(tmp)
        with open(grid_dir + '/Injection_Recovery_Grid' + label + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(y)
    elif energy == True:
        for i in range(len(int_bins)-1):
            tmp = []
            for index in range(len(integral)):
                if integral[index] > int_bins[i] and integral[index] < int_bins[i+1]:
                    tmp.append([integral[index], FWHM[index], Bool[index]])
            x.append(tmp)
        total = 0
        for i in range(len(FWHM_bins)):
            tmp = []
            for cells in x:
                count = 0
                pos = 0
                for flares in cells:
                    if flares[1] > FWHM_bins[i] and flares[1] < FWHM_bins[i+1]:
                        count += 1
                        total += 1
                        if flares[2] == 1:
                            pos += 1
                if i == 16:
                    continue
                if count == 0:
                    tmp.append(np.nan)
                elif count != 0:
                    tmp.append(pos/count * 100)
            y.append(tmp)

        with open(grid_dir + '/Injection_Recovery_Grid_Energy' + label + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(y)



if __name__ == '__main__':
    data = pd.read_csv('/ugrad/whitsett.n/Roman/Flare_Catalogs/Gunter_alpha.csv', skiprows=[1,4])
    catalog = pd.read_csv('/ugrad/whitsett.n/Roman/gaia_results.csv')
    temperature_range1 = [4200, 7000]
    temperature_range2 = [0,4200]
    alphas1 = data.loc[(data['Teff'] >= temperature_range1[0]) & (data['Teff'] <= temperature_range1[1]), 'aFFD'].to_list()
    betas1 = data.loc[(data['Teff'] >= temperature_range1[0]) & (data['Teff'] <= temperature_range1[1]), 'bFFD'].to_list()
    alphas2 = data.loc[(data['Teff'] >= temperature_range2[0]) & (data['Teff'] <= temperature_range2[1]), 'aFFD'].to_list()
    betas2 = data.loc[(data['Teff'] >= temperature_range2[0]) & (data['Teff'] <= temperature_range2[1]), 'bFFD'].to_list()
    FGK = []
    M = []
    for idx, alpha in enumerate(alphas1):
        if np.isnan(alpha) == True:
            continue
        flares = alphas1[idx] * 34 + betas1[idx]
        FGK.append(10**flares)
    for idx, alpha in enumerate(alphas2):
        if np.isnan(alpha) == True:
            continue
        flares = alphas2[idx] * 33 + betas2[idx]
        M.append(10**flares)
    
    gaia_FGK = catalog.loc[(catalog['teff_gspphot'] >= temperature_range1[0]) & (catalog['teff_gspphot'] <= temperature_range1[1]), 'teff_gspphot'].to_list()
    gaia_M = catalog.loc[(catalog['teff_gspphot'] >= temperature_range2[0]) & (catalog['teff_gspphot'] <= temperature_range2[1]), 'teff_gspphot'].to_list()
    total = 0
    FGK.sort()
    FGK.pop(-1)
    total_M = 0
    M.sort()
    print(M)
    for probs in range(int(len(gaia_FGK)*0.05)):
        FGK_probs = np.random.choice(FGK)
        happened = np.random.rand() < FGK_probs/48
        if happened == True:
            total += 1
    for probs in range(int(len(gaia_M)*0.3)):
        M_probs = np.random.choice(M)
        happened_M = np.random.rand() < M_probs/48
        if happened_M == True:
            total_M += 1
    print('FGK flares per observation: ' + str(total))
    print('M flares per observation: ' + str(total_M))
