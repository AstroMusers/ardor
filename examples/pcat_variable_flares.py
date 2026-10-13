#!/usr/bin/env python3
"""Demonstrate variable-count flare fitting with PCAT on simulated photometry.

The example injects 2 to 5 FRED flares into a 0.3-day light curve with a
constant baseline of 1 relative-flux unit and Gaussian measurement noise of
0.004 relative-flux units. These are illustrative synthetic data, not an
observed stellar light curve.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pcat.time_series import evaluate_flare_profile

from ardor.Flares.pcat_flares import posterior_flux_quantiles, tier3_pcat


EXAMPLE_PATH = Path(__file__).resolve().parent
DURATION_DAYS = 0.3  # [day]
CADENCE_MINUTES = 2.0  # [minute]
BASELINE_FLUX = 1.0  # [relative flux]
FLUX_UNCERTAINTY = 0.004  # [relative flux]
REFERENCE_BJD = 2460000.0  # [day]


def simulate_light_curve(seed: int = 0):
    """Return an Ardor Tier 2-style CSV payload and its injected flare times."""

    rng = np.random.default_rng(seed)
    cadence_days = CADENCE_MINUTES / (24.0 * 60.0)  # [day]
    time_days = np.arange(0.5 * cadence_days, DURATION_DAYS, cadence_days)
    number_flares = int(rng.integers(2, 6))
    peak_times = np.sort(rng.uniform(0.04, DURATION_DAYS - 0.04, number_flares))  # [day]
    amplitudes = rng.uniform(0.015, 0.05, number_flares)  # [relative flux]
    rise_times = rng.uniform(1.0, 4.0, number_flares) / (24.0 * 60.0)  # [day]
    decay_times = rng.uniform(5.0, 25.0, number_flares) / (24.0 * 60.0)  # [day]
    injected_profiles = evaluate_flare_profile(
        time_days,
        "flarfred",
        amplitudes,
        peak_times,
        rise_time=rise_times,
        decay_time=decay_times,
    )
    noiseless_flux = BASELINE_FLUX + np.sum(injected_profiles, axis=1)
    uncertainty = np.full(time_days.size, FLUX_UNCERTAINTY)
    observed_flux = rng.normal(noiseless_flux, uncertainty)
    return time_days, observed_flux, uncertainty, peak_times


def save_fit_figures(fit, injected_peak_times, run_name: str) -> list[Path]:
    """Save the full-baseline model fit and flare-count posterior."""

    output_directory = EXAMPLE_PATH / "visuals"
    output_directory.mkdir(parents=True, exist_ok=True)
    quantiles = posterior_flux_quantiles(fit)
    hours = (fit.time_days - fit.time_days[0]) * 24.0  # [hour]

    figure, axes = plt.subplots(
        2, 1, figsize=(8.0, 4.8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    axes[0].errorbar(
        hours,
        fit.flux,
        yerr=fit.uncertainty,
        fmt=".",
        ms=2.5,
        lw=0.5,
        color="#454545",
        label="Simulated Ardor light curve",
    )
    axes[0].fill_between(
        hours, quantiles[0], quantiles[2], color="#007360", alpha=0.22,
        label="PCAT 68% posterior interval",
    )
    axes[0].plot(hours, quantiles[1], color="#007360", lw=1.2,
                 label="PCAT posterior median")
    for index, peak_time in enumerate(injected_peak_times):
        axes[0].axvline(
            (peak_time - fit.time_days[0]) * 24.0,
            color="#A51417",
            ls=":",
            lw=0.8,
            label="Injected flare peaks" if index == 0 else None,
        )
    axes[0].set_ylabel("Relative flux")
    axes[0].legend(loc="upper right", fontsize=8)
    residual = (fit.flux - quantiles[1]) / fit.uncertainty
    axes[1].axhline(0.0, color="#777777", lw=0.8)
    axes[1].plot(hours, residual, ".", color="#454545", ms=2.5)
    axes[1].set_ylabel(r"Residual [$\sigma$]")
    axes[1].set_xlabel("Time [hour]")
    figure.tight_layout()
    fit_path = output_directory / f"pcat_variable_flare_fit_{run_name}.png"
    print(f"Writing to {fit_path}...")
    figure.savefig(fit_path, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    flare_counts = np.asarray(fit.posterior.listpostnumbelem, dtype=int).ravel()
    values, counts = np.unique(flare_counts, return_counts=True)
    figure, axis = plt.subplots(figsize=(4.0, 3.0))
    axis.bar(values, counts / counts.sum(), color="#A51417")
    axis.axvline(injected_peak_times.size, color="#007360", ls="--", lw=1.2,
                 label="Injected flare count")
    axis.set_xlabel("Number of flares")
    axis.set_ylabel("Posterior probability")
    axis.legend(loc="upper right", fontsize=8)
    figure.tight_layout()
    count_path = output_directory / f"pcat_variable_flare_count_{run_name}.png"
    print(f"Writing to {count_path}...")
    figure.savefig(count_path, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return [fit_path, count_path]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numbswep", type=int, default=20_000)
    parser.add_argument("--max-flares", type=int, default=10)
    parser.add_argument("--numbproc", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--smoke", action="store_true", help="Run a short pipeline check.")
    arguments = parser.parse_args()

    time_days, flux, uncertainty, injected_peak_times = simulate_light_curve(arguments.seed)
    run_name = datetime.now().strftime("pcat_flares_%Y%m%d_%H%M%S")
    input_directory = EXAMPLE_PATH / "data" / "inpt" / run_name
    input_directory.mkdir(parents=True, exist_ok=True)
    tier2_path = input_directory / "Flare_1.csv"
    print(f"Writing to {tier2_path}...")
    np.savetxt(
        tier2_path,
        np.column_stack((time_days, flux, uncertainty)),
        delimiter=",",
        header=str(REFERENCE_BJD),
    )

    fits = tier3_pcat(
        input_directory,
        EXAMPLE_PATH / "data" / "outp" / run_name,
        host_name="Synthetic_Star",
        numbswep=1200 if arguments.smoke else arguments.numbswep,
        max_flares=arguments.max_flares,
        numbproc=1 if arguments.smoke else arguments.numbproc,
    )
    save_fit_figures(fits[0], injected_peak_times, run_name)
    sampled_counts = np.asarray(fits[0].posterior.listpostnumbelem).astype(int).ravel()
    print(
        f"Injected {injected_peak_times.size} flares; "
        f"PCAT sampled flare counts {sorted(set(sampled_counts))}."
    )


if __name__ == "__main__":
    main()