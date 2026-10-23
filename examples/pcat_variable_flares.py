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
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pcat.time_series import evaluate_flare_profile
from tdpy.verbosity import print

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
    print(f"Writing to visuals/{fit_path.name}...")
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
    print(f"Writing to visuals/{count_path.name}...")
    figure.savefig(count_path, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return [fit_path, count_path]


def save_fit_bokeh_figures(fit, injected_peak_times, run_name: str) -> list[Path]:
    """Save interactive views of the same variable-count PCAT posterior."""
    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource, HoverTool, Span
    from bokeh.plotting import figure
    from tdpy.plotting import save_bokeh_html

    output_directory = EXAMPLE_PATH / "visuals"
    output_directory.mkdir(parents=True, exist_ok=True)
    quantiles = posterior_flux_quantiles(fit)
    time_hours = (fit.time_days - fit.time_days[0]) * 24.0  # [hour]
    median_flux = quantiles[1]
    source = ColumnDataSource({
        "time_hours": time_hours,
        "flux": fit.flux,
        "uncertainty": fit.uncertainty,
        "lower_flux": quantiles[0],
        "median_flux": median_flux,
        "upper_flux": quantiles[2],
        "residual_sigma": (fit.flux - median_flux) / fit.uncertainty,
    })
    fit_plot = figure(
        title="PCAT flare fit and posterior interval",
        width=850,
        height=420,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_axis_label="Time [hour]",
        y_axis_label="Relative flux",
    )
    fit_plot.varea(
        x="time_hours",
        y1="lower_flux",
        y2="upper_flux",
        source=source,
        fill_color="#007360",
        fill_alpha=0.22,
        legend_label="PCAT 68% posterior interval",
    )
    fit_plot.line(
        "time_hours",
        "median_flux",
        source=source,
        color="#007360",
        line_width=2,
        legend_label="PCAT posterior median",
    )
    fit_plot.scatter(
        "time_hours",
        "flux",
        source=source,
        color="#454545",
        size=5,
        legend_label="Synthetic observed flux",
    )
    fit_plot.add_tools(HoverTool(tooltips=[
        ("Time [hour]", "@time_hours{0.00}"),
        ("Flux", "@flux{0.0000}"),
        ("Uncertainty", "@uncertainty{0.0000}"),
        ("PCAT median", "@median_flux{0.0000}"),
    ]))
    for peak_time in injected_peak_times:
        fit_plot.add_layout(Span(
            location=(peak_time - fit.time_days[0]) * 24.0,
            dimension="height",
            line_color="#A51417",
            line_dash="dotted",
            line_width=1.5,
        ))
    fit_plot.legend.location = "top_right"
    fit_plot.legend.click_policy = "hide"
    fit_plot.grid.visible = False

    residual_plot = figure(
        title="Posterior residuals",
        width=850,
        height=230,
        sizing_mode="stretch_width",
        x_range=fit_plot.x_range,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_axis_label="Time [hour]",
        y_axis_label="Residual [sigma]",
    )
    residual_plot.scatter(
        "time_hours",
        "residual_sigma",
        source=source,
        color="#454545",
        size=5,
    )
    residual_plot.add_layout(Span(
        location=0.0,
        dimension="width",
        line_color="#777777",
        line_width=1,
    ))
    residual_plot.grid.visible = False
    fit_path = output_directory / f"pcat_variable_flare_fit_{run_name}.html"
    fit_relative_path = Path(os.path.relpath(fit_path, Path.cwd()))
    save_bokeh_html(
        column(fit_plot, residual_plot, sizing_mode="stretch_width"),
        fit_relative_path,
        title="Ardor PCAT variable-flare fit",
    )

    flare_counts = np.asarray(fit.posterior.listpostnumbelem, dtype=int).ravel()
    values, counts = np.unique(flare_counts, return_counts=True)
    count_source = ColumnDataSource({
        "flare_count": values,
        "probability": counts / counts.sum(),
    })
    count_plot = figure(
        title="PCAT flare-count posterior",
        width=700,
        height=360,
        x_range=(float(values.min()) - 0.5, float(values.max()) + 0.5),
        y_range=(0.0, 1.0),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_axis_label="Number of flares",
        y_axis_label="Posterior probability",
    )
    count_renderer = count_plot.vbar(
        x="flare_count",
        top="probability",
        width=0.8,
        source=count_source,
        color="#A51417",
    )
    count_plot.add_tools(HoverTool(
        renderers=[count_renderer],
        tooltips=[("Flares", "@flare_count"), ("Probability", "@probability{0.000}")],
    ))
    count_plot.add_layout(Span(
        location=int(np.size(injected_peak_times)),
        dimension="height",
        line_color="#007360",
        line_dash="dashed",
        line_width=2,
    ))
    count_plot.grid.visible = False
    count_path = output_directory / f"pcat_variable_flare_count_{run_name}.html"
    count_relative_path = Path(os.path.relpath(count_path, Path.cwd()))
    save_bokeh_html(
        count_plot,
        count_relative_path,
        title="Ardor PCAT sampled flare count",
    )
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
    try:
        save_fit_bokeh_figures(fits[0], injected_peak_times, run_name)
    except ImportError:
        print("Install Bokeh to create the interactive example figures.")
    sampled_counts = np.asarray(fits[0].posterior.listpostnumbelem).astype(int).ravel()
    print(
        f"Injected {injected_peak_times.size} flares; "
        f"PCAT sampled flare counts {sorted(set(sampled_counts))}."
    )


if __name__ == "__main__":
    main()