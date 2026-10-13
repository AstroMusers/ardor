"""Sampler-independent flare, baseline, and auxiliary-ratio models."""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import CubicSpline, RegularGridInterpolator

from ardor.Flares.aflare import aflare1
from ardor.Roman.flare_color import get_data_file

RATIO_GRID_FILES = {"F087": "F087.csv", "F213": "F213.csv"}


def normalize_ratio_filter(ratio_filter):
    if ratio_filter is None or pd.isna(ratio_filter):
        return "F087"
    normalized = str(ratio_filter).strip().upper()
    if normalized not in RATIO_GRID_FILES:
        allowed = ", ".join(sorted(RATIO_GRID_FILES))
        raise ValueError(f"Unknown ratio_filter '{ratio_filter}'. Expected one of: {allowed}")
    return normalized


def resolve_ratio_grid_path(ratio_filter="F087", ratio_grid_csv=None):
    if ratio_grid_csv is not None:
        return Path(ratio_grid_csv)
    return Path(get_data_file(RATIO_GRID_FILES[normalize_ratio_filter(ratio_filter)]))


def infer_ratio_filter_from_path(ratio_grid_path):
    filename = Path(ratio_grid_path).name.upper()
    for ratio_filter, mapped_name in RATIO_GRID_FILES.items():
        if filename == mapped_name.upper():
            return ratio_filter
    return "CUSTOM"


def load_ratio_grid(grid_csv_path):
    """Load legacy or matrix Teff/E(B-V) CSV grids and build an interpolator."""
    dataframe = pd.read_csv(grid_csv_path)
    if "Teff_K" in dataframe.columns:
        teff_values = dataframe["Teff_K"].to_numpy(dtype=float)
        ebv_columns = [column for column in dataframe.columns if column.startswith("EBV_")]
        if not ebv_columns:
            raise ValueError(f"Ratio grid must include EBV_* columns: {grid_csv_path}")
        ebv_values = np.array([float(column.replace("EBV_", "")) for column in ebv_columns])
        ratio_values = dataframe[ebv_columns].to_numpy(dtype=float)
    else:
        dataframe = pd.read_csv(grid_csv_path, index_col=0)
        try:
            teff_values = dataframe.index.to_numpy(dtype=float)
            ebv_values = np.array(dataframe.columns, dtype=float)
        except ValueError as error:
            raise ValueError(
                "Ratio grid must use either a 'Teff_K' column with 'EBV_*' columns "
                "or a numeric Teff index with numeric EBV headers: "
                f"{grid_csv_path}"
            ) from error
        ratio_values = dataframe.to_numpy(dtype=float)

    interpolator = RegularGridInterpolator(
        (teff_values, ebv_values), ratio_values, bounds_error=False, fill_value=np.nan
    )
    return {
        "path": str(grid_csv_path),
        "teff_values": teff_values,
        "ebv_values": ebv_values,
        "ratio_values": ratio_values,
        "interpolator": interpolator,
    }


def infer_ratio_grid_bounds(ratio_grid):
    return {
        "teff_min": float(np.min(ratio_grid["teff_values"])),
        "teff_max": float(np.max(ratio_grid["teff_values"])),
        "ebv_min": float(np.min(ratio_grid["ebv_values"])),
        "ebv_max": float(np.max(ratio_grid["ebv_values"])),
    }


def baseline_model(t, baseline_type="none", params=None, spline_knots=None):
    if baseline_type == "none":
        return np.zeros_like(t)
    if baseline_type == "sine":
        c0, amp, freq, phase = params
        t0 = 0.5 * (np.min(t) + np.max(t))
        return c0 + amp * np.sin(2.0 * np.pi * freq * (t - t0) + phase)
    if baseline_type == "spline":
        if spline_knots is None:
            spline_knots = np.linspace(np.min(t), np.max(t), len(params))
        return CubicSpline(spline_knots, np.asarray(params), bc_type="natural")(t)
    raise ValueError(f"Unknown baseline_type: {baseline_type}")


def interpolate_ratio(teff, ebv, ratio_grid):
    point = np.array([float(teff), float(ebv)], dtype=float)
    return float(ratio_grid["interpolator"](point))


def split_theta(theta, baseline_type="none", n_spline_knots=5, use_ratio_model=False):
    parts = {
        "tpeak": float(theta[0]),
        "fwhm": float(theta[1]),
        "ampl": float(theta[2]),
        "sigma": float(theta[3]),
    }
    index = 4
    if baseline_type == "sine":
        parts["baseline_params"] = np.asarray(theta[index:index + 4], dtype=float)
        index += 4
    elif baseline_type == "spline":
        parts["baseline_params"] = np.asarray(theta[index:index + n_spline_knots], dtype=float)
        index += n_spline_knots
    else:
        parts["baseline_params"] = np.array([], dtype=float)
    if use_ratio_model:
        parts["teff"] = float(theta[index])
        parts["ebv"] = float(theta[index + 1])
    return parts


def model_flux(theta, t, baseline_type="none", spline_knots=None, n_spline_knots=5, use_ratio_model=False):
    parts = split_theta(
        theta,
        baseline_type=baseline_type,
        n_spline_knots=n_spline_knots,
        use_ratio_model=use_ratio_model,
    )
    if baseline_type == "none":
        baseline = np.zeros_like(t)
    elif baseline_type in {"sine", "spline"}:
        baseline = baseline_model(
            t, baseline_type, parts["baseline_params"], spline_knots=spline_knots
        )
    else:
        raise ValueError(f"Unknown baseline_type: {baseline_type}")
    flare = aflare1(t, parts["tpeak"], parts["fwhm"], parts["ampl"])
    return 1.0 + flare + baseline


def ratio_model_value(
    theta,
    ratio_time,
    ratio_grid,
    baseline_type="none",
    n_spline_knots=5,
    flare_quiescent_flux=None,
):
    parts = split_theta(theta, baseline_type, n_spline_knots, use_ratio_model=True)
    expected_ratio = interpolate_ratio(parts["teff"], parts["ebv"], ratio_grid)
    if not np.isfinite(expected_ratio):
        return np.nan, np.nan
    flare_value = float(
        aflare1(np.array([ratio_time]), parts["tpeak"], parts["fwhm"], parts["ampl"])[0]
    )
    if flare_quiescent_flux is not None:
        flare_value *= float(flare_quiescent_flux)
    return expected_ratio * flare_value, expected_ratio


def log_likelihood(
    theta,
    t,
    flux,
    yerr,
    baseline_type="none",
    spline_knots=None,
    n_spline_knots=5,
    ratio_data=None,
    ratio_grid=None,
    flare_quiescent_flux=None,
):
    sigma = theta[3]
    if not np.isfinite(sigma) or sigma <= 0.0:
        return -np.inf
    model = model_flux(
        theta,
        t,
        baseline_type=baseline_type,
        spline_knots=spline_knots,
        n_spline_knots=n_spline_knots,
        use_ratio_model=ratio_data is not None and ratio_grid is not None,
    )
    variance = yerr**2 + sigma**2
    result = -0.5 * np.sum((flux - model) ** 2 / variance + np.log(2.0 * np.pi * variance))
    if ratio_data is not None:
        if ratio_grid is None:
            raise ValueError("ratio_grid must be provided when ratio_data is used")
        auxiliary_model, _ = ratio_model_value(
            theta,
            ratio_data["time"],
            ratio_grid,
            baseline_type=baseline_type,
            n_spline_knots=n_spline_knots,
            flare_quiescent_flux=(
                flare_quiescent_flux
                if ratio_data.get("input_mode") == "absolute"
                else None
            ),
        )
        if not np.isfinite(auxiliary_model):
            return -np.inf
        auxiliary_variance = ratio_data["err"] ** 2
        result -= 0.5 * (
            (ratio_data["value"] - auxiliary_model) ** 2 / auxiliary_variance
            + np.log(2.0 * np.pi * auxiliary_variance)
        )
    return result