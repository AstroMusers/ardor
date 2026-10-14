"""Very basic dynesty nested-sampling example for a single aflare1 flare model.

Model:
    flux(t) = 1 + aflare1(t, tpeak, fwhm, ampl) + baseline(t)

Optional auxiliary constraint:
    aux_value(t_aux) = ratio(Teff, E(B-V)) * flare_excess(t_aux)

Baseline options:
    - none
    - sine
    - spline (cubic spline through free knot values)

This script reads flare/ratio CSV data, defines a Gaussian likelihood, and
runs dynesty to infer flare, baseline, and intrinsic scatter parameters.
"""

import argparse
from functools import partial
import multiprocessing as mp
from pathlib import Path
from ardor.Flares.aflare import aflare1
from ardor.Roman.flare_model import (
    RATIO_GRID_FILES,
    baseline_model,
    ensure_ratio_param_specs,
    infer_ratio_filter_from_path,
    infer_ratio_grid_bounds,
    interpolate_ratio,
    load_flare_data_csv,
    load_ratio_data_csv,
    load_ratio_grid,
    log_likelihood,
    model_flux,
    normalize_ratio_filter,
    parameter_names,
    ratio_model_value,
    resolve_ratio_grid_path,
    split_theta,
    truth_vector_from_names,
)
import numpy as np
import pandas as pd
from scipy.special import erfcinv
from tdpy.numerics import (
    build_kde_prior_from_samples as _build_kde_prior,
    prepare_kde_prior_spec as _prepare_kde_prior_spec,
    prepare_parameter_prior_specs as _prepare_param_specs,
    weighted_median as _weighted_median,
    weighted_quantile as _weighted_quantile,
)


try:
    from dynesty import NestedSampler
except ImportError as exc:
    raise ImportError(
        "dynesty is required for this script. Install with: pip install dynesty"
    ) from exc

def build_kde_prior_from_samples(
    samples,
    low=None,
    high=None,
    bandwidth=None,
    grid_size=4096,
):
    """Build an inverse-CDF representation of a 1D KDE prior."""
    if grid_size < 256:
        raise ValueError("KDE prior grid_size must be at least 256.")
    prior = _build_kde_prior(samples, low, high, bandwidth, grid_size, normalization="rectangle")
    prior["x_grid"] = prior.pop("cdf_x_grid")
    prior.pop("pdf_grid")
    return prior


def prepare_custom_prior_spec(spec, custom_prior_samples=None):
    """Resolve and precompute custom prior metadata for fast prior transforms."""
    return _prepare_kde_prior_spec(
        spec, custom_prior_samples, prior_builder=build_kde_prior_from_samples,
        pipeline_name="nested_sampling_pipeline()",
    )


def prepare_param_specs(param_specs, custom_prior_samples=None):
    """Precompute any custom prior metadata needed before sampling."""
    return _prepare_param_specs(
        param_specs, custom_prior_samples, prepare_custom_prior=prepare_custom_prior_spec,
    )


def apply_prior_transform_value(
    u_i,
    prior_type,
    low=None,
    high=None,
    mean=None,
    std=None,
    fixed_value=None,
    cdf_grid=None,
    x_grid=None,
):
    """Apply a single prior transformation based on type.

    Parameters:
    - u_i : float in [0, 1] from unit cube
    - prior_type : str ('none', 'uniform', 'gaussian', 'fixed', or 'custom')
    - For 'uniform': low, high define the range
    - For 'gaussian': mean, std define the distribution (unbounded)
    """
    if prior_type == "none":
        return u_i

    if prior_type == "uniform":
        if low is None or high is None:
            raise ValueError(f"Uniform prior requires low and high bounds.")
        return low + (high - low) * u_i

    if prior_type == "gaussian":
        if mean is None or std is None:
            raise ValueError(f"Gaussian prior requires mean and std.")
        # Map [0,1] -> [-inf, inf] via inverse normal CDF
        z = np.sqrt(2.0) * erfcinv(2.0 * u_i)
        return mean + std * z

    if prior_type == "fixed":
        if fixed_value is None:
            raise ValueError("Fixed prior requires fixed_value.")
        return float(fixed_value)

    if prior_type == "custom":
        if cdf_grid is None or x_grid is None:
            raise ValueError("Custom prior requires precomputed cdf_grid and x_grid.")
        u_eval = float(np.clip(u_i, 0.0, 1.0))
        return float(np.interp(u_eval, cdf_grid, x_grid))

    raise ValueError(f"Unknown prior_type: {prior_type}")


def prior_transform_from_csv(u, param_specs):
    """Build prior transformation from a list of parameter specifications.

    Parameters:
    - u : array of shape (ndim,), uniform random in [0, 1]
    - param_specs : list of dicts, each with keys such as:
        {'name', 'prior_type', 'low', 'high', 'mean', 'std', 'fixed_value'}
    """
    params = []
    for i, spec in enumerate(param_specs):
        name = spec["name"]
        prior_type = spec.get("prior_type", "uniform").lower()
        low = spec.get("low")
        high = spec.get("high")
        mean = spec.get("mean")
        std = spec.get("std")
        fixed_value = spec.get("fixed_value")
        cdf_grid = spec.get("cdf_grid")
        x_grid = spec.get("x_grid")

        val = apply_prior_transform_value(
            u[i],
            prior_type,
            low=low,
            high=high,
            mean=mean,
            std=std,
            fixed_value=fixed_value,
            cdf_grid=cdf_grid,
            x_grid=x_grid,
        )
        params.append(val)

    return np.array(params)


def prior_transform(
    u,
    tmin,
    tmax,
    baseline_type="none",
    n_spline_knots=5,
    param_specs=None,
    use_ratio_model=False,
    ratio_grid=None,
):
    """Map unit-cube parameters u in [0,1]^N to physical parameters.

    If param_specs is provided, use CSV-based priors. Otherwise, use defaults.
    """
    if param_specs is not None:
        return prior_transform_from_csv(u, param_specs)

    # Default hard-coded priors (backward compatibility)
    tpeak = tmin + (tmax - tmin) * u[0]
    fwhm = 0.002 + 0.078 * u[1]
    ampl = 0.001 + 0.399 * u[2]
    sigma = 1.0e-5 + 0.01999 * u[3]

    params = [tpeak, fwhm, ampl, sigma]

    if baseline_type == "sine":
        base_c0 = -0.02 + 0.04 * u[4]
        base_amp = 0.0 + 0.03 * u[5]
        base_freq = 0.5 + 7.5 * u[6]
        base_phase = -np.pi + 2.0 * np.pi * u[7]
        params.extend([base_c0, base_amp, base_freq, base_phase])
    elif baseline_type == "spline":
        for i in range(n_spline_knots):
            params.append(-0.02 + 0.04 * u[4 + i])

    if use_ratio_model:
        if ratio_grid is None:
            raise ValueError("ratio_grid is required when use_ratio_model=True")
        bounds = infer_ratio_grid_bounds(ratio_grid)
        start = len(params)
        teff = bounds["teff_min"] + (bounds["teff_max"] - bounds["teff_min"]) * u[start]
        ebv = bounds["ebv_min"] + (bounds["ebv_max"] - bounds["ebv_min"]) * u[start + 1]
        params.extend([teff, ebv])

    return np.array(params)


def weighted_quantile(values, quantiles, sample_weight):
    """Compute weighted quantiles for posterior summaries."""
    return _weighted_quantile(values, quantiles, sample_weight=sample_weight)


def summarize_posterior(samples, weights, names):
    """Print median and 68% credible interval for each parameter."""
    q16, q50, q84 = 0.16, 0.50, 0.84

    print("\nPosterior summary (median +/- 1 sigma):")
    for i, name in enumerate(names):
        vals = samples[:, i]
        lo, med, hi = weighted_quantile(vals, [q16, q50, q84], weights)
        print(f"  {name:>5s} = {med: .4f} (+{hi - med: .4f}, -{med - lo: .4f})")


def weighted_median_params(samples, weights):
    """Return weighted median parameter vector."""
    return _weighted_median(samples, sample_weight=weights, axis=0)


def make_corner_plot(samples, weights, names, truths=None, output_path="ns_corner.png"):
    """Create and save a corner plot of posterior samples.

    Requires the `corner` package and matplotlib.
    """
    try:
        import matplotlib.pyplot as plt
        import corner
    except ImportError:
        print(
            "\nCorner plot not created: install plotting deps with "
            "`pip install corner matplotlib`."
        )
        return

    fig = corner.corner(
        samples,
        labels=names,
        truths=truths,
        quantiles=[0.16, 0.50, 0.84],
        show_titles=True,
        title_fmt=".3f",
        weights=weights,
        smooth=1.0,
    )
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Corner plot saved to: {output_path}")


def make_data_model_plot(
    t,
    flux,
    yerr,
    samples,
    weights,
    baseline_type="none",
    n_spline_knots=5,
    spline_knots=None,
    output_path="ns_data_model.png",
    n_draws=80,
    seed=11,
    ratio_data=None,
    ratio_grid=None,
    flare_quiescent_flux=None,
):
    """Plot data, weighted-median model, and posterior sample models."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nData/model plot not created: install matplotlib with `pip install matplotlib`.")
        return

    tgrid = np.linspace(np.min(t), np.max(t), 600)

    rng = np.random.default_rng(seed)
    n = min(n_draws, samples.shape[0])
    draw_idx = rng.choice(samples.shape[0], size=n, replace=True, p=weights / np.sum(weights))

    if ratio_data is not None and ratio_grid is not None:
        fig, axes = plt.subplots(
            1,
            3,
            figsize=(12.6, 5.0),
            gridspec_kw={"width_ratios": [3.3, 1.7, 1.4]},
        )
        ax = axes[0]
        aux_ax = axes[1]
        ratio_ax = axes[2]
    else:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        aux_ax = None
        ratio_ax = None

    for idx in draw_idx:
        ax.plot(
            tgrid,
            model_flux(
                samples[idx],
                tgrid,
                baseline_type=baseline_type,
                spline_knots=spline_knots,
                n_spline_knots=n_spline_knots,
                use_ratio_model=ratio_data is not None,
            ),
            color="tab:blue",
            alpha=0.10,
            lw=1.0,
        )

    theta_best = weighted_median_params(samples, weights)
    ax.plot(
        tgrid,
        model_flux(
            theta_best,
            tgrid,
            baseline_type=baseline_type,
            spline_knots=spline_knots,
            n_spline_knots=n_spline_knots,
            use_ratio_model=ratio_data is not None,
        ),
        color="tab:red",
        lw=2.2,
        label="Best-fit (weighted median)",
    )

    tpeak_best, fwhm_best, ampl_best = theta_best[:3]
    flare_only = 1.0 + aflare1(tgrid, tpeak_best, fwhm_best, ampl_best)
    ax.plot(tgrid, flare_only, color="tab:orange", lw=1.8, ls="--", label="Flare only")

    ax.errorbar(
        t,
        flux,
        yerr=yerr,
        fmt="o",
        ms=4,
        alpha=0.85,
        color="black",
        ecolor="gray",
        elinewidth=1.0,
        capsize=2,
        label="Observed data",
    )

    if ratio_data is not None:
        ax.axvline(ratio_data["time"], color="tab:green", lw=1.5, ls=":", label="Auxiliary point time")

        if aux_ax is not None:
            aux_predictions = []
            ratio_predictions = []
            for idx in draw_idx:
                aux_value, ratio_value = ratio_model_value(
                    samples[idx],
                    ratio_data["time"],
                    ratio_grid,
                    baseline_type=baseline_type,
                    n_spline_knots=n_spline_knots,
                    flare_quiescent_flux=flare_quiescent_flux if ratio_data.get("input_mode") == "absolute" else None,
                )
                aux_predictions.append(aux_value)
                ratio_predictions.append(ratio_value)

            aux_predictions = np.asarray(aux_predictions, dtype=float)
            aux_predictions = aux_predictions[np.isfinite(aux_predictions)]
            ratio_predictions = np.asarray(ratio_predictions, dtype=float)
            ratio_predictions = ratio_predictions[np.isfinite(ratio_predictions)]

            best_aux_value, best_ratio_value = ratio_model_value(
                theta_best,
                ratio_data["time"],
                ratio_grid,
                baseline_type=baseline_type,
                n_spline_knots=n_spline_knots,
                flare_quiescent_flux=flare_quiescent_flux if ratio_data.get("input_mode") == "absolute" else None,
            )

            jitter = rng.uniform(-0.12, 0.12, size=aux_predictions.size) if aux_predictions.size > 0 else np.array([])
            if aux_predictions.size > 0:
                aux_ax.scatter(
                    aux_predictions,
                    jitter,
                    s=18,
                    color="tab:blue",
                    alpha=0.25,
                    edgecolors="none",
                    label="Posterior samples",
                )

                q16, q50, q84 = weighted_quantile(aux_predictions, [0.16, 0.50, 0.84], np.ones_like(aux_predictions))
                aux_ax.axvspan(q16, q84, color="tab:blue", alpha=0.10)
                aux_ax.axvline(q50, color="tab:blue", lw=1.5, ls="--", label="Sample median")

            aux_ax.axvspan(
                ratio_data["value"] - ratio_data["err"],
                ratio_data["value"] + ratio_data["err"],
                color="tab:red",
                alpha=0.15,
            )
            aux_ax.axvline(ratio_data["value"], color="tab:red", lw=2.0, label="Observed point")
            if np.isfinite(best_aux_value):
                aux_ax.axvline(best_aux_value, color="black", lw=1.8, ls="-.", label="Best-fit prediction")

            aux_ax.set_ylim(-0.2, 0.2)
            aux_ax.set_yticks([])
            if ratio_data.get("input_mode") == "absolute":
                aux_ax.set_xlabel("Auxiliary Flux (absolute)")
            else:
                aux_ax.set_xlabel("Auxiliary Value")
            aux_ax.set_title("Auxiliary Posterior")
            aux_ax.grid(alpha=0.25, axis="x")
            aux_ax.legend(loc="best", fontsize=8)

            if ratio_ax is not None and ratio_predictions.size > 0:
                ratio_jitter = rng.uniform(-0.12, 0.12, size=ratio_predictions.size)
                ratio_ax.scatter(
                    ratio_predictions,
                    ratio_jitter,
                    s=18,
                    color="tab:purple",
                    alpha=0.28,
                    edgecolors="none",
                    label="Posterior samples",
                )
                rq16, rq50, rq84 = weighted_quantile(
                    ratio_predictions,
                    [0.16, 0.50, 0.84],
                    np.ones_like(ratio_predictions),
                )
                ratio_ax.axvspan(rq16, rq84, color="tab:purple", alpha=0.10)
                ratio_ax.axvline(rq50, color="tab:purple", lw=1.5, ls="--", label="Sample median")
                if np.isfinite(best_ratio_value):
                    ratio_ax.axvline(best_ratio_value, color="black", lw=1.8, ls="-.", label="Best-fit")

                flare_model_at_aux = float(
                    aflare1(
                        np.array([ratio_data["time"]]),
                        tpeak_best,
                        fwhm_best,
                        ampl_best,
                    )[0]
                )
                if ratio_data.get("input_mode") == "absolute" and flare_quiescent_flux is not None:
                    flare_model_at_aux *= float(flare_quiescent_flux)

                if np.isfinite(flare_model_at_aux) and abs(flare_model_at_aux) > 0.0:
                    ratio_obs = ratio_data["value"] / flare_model_at_aux
                    ratio_obs_err = ratio_data["err"] / abs(flare_model_at_aux)
                    if np.isfinite(ratio_obs) and np.isfinite(ratio_obs_err):
                        ratio_ax.axvspan(
                            ratio_obs - ratio_obs_err,
                            ratio_obs + ratio_obs_err,
                            color="tab:red",
                            alpha=0.12,
                        )
                        ratio_ax.axvline(ratio_obs, color="tab:red", lw=2.0, label="Observed-implied")

                ratio_ax.set_ylim(-0.2, 0.2)
                ratio_ax.set_yticks([])
                ratio_ax.set_xlabel("Ratio")
                ratio_ax.set_title("Ratio Posterior")
                ratio_ax.grid(alpha=0.25, axis="x")
                ratio_ax.legend(loc="best", fontsize=8)

    ax.set_xlabel("Time")
    ax.set_ylabel("Relative Flux")
    ax.legend(loc="best")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    print(f"Data/model plot saved to: {output_path}")


def run_nested_sampling(
    t,
    flux,
    yerr,
    baseline_type="none",
    n_spline_knots=5,
    spline_knots=None,
    nlive=200,
    dlogz=0.5,
    maxiter=None,
    n_cores=1,
    param_specs=None,
    ratio_data=None,
    ratio_grid=None,
    flare_quiescent_flux=None,
):
    """Run dynesty NestedSampler on the aflare1 model with baseline."""

    tmin = float(np.min(t))
    tmax = float(np.max(t))

    if baseline_type == "spline" and spline_knots is None:
        spline_knots = np.linspace(tmin, tmax, n_spline_knots)

    loglike = partial(
        log_likelihood,
        t=t,
        flux=flux,
        yerr=yerr,
        baseline_type=baseline_type,
        spline_knots=spline_knots,
        n_spline_knots=n_spline_knots,
        ratio_data=ratio_data,
        ratio_grid=ratio_grid,
        flare_quiescent_flux=flare_quiescent_flux,
    )

    if param_specs is not None:
        ptform = partial(prior_transform, tmin=tmin, tmax=tmax, param_specs=param_specs)
    else:
        ptform = partial(
            prior_transform,
            tmin=tmin,
            tmax=tmax,
            baseline_type=baseline_type,
            n_spline_knots=n_spline_knots,
            use_ratio_model=ratio_data is not None and ratio_grid is not None,
            ratio_grid=ratio_grid,
        )

    if param_specs is not None:
        ndim = len(param_specs)
    else:
        ndim = len(
            parameter_names(
                baseline_type=baseline_type,
                n_spline_knots=n_spline_knots,
                use_ratio_model=ratio_data is not None and ratio_grid is not None,
            )
        )

    n_cores = max(1, int(n_cores))
    if n_cores == 1:
        sampler = NestedSampler(loglike, ptform, ndim=ndim, nlive=nlive)
        sampler.run_nested(dlogz=dlogz, maxiter=maxiter, print_progress=True)
    else:
        # Use multiprocessing to evaluate likelihoods in parallel.
        with mp.Pool(processes=n_cores) as pool:
            sampler = NestedSampler(
                loglike,
                ptform,
                ndim=ndim,
                nlive=nlive,
                pool=pool,
                queue_size=n_cores,
            )
            sampler.run_nested(dlogz=dlogz, maxiter=maxiter, print_progress=True)

    return sampler.results


def load_settings_csv(filepath):
    """Load settings from a CSV file.

    Expected columns: baseline_model, n_spline_knots, nlive, dlogz, maxiter, n_cores,
    flare_data_csv, use_ratio_model, ratio_filter, ratio_grid_csv, ratio_data_csv,
    aux_quiescent_flux, flare_quiescent_flux
    Returns a dictionary with these settings.
    """
    import pandas as pd
    
    df = pd.read_csv(filepath)
    if len(df) != 1:
        raise ValueError(f"settings.csv must have exactly 1 row, got {len(df)}")
    
    row = df.iloc[0]

    def maybe_bool(value, default=False):
        if pd.isna(value):
            return default
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "y"}
        return bool(value)

    settings = {
        "baseline_model": row.get("baseline_model", "none"),
        "n_spline_knots": int(row.get("n_spline_knots", 5)),
        "nlive": int(row.get("nlive", 200)),
        "dlogz": float(row.get("dlogz", 0.5)),
        "maxiter": int(row.get("maxiter", None)) if pd.notna(row.get("maxiter")) else None,
        "n_cores": max(1, int(row.get("n_cores", 1))) if pd.notna(row.get("n_cores", 1)) else 1,
        "flare_data_csv": row.get("flare_data_csv", None) if pd.notna(row.get("flare_data_csv", None)) else None,
        "use_ratio_model": maybe_bool(row.get("use_ratio_model", False), default=False),
        "ratio_filter": normalize_ratio_filter(row.get("ratio_filter", "F087")),
        "ratio_grid_csv": row.get("ratio_grid_csv", None) if pd.notna(row.get("ratio_grid_csv", None)) else None,
        "ratio_data_csv": row.get("ratio_data_csv", None) if pd.notna(row.get("ratio_data_csv", None)) else None,
        "aux_quiescent_flux": float(row.get("aux_quiescent_flux")) if pd.notna(row.get("aux_quiescent_flux")) else None,
        "flare_quiescent_flux": float(row.get("flare_quiescent_flux")) if pd.notna(row.get("flare_quiescent_flux")) else None,
    }
    return settings


def load_params_csv(filepath):
    """Load parameter specifications from a CSV file.
    
    Expected columns: name, prior_type, low, high, mean, std, fixed_value(optional),
    bandwidth(optional), grid_size(optional), custom_kind(optional)
    prior_type must be one of: 'none', 'uniform', 'gaussian', 'fixed', 'custom'
    
    Returns a list of parameter specification dicts.
    """
    import pandas as pd
    
    df = pd.read_csv(filepath)
    param_specs = []
    
    for idx, row in df.iterrows():
        spec = {
            'name': row['name'],
            'prior_type': row['prior_type'].strip().lower(),
        }
        
        if spec['prior_type'] not in ['none', 'uniform', 'gaussian', 'fixed', 'custom']:
            raise ValueError(f"Unknown prior_type: {spec['prior_type']} at row {idx}")
        
        if spec['prior_type'] == 'uniform':
            spec['low'] = float(row['low'])
            spec['high'] = float(row['high'])
        elif spec['prior_type'] == 'gaussian':
            spec['mean'] = float(row['mean'])
            spec['std'] = float(row['std'])
        elif spec['prior_type'] == 'fixed':
            fixed_value = None
            if 'fixed_value' in row.index and pd.notna(row['fixed_value']):
                fixed_value = float(row['fixed_value'])
            elif 'mean' in row.index and pd.notna(row['mean']):
                fixed_value = float(row['mean'])
            elif 'low' in row.index and pd.notna(row['low']):
                fixed_value = float(row['low'])
            elif 'high' in row.index and pd.notna(row['high']):
                fixed_value = float(row['high'])
            if fixed_value is None:
                raise ValueError(
                    f"Fixed prior requires fixed_value (or mean/low/high fallback) at row {idx}"
                )
            spec['fixed_value'] = fixed_value
        elif spec['prior_type'] == 'custom':
            custom_kind = 'kde'
            if 'custom_kind' in row.index and pd.notna(row['custom_kind']):
                custom_kind = str(row['custom_kind']).strip().lower()
            elif 'prior_kind' in row.index and pd.notna(row['prior_kind']):
                custom_kind = str(row['prior_kind']).strip().lower()
            spec['custom_kind'] = custom_kind
            if 'bandwidth' in row.index and pd.notna(row['bandwidth']):
                spec['bandwidth'] = float(row['bandwidth'])
            if 'grid_size' in row.index and pd.notna(row['grid_size']):
                spec['grid_size'] = int(row['grid_size'])
            if 'low' in row.index and pd.notna(row['low']):
                spec['low'] = float(row['low'])
            if 'high' in row.index and pd.notna(row['high']):
                spec['high'] = float(row['high'])
        
        param_specs.append(spec)
    
    return param_specs


def save_settings_csv(filepath, settings_dict):
    """Save settings dictionary to a CSV file.
    
    Creates a single-row CSV with all settings columns.
    """
    import pandas as pd
    
    df = pd.DataFrame([settings_dict])
    df.to_csv(filepath, index=False)
    print(f"Settings CSV saved to: {filepath}")


def save_params_csv(filepath, param_specs):
    """Save parameter specifications to a CSV file.
    
    Creates a CSV with columns: name, prior_type, low, high, mean, std,
    fixed_value, bandwidth, grid_size, custom_kind
    """
    import pandas as pd
    
    rows = []
    for spec in param_specs:
        row = {
            'name': spec['name'],
            'prior_type': spec['prior_type'],
            'low': spec.get('low', ''),
            'high': spec.get('high', ''),
            'mean': spec.get('mean', ''),
            'std': spec.get('std', ''),
            'fixed_value': spec.get('fixed_value', ''),
            'bandwidth': spec.get('bandwidth', ''),
            'grid_size': spec.get('grid_size', ''),
            'custom_kind': spec.get('custom_kind', ''),
        }
        rows.append(row)
    
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)
    print(f"Parameters CSV saved to: {filepath}")


def nested_sampling_pipeline(
    flare_data_csv,
    params_csv=None,
    settings_csv=None,
    custom_prior_samples=None,
    ratio_data_csv=None,
    ratio_grid_csv=None,
    baseline_model="none",
    n_spline_knots=5,
    nlive=200,
    dlogz=0.5,
    maxiter=None,
    n_cores=1,
    use_ratio_model=False,
    ratio_filter="F087",
    aux_quiescent_flux=None,
    flare_quiescent_flux=None,
    output_dir=".",
    plot_prefix="ns",
    update_settings_csv=False,
    update_params_csv=False,
):
    """
    Main pipeline for nested-sampling flare and auxiliary data analysis.
    
    Loads CSV configuration and data, runs nested sampling, generates plots and summaries.
    
    Parameters:
    - flare_data_csv: Path to flare light-curve CSV (time, flux, yerr)
    - params_csv: Path to parameter priors CSV (name, prior_type, low, high, mean, std)
    - settings_csv: Path to settings CSV (overrides other kwargs if provided)
        - custom_prior_samples: Dict mapping parameter names to sample arrays or dicts
            with keys like samples, low, high, bandwidth, grid_size, custom_kind.
    - ratio_data_csv: Path to single auxiliary point CSV (time, value, err) or (time, flux_abs, err_abs)
    - ratio_grid_csv: Path to Teff/E(B-V) ratio grid CSV
    - baseline_model: "none", "sine", or "spline"
    - n_spline_knots: Number of spline knots if baseline_model="spline"
    - nlive: Number of live points for dynesty
    - dlogz: Stopping criterion on remaining log-evidence
    - maxiter: Optional hard cap on nested-sampling iterations
    - n_cores: Number of CPU cores for parallel likelihood evaluation
    - use_ratio_model: Enable Teff/E(B-V) ratio model tied to auxiliary point
    - ratio_filter: Named ratio grid ("F087" or "F213")
    - aux_quiescent_flux: Auxiliary-band quiescent absolute flux (required for flux_abs/err_abs mode)
    - flare_quiescent_flux: Flare-band quiescent absolute flux (required for flux_abs/err_abs mode)
    - output_dir: Directory for output plots
    - plot_prefix: Prefix for output plot filenames
    - update_settings_csv: If True, write the final settings dict back to settings_csv file
    - update_params_csv: If True, write the final param_specs back to params_csv file
    
    Returns: (results, samples, weights, names)
    """
    # Load settings from CSV if provided (overrides kwargs)
    if settings_csv is not None:
        settings = load_settings_csv(settings_csv)
        baseline_model = settings.get("baseline_model", baseline_model)
        n_spline_knots = settings.get("n_spline_knots", n_spline_knots)
        nlive = settings.get("nlive", nlive)
        dlogz = settings.get("dlogz", dlogz)
        maxiter = settings.get("maxiter", maxiter)
        n_cores = settings.get("n_cores", n_cores)
        flare_data_csv = settings.get("flare_data_csv", flare_data_csv)
        use_ratio_model = settings.get("use_ratio_model", use_ratio_model)
        ratio_filter = settings.get("ratio_filter", ratio_filter)
        # Preserve explicitly passed kwargs: only fall back to settings.csv values
        # when the caller did not provide a path or quiescent flux.
        if ratio_grid_csv is None and settings.get("ratio_grid_csv") is not None:
            ratio_grid_csv = settings["ratio_grid_csv"]
        if ratio_data_csv is None and settings.get("ratio_data_csv") is not None:
            ratio_data_csv = settings["ratio_data_csv"]
        if aux_quiescent_flux is None and settings.get("aux_quiescent_flux") is not None:
            aux_quiescent_flux = settings["aux_quiescent_flux"]
        if flare_quiescent_flux is None and settings.get("flare_quiescent_flux") is not None:
            flare_quiescent_flux = settings["flare_quiescent_flux"]
    
    # Load parameter specifications from CSV if provided
    param_specs = None
    if params_csv is not None:
        param_specs = load_params_csv(params_csv)
        param_specs = prepare_param_specs(param_specs, custom_prior_samples=custom_prior_samples)
    
    # Load and prepare ratio grid if using ratio model
    ratio_grid = None
    if use_ratio_model:
        ratio_grid_path = resolve_ratio_grid_path(ratio_filter=ratio_filter, ratio_grid_csv=ratio_grid_csv)
        ratio_grid = load_ratio_grid(ratio_grid_path)
        if param_specs is not None:
            param_specs = ensure_ratio_param_specs(param_specs, ratio_grid)
    
    # Load flare light-curve data
    if flare_data_csv is None:
        raise ValueError("flare_data_csv is required")
    t, flux, yerr = load_flare_data_csv(flare_data_csv)
    
    # Load auxiliary ratio data if using ratio model
    ratio_data = None
    if use_ratio_model and ratio_data_csv is not None:
        ratio_data = load_ratio_data_csv(
            ratio_data_csv,
            aux_quiescent_flux=aux_quiescent_flux,
            flare_quiescent_flux=flare_quiescent_flux,
        )
    
    # Build metadata for plot
    truth = {"baseline_type": baseline_model}
    spline_knots = None
    if use_ratio_model and ratio_data is not None:
        truth["ratio_filter"] = infer_ratio_filter_from_path(ratio_grid["path"])
        if aux_quiescent_flux is not None:
            truth["aux_quiescent_flux"] = float(aux_quiescent_flux)
        if flare_quiescent_flux is not None:
            truth["flare_quiescent_flux"] = float(flare_quiescent_flux)
    
    print("Model/input metadata:")
    print(truth)
    if ratio_data is not None:
        print("Auxiliary single-point data:")
        print(ratio_data)
    
    print(
        f"Running nested sampling with baseline={baseline_model}, "
        f"nlive={nlive}, dlogz={dlogz}, maxiter={maxiter}, "
        f"n_cores={n_cores}, use_ratio_model={use_ratio_model}"
    )
    
    # Optionally save settings to CSV before running
    if update_settings_csv and settings_csv is not None:
        settings_dict = {
            "baseline_model": baseline_model,
            "n_spline_knots": n_spline_knots,
            "nlive": nlive,
            "dlogz": dlogz,
            "maxiter": maxiter if maxiter is not None else "",
            "n_cores": n_cores,
            "flare_data_csv": flare_data_csv,
            "use_ratio_model": use_ratio_model,
            "ratio_filter": ratio_filter,
            "ratio_grid_csv": ratio_grid_csv if ratio_grid_csv is not None else "",
            "ratio_data_csv": ratio_data_csv if ratio_data_csv is not None else "",
            "aux_quiescent_flux": aux_quiescent_flux if aux_quiescent_flux is not None else "",
            "flare_quiescent_flux": flare_quiescent_flux if flare_quiescent_flux is not None else "",
        }
        save_settings_csv(settings_csv, settings_dict)
    
    # Optionally save params to CSV before running
    if update_params_csv and params_csv is not None and param_specs is not None:
        save_params_csv(params_csv, param_specs)
    
    # Run nested sampling
    results = run_nested_sampling(
        t,
        flux,
        yerr,
        baseline_type=baseline_model,
        n_spline_knots=n_spline_knots,
        spline_knots=spline_knots,
        nlive=nlive,
        dlogz=dlogz,
        maxiter=maxiter,
        n_cores=n_cores,
        param_specs=param_specs,
        ratio_data=ratio_data,
        ratio_grid=ratio_grid,
        flare_quiescent_flux=flare_quiescent_flux,
    )
    
    # Convert log-weights to normalized linear weights
    weights = np.exp(results.logwt - results.logz[-1])
    weights /= np.sum(weights)
    samples = results.samples
    
    # Build parameter names
    if param_specs is not None:
        names = [spec["name"] for spec in param_specs]
    else:
        names = parameter_names(baseline_model, n_spline_knots, use_ratio_model=use_ratio_model)
    
    # For summary/corner plotting, drop parameters explicitly fixed by prior.
    if param_specs is not None:
        plot_indices = [
            i for i, spec in enumerate(param_specs)
            if spec.get("prior_type", "uniform").lower() != "fixed"
        ]
        plot_names = [names[i] for i in plot_indices]
        plot_samples = samples[:, plot_indices] if len(plot_indices) > 0 else np.empty((samples.shape[0], 0))
    else:
        plot_names = names
        plot_samples = samples

    if len(plot_names) > 0:
        truths = truth_vector_from_names(plot_names, truth)
        summarize_posterior(plot_samples, weights, names=plot_names)
        make_corner_plot(
            plot_samples,
            weights,
            names=plot_names,
            truths=truths,
            output_path=f"{output_dir}/{plot_prefix}_corner_{baseline_model}.png",
        )
    else:
        print("All parameters are fixed; skipping posterior summary and corner plot.")

    # Data/model plotting still uses the full parameter vector.
    make_data_model_plot(
        t,
        flux,
        yerr,
        samples,
        weights,
        baseline_type=baseline_model,
        n_spline_knots=n_spline_knots,
        spline_knots=spline_knots,
        output_path=f"{output_dir}/{plot_prefix}_data_model_{baseline_model}.png",
        n_draws=80,
        ratio_data=ratio_data,
        ratio_grid=ratio_grid,
        flare_quiescent_flux=flare_quiescent_flux,
    )
    
    return results, samples, weights, names


def parse_args():
    """Parse command-line options for baseline model configuration."""
    parser = argparse.ArgumentParser(description="dynesty example for aflare1 flare fitting")
    parser.add_argument(
        "--baseline-model",
        choices=["none", "sine", "spline"],
        default="none",
        help="Baseline model to include in the flare fit.",
    )
    parser.add_argument(
        "--n-spline-knots",
        type=int,
        default=5,
        help="Number of spline knots if --baseline-model spline.",
    )
    parser.add_argument(
        "--nlive",
        type=int,
        default=200,
        help="Number of live points for dynesty.",
    )
    parser.add_argument(
        "--dlogz",
        type=float,
        default=0.5,
        help="Stopping criterion on remaining log-evidence.",
    )
    parser.add_argument(
        "--maxiter",
        type=int,
        default=None,
        help="Optional hard cap on nested-sampling iterations.",
    )
    parser.add_argument(
        "--n-cores",
        type=int,
        default=1,
        help="Number of CPU cores for dynesty likelihood evaluation.",
    )
    parser.add_argument(
        "--settings-csv",
        type=str,
        default=None,
        help="Path to CSV file with settings including model choices and optional input-data CSVs.",
    )
    parser.add_argument(
        "--params-csv",
        type=str,
        default=None,
        help="Path to CSV file with parameter priors (name, prior_type, low, high, mean, std).",
    )
    parser.add_argument(
        "--flare-data-csv",
        type=str,
        default=None,
        help="Path to a flare light-curve CSV with columns time,flux,yerr.",
    )
    parser.add_argument(
        "--use-ratio-model",
        action="store_true",
        help="Include a concurrent Teff/E(B-V) ratio model tied to a single auxiliary data point.",
    )
    parser.add_argument(
        "--ratio-filter",
        choices=sorted(RATIO_GRID_FILES),
        default="F087",
        help="Named ratio grid to use for the auxiliary model.",
    )
    parser.add_argument(
        "--ratio-grid-csv",
        type=str,
        default=None,
        help="Path to the Teff/E(B-V) ratio grid CSV. Overrides --ratio-filter if set.",
    )
    parser.add_argument(
        "--ratio-data-csv",
        type=str,
        default=None,
        help=(
            "Path to a one-row auxiliary CSV with either columns time,value,err "
            "or time,flux_abs,err_abs."
        ),
    )
    parser.add_argument(
        "--aux-quiescent-flux",
        type=float,
        default=None,
        help=(
            "Quiescent absolute flux in the auxiliary band. Required when "
            "ratio_data_csv uses columns flux_abs,err_abs."
        ),
    )
    parser.add_argument(
        "--flare-quiescent-flux",
        type=float,
        default=None,
        help=(
            "Quiescent absolute flux in the flare band. Required when "
            "ratio_data_csv uses columns flux_abs,err_abs."
        ),
    )
    return parser.parse_args()


def main():
    """CLI entry point: parse arguments and invoke the nested-sampling pipeline."""
    args = parse_args()
    
    nested_sampling_pipeline(
        flare_data_csv=args.flare_data_csv,
        params_csv=args.params_csv,
        settings_csv=args.settings_csv,
        ratio_data_csv=args.ratio_data_csv,
        ratio_grid_csv=args.ratio_grid_csv,
        baseline_model=args.baseline_model,
        n_spline_knots=args.n_spline_knots,
        nlive=args.nlive,
        dlogz=args.dlogz,
        maxiter=args.maxiter,
        n_cores=args.n_cores,
        use_ratio_model=args.use_ratio_model,
        ratio_filter=args.ratio_filter,
        aux_quiescent_flux=args.aux_quiescent_flux,
        flare_quiescent_flux=args.flare_quiescent_flux,
        update_settings_csv=False,
        update_params_csv=False,
    )


if __name__ == "__main__":
    main()
