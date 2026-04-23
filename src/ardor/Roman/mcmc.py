"""emcee MCMC interface matching the ns.py model and CSV workflow.

Model:
	flux(t) = 1 + aflare1(t, tpeak, fwhm, ampl) + baseline(t)

Optional auxiliary constraint:
	aux_value(t_aux) = ratio(Teff, E(B-V)) * flare_excess(t_aux)

Baseline options:
	- none
	- sine
	- spline (cubic spline through free knot values)
"""

import argparse
import multiprocessing as mp
import sys
from pathlib import Path

from ardor.Roman.flare_color import get_data_file

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, RegularGridInterpolator
from scipy.stats import gaussian_kde


RATIO_GRID_FILES = {
	"F087": "F087.csv",
	"F213": "F213.csv",
}

try:
	import emcee
except ImportError as exc:
	raise ImportError("emcee is required for this script. Install with: pip install emcee") from exc

try:
	from ardor.Flares.Flare import aflare1
except ImportError:
	repo_root = Path(__file__).resolve().parents[1]
	ardor_src = repo_root / "ardor" / "src"
	if ardor_src.exists():
		sys.path.insert(0, str(ardor_src))
	from ardor.Flares.aflare import aflare1


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
		cs = CubicSpline(spline_knots, np.asarray(params), bc_type="natural")
		return cs(t)

	raise ValueError(f"Unknown baseline_type: {baseline_type}")


def load_flare_data_csv(filepath):
	df = pd.read_csv(filepath)
	required = {"time", "flux", "yerr"}
	missing = required.difference(df.columns)
	if missing:
		raise ValueError(f"flare_data.csv is missing required columns: {', '.join(sorted(missing))}")
	if len(df) < 2:
		raise ValueError("flare_data.csv must contain at least 2 rows")

	return (
		df["time"].to_numpy(dtype=float),
		df["flux"].to_numpy(dtype=float),
		df["yerr"].to_numpy(dtype=float),
	)


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
	normalized = normalize_ratio_filter(ratio_filter)
	return Path(get_data_file(RATIO_GRID_FILES[normalized]))


def infer_ratio_filter_from_path(ratio_grid_path):
	filename = Path(ratio_grid_path).name.upper()
	for ratio_filter, mapped_name in RATIO_GRID_FILES.items():
		if filename == mapped_name.upper():
			return ratio_filter
	return "CUSTOM"


def load_ratio_grid(grid_csv_path):
	"""Load the Teff/E(B-V) ratio grid and build a 2D interpolator.

	Expected CSV layout:
	- Teff values in the first column/index.
	- E(B-V) values as raw numeric column headers.
	"""
	df = pd.read_csv(grid_csv_path, index_col=0)
	try:
		teff_values = df.index.to_numpy(dtype=float)
		ebv_values = np.array(df.columns, dtype=float)
	except ValueError as exc:
		raise ValueError(
			"Ratio grid must use a numeric Teff index and numeric E(B-V) column headers: "
			f"{grid_csv_path}"
		) from exc

	ratio_values = df.to_numpy(dtype=float)
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


def build_kde_prior_from_samples(samples, low=None, high=None, bandwidth=None, grid_size=4096):
	samples = np.asarray(samples, dtype=float).reshape(-1)
	samples = samples[np.isfinite(samples)]
	if samples.size < 2:
		raise ValueError("KDE prior requires at least 2 finite samples.")

	if int(grid_size) < 256:
		raise ValueError("KDE prior grid_size must be at least 256.")

	sample_std = float(np.std(samples, ddof=1)) if samples.size > 1 else 0.0
	sample_range = float(np.max(samples) - np.min(samples))
	fallback_width = sample_std if sample_std > 0.0 else max(sample_range, 1.0)
	padding_scale = 3.0 * fallback_width if fallback_width > 0.0 else 1.0

	domain_low = float(low) if low is not None else float(np.min(samples) - padding_scale)
	domain_high = float(high) if high is not None else float(np.max(samples) + padding_scale)
	if not np.isfinite(domain_low) or not np.isfinite(domain_high) or domain_high <= domain_low:
		raise ValueError("KDE prior bounds must define a finite interval with high > low.")

	kde_kwargs = {}
	if bandwidth is not None:
		kde_kwargs["bw_method"] = float(bandwidth)
	kde = gaussian_kde(samples, **kde_kwargs)

	x_grid = np.linspace(domain_low, domain_high, int(grid_size), dtype=float)
	pdf = np.asarray(kde(x_grid), dtype=float)
	pdf = np.clip(pdf, 0.0, None)
	norm_const = float(np.trapz(pdf, x_grid))
	if not np.isfinite(norm_const) or norm_const <= 0.0:
		raise ValueError("KDE prior produced a non-finite or zero integral over the sampling domain.")

	pdf /= norm_const
	dx = float(x_grid[1] - x_grid[0])
	cdf = np.cumsum(pdf) * dx
	cdf[0] = 0.0
	cdf[-1] = 1.0
	unique_cdf, unique_idx = np.unique(cdf, return_index=True)
	if unique_cdf.size < 2:
		raise ValueError("KDE prior CDF is degenerate; broaden bounds or use different samples.")

	return {
		"kind": "kde",
		"x_grid": x_grid,
		"pdf_grid": pdf,
		"cdf_grid": unique_cdf,
		"cdf_x_grid": x_grid[unique_idx],
		"low": domain_low,
		"high": domain_high,
		"bandwidth": None if bandwidth is None else float(bandwidth),
		"grid_size": int(grid_size),
	}


def prepare_custom_prior_spec(spec, custom_prior_samples=None):
	prepared = dict(spec)
	custom_kind = str(prepared.get("custom_kind", prepared.get("prior_kind", "kde"))).strip().lower()
	if custom_kind != "kde":
		raise ValueError(
			f"Unsupported custom prior kind '{custom_kind}' for {prepared.get('name', '?')}."
		)

	if custom_prior_samples is None:
		custom_prior_samples = {}

	param_name = str(prepared.get("name", "")).strip()
	if param_name not in custom_prior_samples:
		raise ValueError(
			f"Custom KDE prior for parameter {prepared.get('name', '?')} requires an entry in "
			"custom_prior_samples passed to mcmc_pipeline()."
		)

	sample_input = custom_prior_samples[param_name]
	if isinstance(sample_input, dict):
		if "samples" not in sample_input:
			raise ValueError(f"custom_prior_samples['{param_name}'] must contain a 'samples' entry.")
		samples = sample_input["samples"]
		bandwidth = sample_input.get("bandwidth", prepared.get("bandwidth"))
		grid_size = sample_input.get("grid_size", prepared.get("grid_size", 4096))
		low = sample_input.get("low", prepared.get("low"))
		high = sample_input.get("high", prepared.get("high"))
		custom_kind = str(sample_input.get("custom_kind", custom_kind)).strip().lower()
	else:
		samples = sample_input
		bandwidth = prepared.get("bandwidth")
		grid_size = prepared.get("grid_size", 4096)
		low = prepared.get("low")
		high = prepared.get("high")

	if custom_kind != "kde":
		raise ValueError(
			f"Unsupported custom prior kind '{custom_kind}' for {prepared.get('name', '?')}."
		)

	prior_data = build_kde_prior_from_samples(
		samples,
		low=low,
		high=high,
		bandwidth=bandwidth,
		grid_size=grid_size,
	)

	prepared["custom_kind"] = custom_kind
	prepared["bandwidth"] = prior_data["bandwidth"]
	prepared["grid_size"] = prior_data["grid_size"]
	prepared["low"] = prior_data["low"]
	prepared["high"] = prior_data["high"]
	prepared["x_grid"] = prior_data["x_grid"]
	prepared["pdf_grid"] = prior_data["pdf_grid"]
	prepared["cdf_grid"] = prior_data["cdf_grid"]
	prepared["cdf_x_grid"] = prior_data["cdf_x_grid"]
	return prepared


def prepare_param_specs(param_specs, custom_prior_samples=None):
	prepared_specs = []
	for spec in param_specs:
		prior_type = str(spec.get("prior_type", "uniform")).strip().lower()
		if prior_type == "custom":
			prepared_specs.append(prepare_custom_prior_spec(spec, custom_prior_samples=custom_prior_samples))
		else:
			prepared_specs.append(dict(spec))
	return prepared_specs


def interpolate_ratio(teff, ebv, ratio_grid):
	point = np.array([float(teff), float(ebv)], dtype=float)
	return float(ratio_grid["interpolator"](point))


def load_ratio_data_csv(filepath, aux_quiescent_flux=None, flare_quiescent_flux=None):
	df = pd.read_csv(filepath)
	if len(df) != 1:
		raise ValueError(f"ratio_data.csv must have exactly 1 row, got {len(df)}")

	row = df.iloc[0]
	if "time" not in row.index:
		raise ValueError("ratio_data.csv is missing required column: time")

	has_relative = all(col in row.index for col in ["value", "err"])
	has_absolute = all(col in row.index for col in ["flux_abs", "err_abs"])

	if has_relative:
		if pd.isna(row["value"]) or pd.isna(row["err"]):
			raise ValueError("ratio_data.csv has missing value/err for relative-input format")
		return {
			"time": float(row["time"]),
			"value": float(row["value"]),
			"err": float(row["err"]),
			"input_mode": "relative",
		}

	if has_absolute:
		if pd.isna(row["flux_abs"]) or pd.isna(row["err_abs"]):
			raise ValueError("ratio_data.csv has missing flux_abs/err_abs for absolute-input format")

		flux_abs = float(row["flux_abs"])
		err_abs = float(row["err_abs"])
		if aux_quiescent_flux is None or flare_quiescent_flux is None:
			raise ValueError(
				"ratio_data.csv uses flux_abs/err_abs. Provide both aux_quiescent_flux "
				"and flare_quiescent_flux as absolute quantities."
			)
		if float(aux_quiescent_flux) <= 0.0:
			raise ValueError("aux_quiescent_flux must be > 0 when using absolute auxiliary flux")
		if float(flare_quiescent_flux) <= 0.0:
			raise ValueError("flare_quiescent_flux must be > 0 when using absolute auxiliary flux")

		aux_quiescent_flux = float(aux_quiescent_flux)
		flare_quiescent_flux = float(flare_quiescent_flux)
		return {
			"time": float(row["time"]),
			"value": flux_abs - aux_quiescent_flux,
			"err": err_abs,
			"input_mode": "absolute",
			"flux_abs": flux_abs,
			"err_abs": err_abs,
			"aux_quiescent_flux": aux_quiescent_flux,
			"flare_quiescent_flux": flare_quiescent_flux,
		}

	raise ValueError(
		"ratio_data.csv must contain either columns (time,value,err) or (time,flux_abs,err_abs)"
	)


def ensure_ratio_param_specs(param_specs, ratio_grid):
	specs = list(param_specs)
	existing_names = {spec["name"].lower() for spec in specs}
	bounds = infer_ratio_grid_bounds(ratio_grid)

	if "teff" not in existing_names:
		specs.append({"name": "teff", "prior_type": "uniform", "low": bounds["teff_min"], "high": bounds["teff_max"]})
	if "ebv" not in existing_names:
		specs.append({"name": "ebv", "prior_type": "uniform", "low": bounds["ebv_min"], "high": bounds["ebv_max"]})
	return specs


def parameter_names(baseline_type="none", n_spline_knots=5, use_ratio_model=False):
	names = ["tpeak", "fwhm", "ampl", "sigma"]
	if baseline_type == "sine":
		names.extend(["base_c0", "base_amp", "base_freq", "base_phase"])
	elif baseline_type == "spline":
		names.extend([f"spline_c{i}" for i in range(n_spline_knots)])
	if use_ratio_model:
		names.extend(["teff", "ebv"])
	return names


def default_param_specs(tmin, tmax, baseline_type="none", n_spline_knots=5, use_ratio_model=False, ratio_grid=None):
	specs = [
		{"name": "tpeak", "prior_type": "uniform", "low": float(tmin), "high": float(tmax)},
		{"name": "fwhm", "prior_type": "uniform", "low": 0.002, "high": 0.08},
		{"name": "ampl", "prior_type": "uniform", "low": 0.001, "high": 0.4},
		{"name": "sigma", "prior_type": "uniform", "low": 1.0e-5, "high": 0.02},
	]

	if baseline_type == "sine":
		specs.extend(
			[
				{"name": "base_c0", "prior_type": "uniform", "low": -0.02, "high": 0.02},
				{"name": "base_amp", "prior_type": "uniform", "low": 0.0, "high": 0.03},
				{"name": "base_freq", "prior_type": "uniform", "low": 0.5, "high": 8.0},
				{"name": "base_phase", "prior_type": "uniform", "low": -np.pi, "high": np.pi},
			]
		)
	elif baseline_type == "spline":
		for i in range(n_spline_knots):
			specs.append({"name": f"spline_c{i}", "prior_type": "uniform", "low": -0.02, "high": 0.02})

	if use_ratio_model:
		if ratio_grid is None:
			raise ValueError("ratio_grid is required when use_ratio_model=True")
		bounds = infer_ratio_grid_bounds(ratio_grid)
		specs.extend(
			[
				{"name": "teff", "prior_type": "uniform", "low": bounds["teff_min"], "high": bounds["teff_max"]},
				{"name": "ebv", "prior_type": "uniform", "low": bounds["ebv_min"], "high": bounds["ebv_max"]},
			]
		)

	return specs


def split_theta(theta, baseline_type="none", n_spline_knots=5, use_ratio_model=False):
	parts = {
		"tpeak": float(theta[0]),
		"fwhm": float(theta[1]),
		"ampl": float(theta[2]),
		"sigma": float(theta[3]),
	}
	idx = 4
	if baseline_type == "sine":
		parts["baseline_params"] = np.asarray(theta[idx : idx + 4], dtype=float)
		idx += 4
	elif baseline_type == "spline":
		parts["baseline_params"] = np.asarray(theta[idx : idx + n_spline_knots], dtype=float)
		idx += n_spline_knots
	else:
		parts["baseline_params"] = np.array([], dtype=float)

	if use_ratio_model:
		parts["teff"] = float(theta[idx])
		parts["ebv"] = float(theta[idx + 1])

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
	elif baseline_type == "sine":
		baseline = baseline_model(t, "sine", parts["baseline_params"], spline_knots=spline_knots)
	elif baseline_type == "spline":
		baseline = baseline_model(t, "spline", parts["baseline_params"], spline_knots=spline_knots)
	else:
		raise ValueError(f"Unknown baseline_type: {baseline_type}")

	return 1.0 + aflare1(t, parts["tpeak"], parts["fwhm"], parts["ampl"]) + baseline


def ratio_model_value(theta, ratio_time, ratio_grid, baseline_type="none", n_spline_knots=5, flare_quiescent_flux=None):
	parts = split_theta(theta, baseline_type=baseline_type, n_spline_knots=n_spline_knots, use_ratio_model=True)
	expected_ratio = interpolate_ratio(parts["teff"], parts["ebv"], ratio_grid)
	if not np.isfinite(expected_ratio):
		return np.nan, np.nan

	flare_value = float(aflare1(np.array([ratio_time]), parts["tpeak"], parts["fwhm"], parts["ampl"])[0])
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
	if sigma <= 0.0:
		return -np.inf

	model = model_flux(
		theta,
		t,
		baseline_type=baseline_type,
		spline_knots=spline_knots,
		n_spline_knots=n_spline_knots,
		use_ratio_model=ratio_data is not None and ratio_grid is not None,
	)
	var = yerr**2 + sigma**2
	loglike = -0.5 * np.sum((flux - model) ** 2 / var + np.log(2.0 * np.pi * var))

	if ratio_data is not None:
		if ratio_grid is None:
			raise ValueError("ratio_grid must be provided when ratio_data is used")
		aux_model, _ = ratio_model_value(
			theta,
			ratio_data["time"],
			ratio_grid,
			baseline_type=baseline_type,
			n_spline_knots=n_spline_knots,
			flare_quiescent_flux=flare_quiescent_flux if ratio_data.get("input_mode") == "absolute" else None,
		)
		if not np.isfinite(aux_model):
			return -np.inf
		aux_var = ratio_data["err"] ** 2
		loglike += -0.5 * ((ratio_data["value"] - aux_model) ** 2 / aux_var + np.log(2.0 * np.pi * aux_var))

	return loglike


def log_prior(theta, param_specs):
	lp = 0.0
	for value, spec in zip(theta, param_specs):
		prior_type = spec.get("prior_type", "uniform").lower()
		if prior_type == "none":
			continue
		if prior_type == "fixed":
			fixed_value = spec.get("fixed_value")
			if fixed_value is None:
				raise ValueError(f"Fixed prior requires fixed_value for {spec.get('name', '?')}")
			if abs(value - float(fixed_value)) > 1.0e-12:
				return -np.inf
			continue
		if prior_type == "uniform":
			low, high = spec.get("low"), spec.get("high")
			if low is None or high is None:
				raise ValueError(f"Uniform prior requires low/high for {spec.get('name', '?')}")
			if value < low or value > high:
				return -np.inf
			continue
		if prior_type == "gaussian":
			mean, std = spec.get("mean"), spec.get("std")
			if mean is None or std is None:
				raise ValueError(f"Gaussian prior requires mean/std for {spec.get('name', '?')}")
			if std <= 0.0:
				raise ValueError(f"Gaussian prior std must be > 0 for {spec.get('name', '?')}")
			z = (value - mean) / std
			lp += -0.5 * z * z - np.log(std) - 0.5 * np.log(2.0 * np.pi)
			continue
		if prior_type == "custom":
			x_grid = spec.get("x_grid")
			pdf_grid = spec.get("pdf_grid")
			low = spec.get("low")
			high = spec.get("high")
			if x_grid is None or pdf_grid is None:
				raise ValueError(f"Custom prior for {spec.get('name', '?')} is missing prepared KDE metadata")
			if low is not None and value < float(low):
				return -np.inf
			if high is not None and value > float(high):
				return -np.inf
			pdf = float(np.interp(float(value), x_grid, pdf_grid, left=0.0, right=0.0))
			if not np.isfinite(pdf) or pdf <= 0.0:
				return -np.inf
			lp += np.log(pdf)
			continue
		raise ValueError(f"Unknown prior_type: {prior_type}")
	return lp


def log_probability(
	theta,
	t,
	flux,
	yerr,
	param_specs,
	baseline_type="none",
	spline_knots=None,
	n_spline_knots=5,
	ratio_data=None,
	ratio_grid=None,
	flare_quiescent_flux=None,
):
	lp = log_prior(theta, param_specs)
	if not np.isfinite(lp):
		return -np.inf
	ll = log_likelihood(
		theta,
		t,
		flux,
		yerr,
		baseline_type=baseline_type,
		spline_knots=spline_knots,
		n_spline_knots=n_spline_knots,
		ratio_data=ratio_data,
		ratio_grid=ratio_grid,
		flare_quiescent_flux=flare_quiescent_flux,
	)
	if not np.isfinite(ll):
		return -np.inf
	return lp + ll


def initialize_walkers(nwalkers, param_specs, rng):
	def _draw_from_custom_prior(spec, size):
		cdf_grid = spec.get("cdf_grid")
		cdf_x_grid = spec.get("cdf_x_grid")
		if cdf_grid is None or cdf_x_grid is None:
			raise ValueError(f"Custom prior for {spec.get('name', '?')} is missing prepared KDE metadata")
		u = rng.uniform(0.0, 1.0, size=int(size))
		return np.interp(u, cdf_grid, cdf_x_grid)

	ndim = len(param_specs)
	p0 = np.zeros((nwalkers, ndim), dtype=float)
	for j, spec in enumerate(param_specs):
		prior_type = spec.get("prior_type", "uniform").lower()
		if prior_type == "fixed":
			fixed_value = spec.get("fixed_value")
			if fixed_value is None:
				raise ValueError(f"Fixed prior requires fixed_value for {spec.get('name', '?')}")
			p0[:, j] = float(fixed_value)
		elif prior_type == "uniform":
			p0[:, j] = rng.uniform(float(spec["low"]), float(spec["high"]), size=nwalkers)
		elif prior_type == "gaussian":
			p0[:, j] = rng.normal(float(spec["mean"]), float(spec["std"]), size=nwalkers)
		elif prior_type == "custom":
			p0[:, j] = _draw_from_custom_prior(spec, nwalkers)
		else:
			p0[:, j] = rng.normal(0.0, 1.0, size=nwalkers)
	return p0


def summarize_posterior(samples, names):
	q16, q50, q84 = np.quantile(samples, [0.16, 0.50, 0.84], axis=0)
	print("\nPosterior summary (median +/- 1 sigma):")
	for i, name in enumerate(names):
		print(f"  {name:>5s} = {q50[i]: .4f} (+{q84[i] - q50[i]: .4f}, -{q50[i] - q16[i]: .4f})")


def weighted_median_params(samples):
	return np.median(samples, axis=0)


def make_corner_plot(samples, names, truths=None, output_path="mcmc_corner.png"):
	try:
		import matplotlib.pyplot as plt
		import corner
	except ImportError:
		print("\nCorner plot not created: install plotting deps with `pip install corner matplotlib`.")
		return

	fig = corner.corner(
		samples,
		labels=names,
		truths=truths,
		quantiles=[0.16, 0.50, 0.84],
		show_titles=True,
		title_fmt=".3f",
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
	baseline_type="none",
	n_spline_knots=5,
	spline_knots=None,
	output_path="mcmc_data_model.png",
	n_draws=80,
	seed=11,
	ratio_data=None,
	ratio_grid=None,
	flare_quiescent_flux=None,
):
	try:
		import matplotlib.pyplot as plt
	except ImportError:
		print("\nData/model plot not created: install matplotlib with `pip install matplotlib`.")
		return

	tgrid = np.linspace(np.min(t), np.max(t), 600)
	rng = np.random.default_rng(seed)
	n = min(n_draws, samples.shape[0])
	draw_idx = rng.choice(samples.shape[0], size=n, replace=False)

	if ratio_data is not None and ratio_grid is not None:
		fig, axes = plt.subplots(1, 3, figsize=(12.6, 5.0), gridspec_kw={"width_ratios": [3.3, 1.7, 1.4]})
		ax, aux_ax, ratio_ax = axes
	else:
		fig, ax = plt.subplots(figsize=(7.5, 5.0))
		aux_ax, ratio_ax = None, None

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

	theta_best = weighted_median_params(samples)
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
		label="Best-fit (median)",
	)

	tpeak_best, fwhm_best, ampl_best = theta_best[:3]
	flare_only = 1.0 + aflare1(tgrid, tpeak_best, fwhm_best, ampl_best)
	ax.plot(tgrid, flare_only, color="tab:orange", lw=1.8, ls="--", label="Flare only")

	ax.errorbar(t, flux, yerr=yerr, fmt="o", ms=4, alpha=0.85, color="black", ecolor="gray", elinewidth=1.0, capsize=2, label="Observed data")

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
				aux_ax.scatter(aux_predictions, jitter, s=18, color="tab:blue", alpha=0.25, edgecolors="none", label="Posterior samples")
				q16, q50, q84 = np.quantile(aux_predictions, [0.16, 0.50, 0.84])
				aux_ax.axvspan(q16, q84, color="tab:blue", alpha=0.10)
				aux_ax.axvline(q50, color="tab:blue", lw=1.5, ls="--", label="Sample median")

			aux_ax.axvspan(ratio_data["value"] - ratio_data["err"], ratio_data["value"] + ratio_data["err"], color="tab:red", alpha=0.15)
			aux_ax.axvline(ratio_data["value"], color="tab:red", lw=2.0, label="Observed point")
			if np.isfinite(best_aux_value):
				aux_ax.axvline(best_aux_value, color="black", lw=1.8, ls="-.", label="Best-fit prediction")

			aux_ax.set_ylim(-0.2, 0.2)
			aux_ax.set_yticks([])
			aux_ax.set_xlabel("Auxiliary Flux (absolute)" if ratio_data.get("input_mode") == "absolute" else "Auxiliary Value")
			aux_ax.set_title("Auxiliary Posterior")
			aux_ax.grid(alpha=0.25, axis="x")
			aux_ax.legend(loc="best", fontsize=8)

			if ratio_ax is not None and ratio_predictions.size > 0:
				ratio_jitter = rng.uniform(-0.12, 0.12, size=ratio_predictions.size)
				ratio_ax.scatter(ratio_predictions, ratio_jitter, s=18, color="tab:purple", alpha=0.28, edgecolors="none", label="Posterior samples")
				rq16, rq50, rq84 = np.quantile(ratio_predictions, [0.16, 0.50, 0.84])
				ratio_ax.axvspan(rq16, rq84, color="tab:purple", alpha=0.10)
				ratio_ax.axvline(rq50, color="tab:purple", lw=1.5, ls="--", label="Sample median")
				if np.isfinite(best_ratio_value):
					ratio_ax.axvline(best_ratio_value, color="black", lw=1.8, ls="-.", label="Best-fit")

				flare_model_at_aux = float(aflare1(np.array([ratio_data["time"]]), tpeak_best, fwhm_best, ampl_best)[0])
				if ratio_data.get("input_mode") == "absolute" and flare_quiescent_flux is not None:
					flare_model_at_aux *= float(flare_quiescent_flux)

				if np.isfinite(flare_model_at_aux) and abs(flare_model_at_aux) > 0.0:
					ratio_obs = ratio_data["value"] / flare_model_at_aux
					ratio_obs_err = ratio_data["err"] / abs(flare_model_at_aux)
					if np.isfinite(ratio_obs) and np.isfinite(ratio_obs_err):
						ratio_ax.axvspan(ratio_obs - ratio_obs_err, ratio_obs + ratio_obs_err, color="tab:red", alpha=0.12)
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


def run_mcmc(
	t,
	flux,
	yerr,
	param_specs,
	baseline_type="none",
	n_spline_knots=5,
	spline_knots=None,
	nwalkers=64,
	nsteps=2000,
	nburn=500,
	n_cores=1,
	ratio_data=None,
	ratio_grid=None,
	flare_quiescent_flux=None,
	seed=123,
):
	tmin = float(np.min(t))
	tmax = float(np.max(t))
	if baseline_type == "spline" and spline_knots is None:
		spline_knots = np.linspace(tmin, tmax, n_spline_knots)

	ndim = len(param_specs)
	nwalkers = max(int(nwalkers), 2 * ndim + 2)
	nsteps = int(nsteps)
	nburn = max(0, min(int(nburn), nsteps - 1))

	rng = np.random.default_rng(seed)
	p0 = initialize_walkers(nwalkers, param_specs, rng)

	sampler_kwargs = dict(
		nwalkers=nwalkers,
		ndim=ndim,
		log_prob_fn=log_probability,
		args=(
			t,
			flux,
			yerr,
			param_specs,
			baseline_type,
			spline_knots,
			n_spline_knots,
			ratio_data,
			ratio_grid,
			flare_quiescent_flux,
		),
	)

	n_cores = max(1, int(n_cores))
	if n_cores == 1:
		sampler = emcee.EnsembleSampler(**sampler_kwargs)
		sampler.run_mcmc(p0, nsteps, progress=True)
	else:
		with mp.Pool(processes=n_cores) as pool:
			sampler = emcee.EnsembleSampler(pool=pool, **sampler_kwargs)
			sampler.run_mcmc(p0, nsteps, progress=True)

	chain = sampler.get_chain(discard=nburn, flat=True)
	log_prob = sampler.get_log_prob(discard=nburn, flat=True)
	return sampler, chain, log_prob


def truth_vector_from_names(names, truth):
	truths = []
	spline_coeffs = truth.get("spline_coeffs", None)
	for name in names:
		if name.startswith("spline_c") and spline_coeffs is not None:
			idx = int(name.replace("spline_c", ""))
			truths.append(float(spline_coeffs[idx]))
		elif name in truth:
			truths.append(float(truth[name]))
		else:
			truths.append(np.nan)
	return truths


def load_settings_csv(filepath):
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

	nwalkers = row.get("nwalkers", row.get("nlive", 64))
	nsteps = row.get("nsteps", row.get("maxiter", 2000))
	nburn = row.get("nburn", 500)

	return {
		"baseline_model": row.get("baseline_model", "none"),
		"n_spline_knots": int(row.get("n_spline_knots", 5)),
		"nwalkers": int(nwalkers) if pd.notna(nwalkers) else 64,
		"nsteps": int(nsteps) if pd.notna(nsteps) else 2000,
		"nburn": int(nburn) if pd.notna(nburn) else 500,
		"n_cores": max(1, int(row.get("n_cores", 1))) if pd.notna(row.get("n_cores", 1)) else 1,
		"flare_data_csv": row.get("flare_data_csv", None) if pd.notna(row.get("flare_data_csv", None)) else None,
		"use_ratio_model": maybe_bool(row.get("use_ratio_model", False), default=False),
		"ratio_filter": normalize_ratio_filter(row.get("ratio_filter", "F087")),
		"ratio_grid_csv": row.get("ratio_grid_csv", None) if pd.notna(row.get("ratio_grid_csv", None)) else None,
		"ratio_data_csv": row.get("ratio_data_csv", None) if pd.notna(row.get("ratio_data_csv", None)) else None,
		"aux_quiescent_flux": float(row.get("aux_quiescent_flux")) if pd.notna(row.get("aux_quiescent_flux")) else None,
		"flare_quiescent_flux": float(row.get("flare_quiescent_flux")) if pd.notna(row.get("flare_quiescent_flux")) else None,
	}


def load_params_csv(filepath):
	df = pd.read_csv(filepath)
	param_specs = []
	for idx, row in df.iterrows():
		spec = {"name": row["name"], "prior_type": row["prior_type"].strip().lower()}
		if spec["prior_type"] not in ["none", "uniform", "gaussian", "fixed", "custom"]:
			raise ValueError(f"Unknown prior_type: {spec['prior_type']} at row {idx}")
		if spec["prior_type"] == "uniform":
			spec["low"] = float(row["low"])
			spec["high"] = float(row["high"])
		elif spec["prior_type"] == "gaussian":
			spec["mean"] = float(row["mean"])
			spec["std"] = float(row["std"])
		elif spec["prior_type"] == "custom":
			custom_kind = "kde"
			if "custom_kind" in row.index and pd.notna(row["custom_kind"]):
				custom_kind = str(row["custom_kind"]).strip().lower()
			elif "prior_kind" in row.index and pd.notna(row["prior_kind"]):
				custom_kind = str(row["prior_kind"]).strip().lower()
			spec["custom_kind"] = custom_kind
			if "bandwidth" in row.index and pd.notna(row["bandwidth"]):
				spec["bandwidth"] = float(row["bandwidth"])
			if "grid_size" in row.index and pd.notna(row["grid_size"]):
				spec["grid_size"] = int(row["grid_size"])
			if "low" in row.index and pd.notna(row["low"]):
				spec["low"] = float(row["low"])
			if "high" in row.index and pd.notna(row["high"]):
				spec["high"] = float(row["high"])
		elif spec["prior_type"] == "fixed":
			fixed_value = None
			if "fixed_value" in row.index and pd.notna(row["fixed_value"]):
				fixed_value = float(row["fixed_value"])
			elif "mean" in row.index and pd.notna(row["mean"]):
				fixed_value = float(row["mean"])
			elif "low" in row.index and pd.notna(row["low"]):
				fixed_value = float(row["low"])
			elif "high" in row.index and pd.notna(row["high"]):
				fixed_value = float(row["high"])
			if fixed_value is None:
				raise ValueError(
					f"Fixed prior requires fixed_value (or mean/low/high fallback) at row {idx}"
				)
			spec["fixed_value"] = fixed_value
		param_specs.append(spec)
	return param_specs


def save_settings_csv(filepath, settings_dict):
	pd.DataFrame([settings_dict]).to_csv(filepath, index=False)
	print(f"Settings CSV saved to: {filepath}")


def save_params_csv(filepath, param_specs):
	rows = []
	for spec in param_specs:
		rows.append(
			{
				"name": spec["name"],
				"prior_type": spec["prior_type"],
				"low": spec.get("low", ""),
				"high": spec.get("high", ""),
				"mean": spec.get("mean", ""),
				"std": spec.get("std", ""),
				"fixed_value": spec.get("fixed_value", ""),
				"bandwidth": spec.get("bandwidth", ""),
				"grid_size": spec.get("grid_size", ""),
				"custom_kind": spec.get("custom_kind", ""),
			}
		)
	pd.DataFrame(rows).to_csv(filepath, index=False)
	print(f"Parameters CSV saved to: {filepath}")


def mcmc_pipeline(
	flare_data_csv,
	params_csv=None,
	settings_csv=None,
	custom_prior_samples=None,
	ratio_data_csv=None,
	ratio_grid_csv=None,
	baseline_model="none",
	n_spline_knots=5,
	nwalkers=64,
	nsteps=2000,
	nburn=500,
	n_cores=1,
	use_ratio_model=False,
	ratio_filter="F087",
	aux_quiescent_flux=None,
	flare_quiescent_flux=None,
	output_dir=".",
	plot_prefix="mcmc",
	update_settings_csv=False,
	update_params_csv=False,
	seed=123,
):
	if settings_csv is not None:
		settings = load_settings_csv(settings_csv)
		baseline_model = settings.get("baseline_model", baseline_model)
		n_spline_knots = settings.get("n_spline_knots", n_spline_knots)
		nwalkers = settings.get("nwalkers", nwalkers)
		nsteps = settings.get("nsteps", nsteps)
		nburn = settings.get("nburn", nburn)
		n_cores = settings.get("n_cores", n_cores)
		flare_data_csv = settings.get("flare_data_csv", flare_data_csv)
		use_ratio_model = settings.get("use_ratio_model", use_ratio_model)
		ratio_filter = settings.get("ratio_filter", ratio_filter)
		if settings.get("ratio_grid_csv") is not None:
			ratio_grid_csv = settings["ratio_grid_csv"]
		if settings.get("ratio_data_csv") is not None:
			ratio_data_csv = settings["ratio_data_csv"]
		if aux_quiescent_flux is None and settings.get("aux_quiescent_flux") is not None:
			aux_quiescent_flux = settings["aux_quiescent_flux"]
		if flare_quiescent_flux is None and settings.get("flare_quiescent_flux") is not None:
			flare_quiescent_flux = settings["flare_quiescent_flux"]

	ratio_grid = None
	if use_ratio_model:
		ratio_grid_path = resolve_ratio_grid_path(ratio_filter=ratio_filter, ratio_grid_csv=ratio_grid_csv)
		ratio_grid = load_ratio_grid(ratio_grid_path)

	if flare_data_csv is None:
		raise ValueError("flare_data_csv is required")
	t, flux, yerr = load_flare_data_csv(flare_data_csv)
	tmin, tmax = float(np.min(t)), float(np.max(t))

	if params_csv is not None:
		param_specs = load_params_csv(params_csv)
	else:
		param_specs = default_param_specs(
			tmin,
			tmax,
			baseline_type=baseline_model,
			n_spline_knots=n_spline_knots,
			use_ratio_model=use_ratio_model,
			ratio_grid=ratio_grid,
		)

	if use_ratio_model and ratio_grid is not None:
		param_specs = ensure_ratio_param_specs(param_specs, ratio_grid)

	param_specs = prepare_param_specs(param_specs, custom_prior_samples=custom_prior_samples)

	ratio_data = None
	if use_ratio_model:
		if ratio_data_csv is None:
			raise ValueError(
				"ratio_data_csv is required when use_ratio_model=true. Provide it in settings.csv or via argument."
			)
		ratio_data = load_ratio_data_csv(
			ratio_data_csv,
			aux_quiescent_flux=aux_quiescent_flux,
			flare_quiescent_flux=flare_quiescent_flux,
		)

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
		f"Running MCMC with baseline={baseline_model}, nwalkers={nwalkers}, "
		f"nsteps={nsteps}, nburn={nburn}, n_cores={n_cores}, use_ratio_model={use_ratio_model}"
	)

	if update_settings_csv and settings_csv is not None:
		settings_dict = {
			"baseline_model": baseline_model,
			"n_spline_knots": n_spline_knots,
			"nwalkers": nwalkers,
			"nsteps": nsteps,
			"nburn": nburn,
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

	if update_params_csv and params_csv is not None:
		save_params_csv(params_csv, param_specs)

	sampler, samples, log_prob = run_mcmc(
		t,
		flux,
		yerr,
		param_specs=param_specs,
		baseline_type=baseline_model,
		n_spline_knots=n_spline_knots,
		spline_knots=spline_knots,
		nwalkers=nwalkers,
		nsteps=nsteps,
		nburn=nburn,
		n_cores=n_cores,
		ratio_data=ratio_data,
		ratio_grid=ratio_grid,
		flare_quiescent_flux=flare_quiescent_flux,
		seed=seed,
	)

	names = [spec["name"] for spec in param_specs]
	plot_indices = [
		i for i, spec in enumerate(param_specs)
		if spec.get("prior_type", "uniform").lower() != "fixed"
	]
	plot_names = [names[i] for i in plot_indices]
	plot_samples = samples[:, plot_indices] if len(plot_indices) > 0 else np.empty((samples.shape[0], 0))

	if len(plot_names) > 0:
		truths = truth_vector_from_names(plot_names, truth)
		summarize_posterior(plot_samples, names=plot_names)
		make_corner_plot(
			plot_samples,
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
		baseline_type=baseline_model,
		n_spline_knots=n_spline_knots,
		spline_knots=spline_knots,
		output_path=f"{output_dir}/{plot_prefix}_data_model_{baseline_model}.png",
		n_draws=80,
		ratio_data=ratio_data,
		ratio_grid=ratio_grid,
		flare_quiescent_flux=flare_quiescent_flux,
	)

	return sampler, samples, log_prob, names


def parse_args():
	parser = argparse.ArgumentParser(description="emcee MCMC for aflare1 flare fitting")
	parser.add_argument("--baseline-model", choices=["none", "sine", "spline"], default="none")
	parser.add_argument("--n-spline-knots", type=int, default=5)
	parser.add_argument("--nwalkers", type=int, default=64)
	parser.add_argument("--nsteps", type=int, default=2000)
	parser.add_argument("--nburn", type=int, default=500)
	parser.add_argument("--n-cores", type=int, default=1)
	parser.add_argument("--settings-csv", type=str, default=None)
	parser.add_argument("--params-csv", type=str, default=None)
	parser.add_argument("--flare-data-csv", type=str, default=None)
	parser.add_argument("--use-ratio-model", action="store_true")
	parser.add_argument("--ratio-filter", choices=sorted(RATIO_GRID_FILES), default="F087")
	parser.add_argument("--ratio-grid-csv", type=str, default=None)
	parser.add_argument("--ratio-data-csv", type=str, default=None)
	parser.add_argument("--aux-quiescent-flux", type=float, default=None)
	parser.add_argument("--flare-quiescent-flux", type=float, default=None)
	parser.add_argument("--seed", type=int, default=123)
	return parser.parse_args()


def main():
	args = parse_args()
	mcmc_pipeline(
		flare_data_csv=args.flare_data_csv,
		params_csv=args.params_csv,
		settings_csv=args.settings_csv,
		ratio_data_csv=args.ratio_data_csv,
		ratio_grid_csv=args.ratio_grid_csv,
		baseline_model=args.baseline_model,
		n_spline_knots=args.n_spline_knots,
		nwalkers=args.nwalkers,
		nsteps=args.nsteps,
		nburn=args.nburn,
		n_cores=args.n_cores,
		use_ratio_model=args.use_ratio_model,
		ratio_filter=args.ratio_filter,
		aux_quiescent_flux=args.aux_quiescent_flux,
		flare_quiescent_flux=args.flare_quiescent_flux,
		seed=args.seed,
		update_settings_csv=False,
		update_params_csv=False,
	)
