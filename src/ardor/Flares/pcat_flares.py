"""Variable-count flare inference with PCAT."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ardor.paths import get_repository_path
from tdpy.verbosity import print


def _log_path(path):
    try:
        return Path(path).resolve().relative_to(get_repository_path()).as_posix()
    except ValueError:
        return Path(path).name


@dataclass
class PCATFlareFit:
    """PCAT posterior and the flux/time scaling needed to interpret it."""

    posterior: object
    time_days: np.ndarray
    flux: np.ndarray
    uncertainty: np.ndarray
    bin_edges: np.ndarray
    time_offset_days: float
    reference_bjd: float | None
    baseline_flux: float
    count_scale: np.ndarray
    run_directory: Path


def _time_bin_edges(time_days: np.ndarray) -> tuple[np.ndarray, float]:
    from tdpy.numerics import midpoint_bin_edges

    if np.any(np.diff(time_days) <= 0.0):
        raise ValueError("time_days must be strictly increasing")
    edges = midpoint_bin_edges(time_days)
    offset = 1.0 - edges[0]  # [day], keeps PCAT's log-grid time axis positive
    return edges + offset, offset


def _flare_fit_configuration(time_days, flux, uncertainty, baseline_flux=None):
    if any(np.ma.is_masked(values) for values in (time_days, flux, uncertainty)):
        raise ValueError("time, flux, and uncertainty must not contain masked observations")
    time_days = np.asarray(time_days, dtype=float)
    flux = np.asarray(flux, dtype=float)
    uncertainty = np.asarray(uncertainty, dtype=float)
    if (time_days.ndim != 1 or flux.shape != time_days.shape
            or uncertainty.shape != time_days.shape or time_days.size < 5):
        raise ValueError("time, flux, and uncertainty must be matching 1D arrays with at least five points")
    if (not np.isfinite(time_days).all() or not np.isfinite(flux).all()
            or not np.isfinite(uncertainty).all() or np.any(flux <= 0.0)
            or np.any(uncertainty <= 0.0)):
        raise ValueError("flux must be positive, and all inputs must be finite with positive uncertainties")

    bin_edges, time_offset = _time_bin_edges(time_days)
    bin_widths = np.diff(bin_edges)
    time_span = float(bin_edges[-1] - bin_edges[0])  # [day]
    cadence = float(np.median(np.diff(time_days)))  # [day]
    if baseline_flux is None:
        baseline_flux = float(np.quantile(flux, 0.2))
    baseline_flux = float(baseline_flux)
    if not np.isfinite(baseline_flux) or baseline_flux <= 0.0:
        raise ValueError("baseline_flux must be finite and positive")

    median_uncertainty = float(np.median(uncertainty))
    amplitude_minimum = max(0.2 * median_uncertainty, 1e-10)
    peak_excess = max(float(np.max(flux - baseline_flux)), median_uncertainty)
    amplitude_maximum = max(2.0 * peak_excess, 10.0 * median_uncertainty)
    rise_minimum = max(0.5 * cadence, 1e-6)  # [day]
    rise_maximum = max(rise_minimum * 2.0, min(0.25 * time_span, 20.0 * cadence))  # [day]
    decay_minimum = max(cadence, rise_minimum)  # [day]
    decay_maximum = max(decay_minimum * 2.0, min(0.9 * time_span, 80.0 * cadence))  # [day]
    limits = {
        "flux": (amplitude_minimum, amplitude_maximum),
        "scalrise": (rise_minimum, rise_maximum),
        "scalfall": (decay_minimum, decay_maximum),
    }
    count_scale = baseline_flux / uncertainty**2
    exposure = count_scale / bin_widths
    return {
        "time_days": time_days,
        "flux": flux,
        "uncertainty": uncertainty,
        "bin_edges": bin_edges,
        "time_offset": time_offset,
        "baseline_flux": baseline_flux,
        "count_scale": count_scale,
        "exposure": exposure,
        "limits": limits,
        "cadence": cadence,
    }


def fit_variable_flares_pcat(
    time_days,
    flux,
    uncertainty,
    run_directory,
    *,
    reference_bjd=None,
    baseline_flux=None,
    numbswep=20_000,
    max_flares=10,
    numbproc=4,
    seed=0,
):
    """Fit an unknown number of FRED flares to an Ardor light-curve segment.

    Flux and uncertainty use the same units, and times are relative days. A
    fixed quiescent baseline is inferred from the lower flux quantile unless
    supplied. The observed errors are encoded as equivalent Poisson exposure
    weights, allowing PCAT's transdimensional flare sampler to use Ardor's
    heteroscedastic light-curve uncertainties.
    """

    try:
        from pcat import sampling
        from pcat.time_series import evaluate_flare_profile, retr_dictpropelemtmpl
    except ImportError as exception:
        raise ImportError(
            "PCAT fitting requires the optional Ardor extra: pip install -e '.[pcat]'"
        ) from exception

    configuration = _flare_fit_configuration(
        time_days, flux, uncertainty, baseline_flux
    )
    if numbswep < 20 or max_flares < 1 or numbproc < 1:
        raise ValueError("numbswep must be at least 20; max_flares and numbproc must be positive")

    run_directory = Path(run_directory)
    if run_directory.exists() and any(run_directory.iterdir()):
        raise FileExistsError(f"PCAT run directory is not empty: {run_directory}")
    input_directory = run_directory / "data" / "inpt"
    input_directory.mkdir(parents=True, exist_ok=True)

    import astropy.io.fits

    count_scale = configuration["count_scale"]
    input_flux = np.rint(configuration["flux"] * count_scale) / count_scale
    for filename, values in (
        (
            "sbrt.fits",
            input_flux[:, None, None, None],
        ),
        (
            "expo.fits",
            configuration["exposure"][:, None, None],
        ),
    ):
        path = input_directory / filename
        print(f"Writing to {_log_path(path)}...")
        astropy.io.fits.writeto(path, np.asarray(values, dtype=float), overwrite=True)

    time_days = configuration["time_days"]
    time_edges = configuration["bin_edges"]
    time_centers = time_days + configuration["time_offset"]
    observed_counts = np.rint(configuration["flux"] * count_scale)
    baseline_counts = configuration["baseline_flux"] * count_scale
    limits = configuration["limits"]
    number_templates = 512
    peak_times = np.geomspace(time_edges[0], time_edges[-1], number_templates)  # [day]
    ones = np.ones(number_templates)
    rise_template = np.full(number_templates, np.sqrt(np.prod(limits["scalrise"])))  # [day]
    decay_template = np.full(number_templates, np.sqrt(np.prod(limits["scalfall"])))  # [day]
    unit_profiles = evaluate_flare_profile(
        time_centers,
        "flarfred",
        ones,
        peak_times,
        rise_time=rise_template,
        decay_time=decay_template,
    )
    birth_proposal = retr_dictpropelemtmpl(
        observed_counts - baseline_counts,
        np.maximum(observed_counts, 1.0),
        unit_profiles.T * count_scale[None, :],
        *limits["flux"],
        index_position=1,
        index_flux=0,
        number_parameters=4,
    )

    number_burn = 3 * numbswep // 4
    number_samples = min(max(numbswep // 40, 16), 2000)
    sampling_result = sampling.sample(
        typeexpr="fire",
        typedata="inpt",
        strgexprsbrt="sbrt.fits",
        typeexpo="file",
        strgexpo="expo.fits",
        binsenerfull=time_edges,
        spectype=["flarfred"],
        spatdisttype=["line"],
        typeelem=["lghtlinevoig"],
        dictfitt={
            "typeelem": ["lghtlinevoig"],
            "spectype": ["flarfred"],
            "sbrtbacknorm":[
                np.full((time_days.size, 1, 1), configuration["baseline_flux"])
            ],
            "listnamediff": ["back0000"],
        },
        limtparaelem=limits,
        maxmgangdata=100.0 / (3600.0 * 180.0 / np.pi),  # [rad], unused for time-series data
        anlytype="spec",
        fittminmnumbelempop0=0,
        fittmaxmnumbelempop0=max_flares,
        inittype="rand",
        typeseed=seed,
        numbproc=numbproc,
        probtran=0.5,
        probspmr=0.3,
        probjump=0.3,
        radispmr=max(0.01, configuration["cadence"]),  # [day]
        stdvpropelemfire=[
            0.05 * limits["flux"][1],
            configuration["cadence"],
            0.1 * (limits["scalrise"][1] - limits["scalrise"][0]),
            0.1 * (limits["scalfall"][1] - limits["scalfall"][0]),
        ],
        booladaptstdp=True,
        boolburntmpr=True,
        factburntmpr=0.9,
        numbswep=numbswep,
        numbburn=number_burn,
        numbsamp=number_samples,
        numbswepplot=max(numbswep // 20, 1),
        **birth_proposal,
        makeanim=False,
        boolmakeplot=False,
        boolmakeplotinit=False,
        booldiag=False,
        typeverb=0,
        pathbase=str(run_directory),
        strgcnfg=run_directory.name,
    )
    from pcat import main as pcat_main

    posterior_path = (
        run_directory
        / "data"
        / "outp"
        / run_directory.name
        / "gdatfinlpost"
    )
    print(f"Reading from {_log_path(posterior_path)}...")
    posterior = pcat_main.readfile(str(posterior_path))
    return PCATFlareFit(
        posterior=posterior,
        time_days=time_days,
        flux=configuration["flux"],
        uncertainty=configuration["uncertainty"],
        bin_edges=time_edges,
        time_offset_days=configuration["time_offset"],
        reference_bjd=reference_bjd,
        baseline_flux=configuration["baseline_flux"],
        count_scale=count_scale,
        run_directory=run_directory,
    )


def fit_tier2_csv_pcat(target_file, run_directory, **fit_options):
    """Fit one Ardor Tier 2 flare CSV with PCAT's variable-count sampler."""

    target_file = Path(target_file)
    print(f"Reading from {_log_path(target_file)}...")
    with target_file.open(encoding="utf-8") as file_object:
        reference_bjd = float(file_object.readline().strip().lstrip("# "))
    data = np.loadtxt(target_file, delimiter=",", comments="#", ndmin=2)
    if data.shape[1] != 3:
        raise ValueError("Tier 2 flare CSV must contain time [day], flux, and uncertainty columns")
    return fit_variable_flares_pcat(
        data[:, 0],
        data[:, 1],
        data[:, 2],
        run_directory,
        reference_bjd=reference_bjd,
        **fit_options,
    )


def tier3_pcat(tier_2_output_dir, tier_3_output_dir, host_name="My_Host",
               numbswep=20_000, max_flares=10, numbproc=4):
    """Fit every Ardor Tier 2 flare cutout with PCAT's variable-count model."""

    tier_2_output_dir = Path(tier_2_output_dir)
    if not tier_2_output_dir.is_dir():
        raise NotADirectoryError(tier_2_output_dir)
    flare_csvs = sorted(
        tier_2_output_dir.glob("*.csv"),
        key=lambda path: path.name,
    )
    if not flare_csvs:
        raise FileNotFoundError(f"No Tier 2 flare CSV files found in {tier_2_output_dir}")

    output_root = Path(tier_3_output_dir) / str(host_name)
    fits = []
    for flare_csv in flare_csvs:
        fits.append(
            fit_tier2_csv_pcat(
                flare_csv,
                output_root / flare_csv.stem,
                numbswep=numbswep,
                max_flares=max_flares,
                numbproc=numbproc,
            )
        )
    return fits


def posterior_flux_quantiles(fit: PCATFlareFit) -> np.ndarray:
    """Return 16th, 50th, and 84th percentile model fluxes at each cadence."""

    models = np.asarray(fit.posterior.listpostcntpmodl, dtype=float)[:, :, 0, 0]
    return np.percentile(models / fit.count_scale[None, :], [16.0, 50.0, 84.0], axis=0)