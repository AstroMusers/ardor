import numpy as np
import pytest

from ardor.Roman.flare_model import (
    baseline_model,
    ensure_ratio_param_specs,
    infer_ratio_grid_bounds,
    load_flare_data_csv,
    load_ratio_data_csv,
    load_ratio_grid,
    log_likelihood,
    model_flux,
    normalize_ratio_filter,
    parameter_names,
    resolve_ratio_grid_path,
    split_theta,
    truth_vector_from_names,
)


def test_shared_flare_csv_loader(tmp_path, capsys):
    filepath = tmp_path / "flare.csv"
    filepath.write_text("flux,yerr,time,extra\n1.0,0.01,0.0,a\n1.2,0.02,0.5,b\n", encoding="utf-8")
    time, flux, uncertainty = load_flare_data_csv(filepath)
    np.testing.assert_array_equal(time, [0.0, 0.5])
    np.testing.assert_array_equal(flux, [1.0, 1.2])
    np.testing.assert_array_equal(uncertainty, [0.01, 0.02])
    assert capsys.readouterr().out == f"Reading from {filepath}...\n"


@pytest.mark.parametrize("contents, message", [
    ("time,flux\n0,1\n1,1\n", "missing required columns: yerr"),
    ("time,flux,yerr\n0,1,0.01\n", "at least 2 rows"),
])
def test_shared_flare_csv_loader_rejects_invalid_input(tmp_path, contents, message):
    filepath = tmp_path / "flare.csv"
    filepath.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_flare_data_csv(filepath)


@pytest.mark.parametrize("sampler_name, dependency", [("mcmc", "emcee"), ("ns", "dynesty")])
def test_roman_samplers_reuse_upstream_kde_and_shared_csv(sampler_name, dependency):
    from importlib import import_module
    from ardor.Flares.aflare import aflare1
    from tdpy.numerics import build_kde_prior_from_samples

    pytest.importorskip(dependency)
    sampler = import_module(f"ardor.Roman.{sampler_name}")
    assert sampler.aflare1 is aflare1
    assert sampler.load_flare_data_csv is load_flare_data_csv
    samples = np.array([-0.5, -0.1, 0.2, 1.0, np.nan])
    normalization = "trapezoid" if sampler_name == "mcmc" else "rectangle"
    expected = build_kde_prior_from_samples(
        samples, low=-1.0, high=2.0, grid_size=256, normalization=normalization,
    )
    if sampler_name == "mcmc":
        expected.pop("n_samples")
    else:
        expected["x_grid"] = expected.pop("cdf_x_grid")
        expected.pop("pdf_grid")
    actual = sampler.build_kde_prior_from_samples(samples, low=-1.0, high=2.0, grid_size=256)
    assert actual.keys() == expected.keys()
    for key in actual:
        if isinstance(expected[key], np.ndarray):
            np.testing.assert_array_equal(actual[key], expected[key])
        else:
            assert actual[key] == expected[key]

    specs = [{"name": "ampl", "prior_type": "custom", "low": -1.0, "high": 2.0}]
    prepared = sampler.prepare_param_specs(specs, custom_prior_samples={"ampl": samples})
    assert "x_grid" in prepared[0]
    assert "x_grid" not in specs[0]


@pytest.mark.parametrize("sampler_name, dependency", [("mcmc", "emcee"), ("ns", "dynesty")])
def test_roman_prior_overrides_preserve_sampler_specific_metadata(sampler_name, dependency):
    from importlib import import_module

    pytest.importorskip(dependency)
    sampler = import_module(f"ardor.Roman.{sampler_name}")
    specs = [{"name": "tpeak", "prior_type": "uniform", "low": 0.0, "high": 1.0},
             {"name": "ampl", "prior_type": " CUSTOM ", "prior_kind": "KDE",
              "low": -1.0, "high": 2.0, "bandwidth": 0.2, "grid_size": 512}]
    original = [dict(spec) for spec in specs]
    samples = np.array([-0.5, 0.0, 0.5, 1.0, np.nan])
    overrides = {"samples": samples, "low": -0.8, "high": 1.5, "bandwidth": 0.4, "grid_size": 256}
    prepared = sampler.prepare_param_specs(specs, custom_prior_samples={"ampl": overrides})
    expected = sampler.build_kde_prior_from_samples(samples, low=-0.8, high=1.5, bandwidth=0.4, grid_size=256)
    fields = {"bandwidth", "grid_size", "low", "high", "x_grid", "cdf_grid"}
    if sampler_name == "mcmc":
        fields.update({"pdf_grid", "cdf_x_grid"})
    assert prepared[1].keys() == specs[1].keys() | fields | {"custom_kind"}
    for name in fields:
        np.testing.assert_array_equal(prepared[1][name], expected[name])
    assert specs == original
    assert prepared[0] == specs[0] and prepared[0] is not specs[0]
    assert np.isnan(samples[-1])
    pipeline_name = "mcmc_pipeline" if sampler_name == "mcmc" else "nested_sampling_pipeline"
    with pytest.raises(ValueError, match=pipeline_name):
        sampler.prepare_custom_prior_spec({"name": "missing"})


def test_roman_posterior_summaries_use_shared_weighted_statistics(capsys):
    from importlib import import_module
    from tdpy.numerics import weighted_median, weighted_quantile

    pytest.importorskip("dynesty")
    sampler = import_module("ardor.Roman.ns")
    assert sampler._weighted_quantile is weighted_quantile
    assert sampler._weighted_median is weighted_median
    samples = np.array([[10.0, 1.0], [0.0, 4.0], [20.0, 2.0], [30.0, 8.0]])
    weights = np.array([1.0, 4.0, 2.0, 3.0])
    np.testing.assert_allclose(sampler.weighted_median_params(samples, weights), [10.0, 3.0])
    np.testing.assert_allclose(
        sampler.weighted_quantile(samples[:, 0], [0.16, 0.5, 0.84], sample_weight=weights),
        [0.0, 10.0, 20.0 + 10.0 * 0.14 / 0.3],
    )
    sampler.summarize_posterior(samples, weights, ["first", "second"])
    output = capsys.readouterr().out
    assert "10.0000" in output and "3.0000" in output


def test_roman_posterior_summaries_preserve_weight_validation():
    from importlib import import_module

    pytest.importorskip("dynesty")
    sampler = import_module("ardor.Roman.ns")
    with pytest.raises(ValueError, match="positive sum"):
        sampler.weighted_quantile([1.0, 2.0], [0.5], [0.0, 0.0])
    with pytest.raises(ValueError, match="non-negative"):
        sampler.weighted_median_params(np.array([[1.0], [2.0]]), np.array([1.0, -1.0]))


@pytest.mark.parametrize("sampler_name, dependency", [("mcmc", "emcee"), ("ns", "dynesty")])
@pytest.mark.parametrize("entry_point", ["module", "script"])
def test_roman_sampler_cli_without_import_path_mutation(tmp_path, sampler_name, dependency, entry_point):
    from importlib import import_module
    import subprocess
    import sys

    pytest.importorskip(dependency)
    module_name = f"ardor.Roman.{sampler_name}"
    sampler = import_module(module_name)
    arguments = ["-m", module_name] if entry_point == "module" else [sampler.__file__]
    result = subprocess.run(
        [sys.executable, *arguments, "--help"], cwd=tmp_path,
        capture_output=True, text=True, check=True,
    )
    if sampler_name == "ns":
        assert "usage:" in result.stdout
        assert "--help" in result.stdout
    else:
        assert result.stdout == ""


@pytest.mark.parametrize("sampler_name, dependency", [("mcmc", "emcee"), ("ns", "dynesty")])
def test_roman_sampler_argument_parser(tmp_path, sampler_name, dependency):
    import subprocess
    import sys

    pytest.importorskip(dependency)
    result = subprocess.run(
        [sys.executable, "-c", f"from ardor.Roman.{sampler_name} import parse_args; parse_args()", "--help"],
        cwd=tmp_path, capture_output=True, text=True, check=True,
    )
    assert "usage:" in result.stdout
    assert "--help" in result.stdout


def test_baseline_models_and_parameter_split():
    time = np.linspace(-1.0, 1.0, 11)
    np.testing.assert_array_equal(baseline_model(time), np.zeros_like(time))
    np.testing.assert_allclose(
        baseline_model(time, "sine", [0.1, 0.02, 1.0, 0.0]),
        0.1 + 0.02 * np.sin(2.0 * np.pi * time),
    )
    assert np.isfinite(baseline_model(time, "spline", [0.0, 0.1, 0.0])).all()

    parts = split_theta(
        [0.0, 0.1, 0.2, 0.01, 0.3, 0.02, 1.0, 0.0, 5000.0, 0.1],
        "sine",
        use_ratio_model=True,
    )
    np.testing.assert_array_equal(parts["baseline_params"], [0.3, 0.02, 1.0, 0.0])
    assert parts["teff"] == 5000.0
    assert parts["ebv"] == 0.1


def test_shared_flare_model_and_likelihood_reject_invalid_scatter():
    time = np.linspace(-0.05, 0.05, 51)
    theta = np.array([0.0, 0.02, 0.1, 0.003])
    model = model_flux(theta, time)

    assert model.shape == time.shape
    assert np.isfinite(model).all()
    assert np.isfinite(log_likelihood(theta, time, model, np.full(time.size, 0.005)))
    assert log_likelihood(theta * np.array([1.0, 1.0, 1.0, 0.0]), time, model, np.full(time.size, 0.005)) == -np.inf


def test_ratio_grid_loader_accepts_legacy_and_matrix_csv(tmp_path):
    legacy_path = tmp_path / "legacy.csv"
    legacy_path.write_text(
        "Teff_K,EBV_0.0,EBV_0.1\n5000,1.0,2.0\n6000,3.0,4.0\n",
        encoding="utf-8",
    )
    matrix_path = tmp_path / "matrix.csv"
    matrix_path.write_text(
        ",0.0,0.1\n5000,1.0,2.0\n6000,3.0,4.0\n",
        encoding="utf-8",
    )

    legacy = load_ratio_grid(legacy_path)
    matrix = load_ratio_grid(matrix_path)

    np.testing.assert_array_equal(legacy["ratio_values"], matrix["ratio_values"])
    assert legacy["interpolator"]((5500.0, 0.05)) == 2.5
    assert matrix["interpolator"]((5500.0, 0.05)) == 2.5
    assert infer_ratio_grid_bounds(legacy) == {
        "teff_min": 5000.0, "teff_max": 6000.0, "ebv_min": 0.0, "ebv_max": 0.1
    }


def test_ratio_grid_filter_resolution():
    assert normalize_ratio_filter(None) == "F087"
    assert normalize_ratio_filter(" f213 ") == "F213"
    assert resolve_ratio_grid_path("F087", "custom.csv").as_posix() == "custom.csv"


def test_shared_ratio_sampler_helpers(tmp_path):
    relative_path = tmp_path / "relative.csv"
    relative_path.write_text("time,value,err\n0.5,0.2,0.03\n", encoding="utf-8")
    relative = load_ratio_data_csv(relative_path)
    assert relative == {"time": 0.5, "value": 0.2, "err": 0.03, "input_mode": "relative"}

    absolute_path = tmp_path / "absolute.csv"
    absolute_path.write_text("time,flux_abs,err_abs\n0.5,1.2,0.03\n", encoding="utf-8")
    absolute = load_ratio_data_csv(absolute_path, aux_quiescent_flux=1.0, flare_quiescent_flux=2.0)
    assert np.isclose(absolute["value"], 0.2)
    assert absolute["input_mode"] == "absolute"

    ratio_grid = {
        "teff_values": np.array([5000.0, 6000.0]),  # [K]
        "ebv_values": np.array([0.0, 0.2]),  # [mag]
    }
    specs = ensure_ratio_param_specs([], ratio_grid)
    assert [spec["name"] for spec in specs] == ["teff", "ebv"]
    partial_specs = [{"name": "Teff", "prior_type": "fixed", "fixed_value": 5500.0}]
    completed_specs = ensure_ratio_param_specs(partial_specs, ratio_grid)
    assert [spec["name"].lower() for spec in completed_specs] == ["teff", "ebv"]
    assert len(partial_specs) == 1
    assert parameter_names("spline", 2, use_ratio_model=True) == [
        "tpeak", "fwhm", "ampl", "sigma", "spline_c0", "spline_c1", "teff", "ebv"
    ]
    truths = truth_vector_from_names(["spline_c1", "ebv"], {"spline_coeffs": [0.1, 0.2]})
    assert truths[0] == 0.2
    assert np.isnan(truths[1])