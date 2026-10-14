import importlib
import sys

import numpy as np
import pytest
from astropy.io import fits

from ardor.Flares import pcat_flares


def test_flare_bin_edges_reuse_upstream_midpoints_and_preserve_positive_axis(monkeypatch):
    from unittest.mock import Mock
    import tdpy.numerics

    shared_edges = Mock(wraps=tdpy.numerics.midpoint_bin_edges)
    monkeypatch.setattr(tdpy.numerics, "midpoint_bin_edges", shared_edges)
    time_days = np.array([-2.0, -1.0, 3.0, 5.0, 8.0])  # [day]
    edges, offset = pcat_flares._time_bin_edges(time_days)
    original_edges = np.array([-2.5, -1.5, 1.0, 4.0, 6.5, 9.5])  # [day]
    assert offset == 3.5
    np.testing.assert_array_equal(edges, original_edges + offset)
    assert edges[0] == 1.0
    np.testing.assert_array_equal(np.diff(edges), np.diff(original_edges))
    shared_edges.assert_called_once()
    np.testing.assert_array_equal(shared_edges.call_args.args[0], time_days)


def test_allesfitter_is_not_required_to_import_pcat_backend(monkeypatch):
    monkeypatch.setitem(sys.modules, "allesfitter", None)

    module = importlib.import_module("ardor.Flares.allesfitter_priors")

    assert callable(module.model_compare)


def test_variable_flare_fit_uses_pcat_transdimensional_model(tmp_path, monkeypatch):
    from pcat import sampling
    from pcat import main as pcat_main

    captured = {}
    posterior = object()
    run_directory = tmp_path / "Flare_1"

    def fake_sample(**configuration):
        captured.update(configuration)
        return object()

    monkeypatch.setattr(sampling, "sample", fake_sample)
    monkeypatch.setattr(pcat_main, "readfile", lambda path: posterior)
    time_days = np.linspace(-0.04, 0.04, 41)
    flux = 1.0 + 0.025 * np.exp(-0.5 * (time_days / 0.006) ** 2)
    uncertainty = np.full(time_days.size, 0.003)
    fit = pcat_flares.fit_variable_flares_pcat(
        time_days,
        flux,
        uncertainty,
        run_directory,
        reference_bjd=2460000.0,
        numbswep=100,
        max_flares=7,
        numbproc=1,
    )

    assert fit.posterior is posterior
    assert fit.reference_bjd == 2460000.0
    assert captured["typeexpr"] == "fire"
    assert captured["spectype"] == ["flarfred"]
    assert captured["fittminmnumbelempop0"] == 0
    assert captured["fittmaxmnumbelempop0"] == 7
    assert captured["strgcnfg"] == "Flare_1"
    assert captured["binsenerfull"][0] > 0.0
    counts = fits.getdata(run_directory / "data" / "inpt" / "sbrt.fits")
    exposures = fits.getdata(run_directory / "data" / "inpt" / "expo.fits")
    assert counts.shape == (time_days.size, 1, 1, 1)
    assert exposures.shape == (time_days.size, 1, 1)
    assert np.isfinite(exposures).all() and np.all(exposures > 0.0)


def test_planck_law_and_integrator_yield_positive_flux():
    from ardor.Utils.planck_law import planck_law, planck_integrator

    density = planck_law(np.array([1e-5, 2e-5]), 5778.0)
    assert np.all(density > 0.0)
    assert np.all(np.isfinite(density))

    integral = planck_integrator(1e-5, 2e-5, 5778.0)
    assert integral > 0.0
    assert np.isfinite(integral)


def test_tier2_csv_reader_passes_reference_bjd_and_flux_data(tmp_path, monkeypatch):
    target_file = tmp_path / "Flare_3.csv"
    print(f"Writing to {target_file}...")
    np.savetxt(
        target_file,
        np.array([[0.0, 1.02, 0.01], [0.01, 1.01, 0.01], [0.02, 1.0, 0.01]]),
        delimiter=",",
        header="2460001.5",
    )
    captured = {}

    def fake_fit(time_days, flux, uncertainty, run_directory, **kwargs):
        captured.update(
            time_days=time_days,
            flux=flux,
            uncertainty=uncertainty,
            run_directory=run_directory,
            **kwargs,
        )
        return "fit"

    monkeypatch.setattr(pcat_flares, "fit_variable_flares_pcat", fake_fit)

    result = pcat_flares.fit_tier2_csv_pcat(target_file, tmp_path / "run")

    assert result == "fit"
    assert captured["reference_bjd"] == 2460001.5
    np.testing.assert_allclose(captured["time_days"], [0.0, 0.01, 0.02])
    np.testing.assert_allclose(captured["flux"], [1.02, 1.01, 1.0])
    np.testing.assert_allclose(captured["uncertainty"], 0.01)


def test_variable_flare_fit_rejects_nonpositive_flux_and_reuses_no_output(tmp_path):
    time_days = np.linspace(0.0, 0.1, 10)
    uncertainty = np.full(time_days.size, 0.01)
    with pytest.raises(ValueError, match="flux must be positive"):
        pcat_flares._flare_fit_configuration(
            time_days, np.r_[0.0, np.ones(time_days.size - 1)], uncertainty
        )

    run_directory = tmp_path / "occupied"
    run_directory.mkdir()
    (run_directory / "keep.txt").write_text("existing run", encoding="utf-8")
    with pytest.raises(FileExistsError, match="not empty"):
        pcat_flares.fit_variable_flares_pcat(
            time_days,
            np.ones(time_days.size),
            uncertainty,
            run_directory,
            numbswep=100,
            numbproc=1,
        )