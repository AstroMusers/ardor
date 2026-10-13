import numpy as np

from ardor.Roman.flare_model import (
    baseline_model,
    infer_ratio_grid_bounds,
    load_ratio_grid,
    log_likelihood,
    model_flux,
    normalize_ratio_filter,
    resolve_ratio_grid_path,
    split_theta,
)


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