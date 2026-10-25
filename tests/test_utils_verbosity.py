import pandas as pd

from ardor.Utils.Utils import transfer_columns_between_files


def test_transfer_columns_respects_tdpy_verbosity(tmp_path, capsys, monkeypatch):
    source_path = tmp_path / "source.csv"
    destination_path = tmp_path / "destination.csv"
    output_path = tmp_path / "output.csv"
    pd.DataFrame({"target": ["A", "B"], "radius": [1.1, 2.2]}).to_csv(source_path, index=False)
    pd.DataFrame({"target": ["A", "B"]}).to_csv(destination_path, index=False)

    monkeypatch.setenv("TDPY_VERBOSITY", "0")
    quiet_result = transfer_columns_between_files(
        source_path, destination_path, "target", ["radius"], output_path
    )
    assert capsys.readouterr().out == ""
    pd.testing.assert_series_equal(quiet_result["radius"], pd.Series([1.1, 2.2], name="radius"))

    monkeypatch.setenv("TDPY_VERBOSITY", "1")
    transfer_columns_between_files(
        source_path, destination_path, "target", ["radius"], output_path
    )
    assert capsys.readouterr().out == f"\u2705 Data transferred. Saved to {output_path}\n"
