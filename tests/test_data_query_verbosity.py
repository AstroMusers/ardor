from ardor.Data_Query.Data_Query import Bulk_TESS_lc_Query


def test_bulk_tess_query_respects_tdpy_verbosity(tmp_path, capsys, monkeypatch):
    download_directory = tmp_path / "downloads"
    download_directory.mkdir()
    (download_directory / "target-a").mkdir()
    arguments = ([1.0], [2.0], [123], download_directory, ["target-a"])

    monkeypatch.setenv("TDPY_VERBOSITY", "0")
    Bulk_TESS_lc_Query(*arguments)
    assert capsys.readouterr().out == ""

    (download_directory / "target-a").mkdir()
    monkeypatch.setenv("TDPY_VERBOSITY", "1")
    Bulk_TESS_lc_Query(*arguments)
    output = capsys.readouterr().out
    assert "The directory already exists!" in output
    assert "The already existing directories are:" in output
    assert "The potential undownloaded files are: []" in output
