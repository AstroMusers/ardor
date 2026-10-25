from ardor.Data_Query.Catalog import Catalog


def test_catalog_listing_respects_tdpy_verbosity(tmp_path, capsys, monkeypatch):
    (tmp_path / "planets.csv").write_text("name\nExample b\n", encoding="utf-8")
    catalog = Catalog(tmp_path)

    monkeypatch.setenv("TDPY_VERBOSITY", "0")
    quiet_catalogs = catalog.list_available_catalogs()
    assert quiet_catalogs == ["planets.csv"]
    assert capsys.readouterr().out == ""

    monkeypatch.setenv("TDPY_VERBOSITY", "1")
    verbose_catalogs = catalog.list_available_catalogs()
    assert verbose_catalogs == ["planets.csv"]
    assert capsys.readouterr().out == "\U0001f4cb Available catalogs: 1\n   - planets.csv\n"
