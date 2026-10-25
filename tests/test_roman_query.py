from astropy.table import Table

from ardor.Roman import query


def test_query_gaia_field_respects_tdpy_verbosity(tmp_path, monkeypatch, capsys):
    class Job:
        def get_results(self):
            return Table({"source_id": [123]})

    monkeypatch.setattr(query.Gaia, "launch_job_async", lambda statement: Job())
    output_path = tmp_path / "gaia_results.csv"

    monkeypatch.setenv("TDPY_VERBOSITY", "0")
    quiet_result = query.query_gaia_field(output_path, N=4)
    assert capsys.readouterr().out == ""
    assert quiet_result["source_id"].tolist() == [123]
    assert output_path.is_file()

    monkeypatch.setenv("TDPY_VERBOSITY", "1")
    query.query_gaia_field(output_path, N=4)
    assert capsys.readouterr().out == f"Writing to {output_path}...\n"
