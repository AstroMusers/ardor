import importlib.util
import os
from pathlib import Path

PATHS_FILE = Path(__file__).resolve().parents[1] / "src" / "ardor" / "paths.py"


def load_paths():
    specification = importlib.util.spec_from_file_location("ardor_paths", PATHS_FILE)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_repository_path_defaults_to_checkout(monkeypatch):
    monkeypatch.delenv("ARDOR_PATH", raising=False)
    module = load_paths()

    expected_path = Path(module.__file__).resolve().parents[2]
    assert os.environ["ARDOR_PATH"] == str(expected_path)
    assert module.get_repository_path() == expected_path
    assert module.get_data_path() == expected_path / "data"
    assert module.get_visuals_path() == expected_path / "visuals"


def test_repository_path_honors_environment_override(tmp_path, monkeypatch):
    monkeypatch.setenv("ARDOR_PATH", str(tmp_path))
    paths = load_paths()

    assert paths.get_repository_path() == tmp_path.resolve()