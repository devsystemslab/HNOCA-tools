from importlib.resources import files

import pytest
import yaml


@pytest.fixture
def marker_h1():
    data_path = files("tests.data") / "snapseed_markers_h1.yaml"
    return yaml.safe_load(data_path.read_text())


@pytest.fixture
def marker_h2():
    data_path = files("tests.data") / "snapseed_markers_h2.yaml"
    return yaml.safe_load(data_path.read_text())


@pytest.fixture
def marker_h3():
    data_path = files("tests.data") / "snapseed_markers_h3.yaml"
    return yaml.safe_load(data_path.read_text())
