from importlib.resources import files

import pytest
import yaml


@pytest.fixture
def marker_hierarchy():
    # Here, "tests.data" is the package corresponding to tests/data (thanks to __init__.py files)
    data_path = files("tests.data") / "snapseed_markers_hierachy.yaml"
    return yaml.safe_load(data_path.read_text())


@pytest.fixture
def marker_flat():
    data_path = files("tests.data") / "snapseed_markers_flat.yaml"
    return yaml.safe_load(data_path.read_text())
