import pytest

from hnoca.utils.check import CHECKERS, ModuleChecker


def test_import_name_differs_from_distribution():
    # scvi-tools is imported as `scvi`
    assert CHECKERS["scvi-tools"].name == "scvi"
    assert CHECKERS["scvi-tools"].package_name == "scvi-tools"


def test_no_version_lookup_without_vmin():
    # e.g. cuML is distributed as `cuml-cu12`; without vmin the distribution name must not matter
    ModuleChecker("json", package_name="not-a-distribution").check()


def test_missing_module():
    with pytest.raises(RuntimeError, match="relies on not_a_module"):
        ModuleChecker("not_a_module").check()
