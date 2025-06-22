r"""Dependency checking

Adapted from GLUE: https://github.com/gao-lab/GLUE
"""

import importlib
import types
from abc import abstractmethod
from importlib.metadata import version

from packaging.version import parse


class Checker:
    r"""
    Checks availability and version of a dependency

    Parameters
    ----------
    name
        Name of the dependency
    vmin
        Minimal required version
    install_hint
        Install hint message to be printed if dependency is unavailable
    """

    def __init__(self, name: str, vmin: str | None = None, install_hint: str | None = None) -> None:
        self.name = name
        self.vmin = parse(vmin) if vmin else vmin
        vreq = f" (>={self.vmin})" if self.vmin else ""
        self.vreq_hint = f"This function relies on {self.name}{vreq}."
        self.install_hint = install_hint

    @abstractmethod
    def check(self) -> None:
        r"""Check availability and version"""
        raise NotImplementedError  # pragma: no cover


class ModuleChecker(Checker):
    r"""
    Checks availability and version of a Python module dependency

    Parameters
    ----------
    name
        Name of the dependency
    vmin
        Minimal required version
    install_hint
        Install hint message to be printed if dependency is unavailable
    """

    def __init__(
        self, name: str, package_name: str | None = None, vmin: str | None = None, install_hint: str | None = None
    ) -> None:
        super().__init__(name, vmin, install_hint)
        self.package_name = package_name or name

    def check(self) -> None:  # noqa: D102
        try:
            importlib.import_module(self.name)
        except ModuleNotFoundError as e:
            raise RuntimeError(" ".join(filter(None, [self.vreq_hint, self.install_hint]))) from e
        v = parse(version(self.package_name))
        if self.vmin and v < self.vmin:
            raise RuntimeError(
                " ".join(
                    [
                        self.vreq_hint,
                        f"Detected version is {v}.",
                        "Please install a newer version.",
                        self.install_hint or "",
                    ]
                )
            )


INSTALL_HINTS = types.SimpleNamespace(
    scvi="Please install scvi-tools, either directly or via the corresponding hnoca extra: `pip install 'hnoca[mapping]'`.",
    scarches="Please install scarches, either directly or via the corresponding hnoca extra: `pip install 'hnoca[mapping]'`.",
    decoupler="Please install decoupler, either directly or via the corresponding hnoca extra: `pip install 'hnoca[stats]'`.",
    cuml="Please install cuML from rapids: https://docs.rapids.ai/install/. ",
)


CHECKERS = {
    "scvi-tools": ModuleChecker("scvi", vmin="1.2", install_hint=INSTALL_HINTS.scvi),
    "scarches": ModuleChecker("scarches", vmin="0.6.1", install_hint=INSTALL_HINTS.scarches),
    "decoupler": ModuleChecker("decoupler", vmin="1.6", install_hint=INSTALL_HINTS.decoupler),
    "cuml": ModuleChecker("cuml", vmin=None, install_hint=INSTALL_HINTS.cuml),
}


def check_deps(*args) -> None:
    r"""
    Check whether certain dependencies are installed

    Parameters
    ----------
    args
        A list of dependencies to check
    """
    for item in args:
        CHECKERS[item].check()
