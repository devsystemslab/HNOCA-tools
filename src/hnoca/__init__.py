from importlib.metadata import version

from . import snapseed, stats

__all__ = ["snapseed", "stats"]

__version__ = version("hnoca")
