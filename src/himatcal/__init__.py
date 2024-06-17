"""Init data for himatcal package."""

from __future__ import annotations

from importlib.metadata import version

from himatcal.settings import HimatcalSettings

# load the version from the pyproject.toml file
__version__ = version("himatcal")

# load the settings
SETTINGS = HimatcalSettings()
