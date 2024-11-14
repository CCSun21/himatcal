"""himatcal package."""

from __future__ import annotations

import logging
from importlib.metadata import version

from himatcal.settings import HimatcalSettings

# load the version from the pyproject.toml file
__version__ = version("himatcal")

# load the settings
SETTINGS = HimatcalSettings()

# set logging info
logging.basicConfig(
    level=logging.DEBUG if SETTINGS.DEBUG else logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
