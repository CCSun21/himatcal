"""Settings for himatcal"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pydantic import Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

if TYPE_CHECKING:
    from typing import Any

_DEFAULT_CONFIG_FILE_PATH = Path("~", ".himatcal.yaml").expanduser().resolve()


class HimatcalSettings(BaseSettings):
    """
    Settings for himatcal
    """

    model_config = SettingsConfigDict(
        env_prefix="himatcal_",
        env_nested_delimiter="__",
        env_parse_none_str="None",
        extra="forbid",
        validate_assignment=True,
    )

    CONFIG_FILE: Path = Field(
        _DEFAULT_CONFIG_FILE_PATH, description=("Path to the YAML config file")
    )

    # ---------------------------
    # external software Settings
    # ---------------------------
    XTB_EXE_PATH: Path = Field(Path(), description=("Path to the xtb executable"))

    CREST_EXE_PATH: Path = Field(Path(), description=("Path to the crest executable"))

    CREST_EXE_PATH_V3: Path = Field(
        Path(), description=("Path to the crest executable")
    )

    ISOSTAT_EXE_PATH: Path = Field(
        Path(), description=("Path to the isostat executable")
    )

    FORMCHK_PATH: Path = Field(Path(), description=("Path to the formchk executable"))

    MULTIWFN_PATH: Path = Field(Path(), description=("Path to the Multiwfn executable"))

    SOBTOP_PATH: Path = Field(Path(), description=("Path to the sobtop executable"))

    OBABEL_PATH: Path = Field(Path(), description=("Path to the obabel executable"))

    GROMACS_PATH: Path = Field(Path(), description=("Path to the gromacs executable"))

    ORCA_PATH: Path = Field(Path(), description=("Path to the ORCA executable"))

    # ---------------------------
    # external modle Settings
    # ---------------------------
    MABASIS_PATH: Path = Field(Path(), description=("Path to the mabasis"))

    AIMNET2_MODEL_PATH: Path = Field(Path(), description=("Path to the aimnet2 model"))

    # ---------------------------
    # external api Settings
    # ---------------------------

    MAPI_KEY: str = Field("", description=("API key for the Materials Project"))

    CHEMSPIDER_API_KEY: str = Field("", description=("API key for ChemSpider"))

    # ---------------------------
    # external key file
    # ---------------------------
    SGHPC1_KEY_PATH: Path = Field("", description=("Path to the ssh key file of sghpc1"))
    YEESUAN_KEY_PATH: Path = Field("", description=("Path to the ssh key file of yeesuan"))
    XMU_KEY_PATH: Path = Field("", description=("Path to the ssh key file of XMU"))

    # ---------------------------
    # MongoDB URI
    # ---------------------------
    MONGODB_URI: str = Field(
        "", description=("MongoDB URI for the database connection")
    )

    @model_validator(mode="before")
    @classmethod
    def load_user_setting(cls, settings: dict[str, Any]) -> dict[str, Any]:
        """Load user settings"""
        return _type_handler(_use_custom_config_settings(settings))

    # ---------------------------
    # logging
    # ---------------------------
    DEBUG: bool = Field(
        False, description=("Whether to enable debug logging"))

def _use_custom_config_settings(settings: dict[str, Any]) -> dict[str, Any]:
    """Use custom settings from the config file"""
    from monty.serialization import loadfn

    config_file_path = (
        Path(settings.get("CONFIG_FILE", _DEFAULT_CONFIG_FILE_PATH))
        .expanduser()
        .resolve()
    )

    new_settings = {}  # type: dict
    if config_file_path.exists() and config_file_path.stat().st_size > 0:
        new_settings |= loadfn(config_file_path)

    new_settings.update(settings)
    return new_settings


def _type_handler(settings: dict[str, Any]) -> dict[str, Any]:
    """
    Convert common strings to their proper types.

    Parameters
    ----------
    settings : dict
        Initial settings.

    Returns
    -------
    dict
        Updated settings.
    """
    for key, value in settings.items():
        if isinstance(value, str):
            if value.lower() in {"null", "none"}:
                settings[key] = None
            elif value.lower() in {"true", "false"}:
                settings[key] = value.lower() == "true"

    return settings
