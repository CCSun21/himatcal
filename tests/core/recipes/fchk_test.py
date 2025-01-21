from __future__ import annotations

from himatcal import SETTINGS
from himatcal.recipes.electrolyte.sol_stru._base import formchk

formchkpath = str(SETTINGS.FORMCHK_PATH)


def test_formchk(mocker):
    # Mock the subprocess.run method
    mock_run = mocker.patch("subprocess.run")

    # Define the test input
    chk_file = "/home/suncc/Code/pub/himatcal/examples/EleCal/GAFF/DJN-EC-gas.chk"

    # Call the function
    formchk(chk_file)

    # Assert that subprocess.run was called with the correct arguments
    mock_run.assert_called_once_with(
        f"{formchkpath} {chk_file}", shell=True, check=True
    )
