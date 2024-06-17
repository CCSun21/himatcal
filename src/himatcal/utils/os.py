"""os utilities for himatcal"""

import os
import re

from monty.os import makedirs_p


def labeled_dir(main_workdir, label):
    """
    Create a new folder in the main_workdir with the label.

    Args:

        main_workdir (str): The main working directory.
        label (str): The label of the folder.

    Returns:

        folder_path (str): The path of the new folder.
    """
    # Get the folder names in main_workdir
    folder_names = [
        name
        for name in os.listdir(main_workdir)
        if os.path.isdir(os.path.join(main_workdir, name))
    ]
    numbers = [
        int(re.search(r"\d+", name).group())
        for name in folder_names
        if re.search(r"\d+", name)
    ]
    new_number = max(numbers) + 1 if numbers else 1
    # Create new folder
    folder_name = f"{new_number:02d}.{label}"
    folder_path = os.path.join(main_workdir, folder_name)
    makedirs_p(folder_path)
    print(f"Created new folder: {folder_path}")
    return folder_path


def get_chg_mult(molname):
    """Get the label, chg and mult from the name of a molecule, format: {label}-c{charge}s{mult}"""
    pattern = r"(.*?)-c(n?\d)s(\d+)"
    match = re.match(pattern, molname)
    if match:
        label, chg, mult = match.groups()
        chg = "-" + chg[1:] if chg.startswith("n") else chg
        return label, int(chg), int(mult)
    else:
        return None, None, None
    
def write_chg_mult_label(label, chg, mult):
    """Write the label, chg and mult to a string, format: {label}-c{charge}s{mult}"""
    if chg < 0:
        chg = f"n{abs(chg)}"
    return f"{label}-c{chg}s{mult}"
