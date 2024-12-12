from __future__ import annotations

from monty.serialization import loadfn


def load_quacc_result(file_path: str) -> dict:
    """Load and return the contents of a QUACC result file.

    This function reads a specified file and loads its contents into a dictionary format.
    It utilizes a loading function to handle the file reading process.

    Args:
        file_path (str): The path to the QUACC result file to be loaded.

    Returns:
        dict: The contents of the loaded QUACC result file.

    Examples:
        result = load_quacc_result("path/to/result_file.json")
    """
    return loadfn(file_path)
