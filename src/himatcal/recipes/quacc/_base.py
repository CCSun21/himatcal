from __future__ import annotations

from pathlib import Path

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


def clear_quacc_cache(path: Path | None = None) -> None:
    if path is None:
        path = Path.cwd()
    """Clear the QUACC cache directory.

    This function clears the QUACC cache directory by removing all files and directories within it.
    The cache directory is specified by the provided path.

    Args:
        path (Path, optional): The path to the QUACC cache directory. Defaults to the current working directory.

    Examples:
        clear_quacc_cache(Path("path/to/cache_dir"))
    """
    directory_patterns = ["quacc-*", "tmp*"]

    for pattern in directory_patterns:
        for dir_path in path.glob(pattern):
            if dir_path.is_dir():
                for file in dir_path.iterdir():
                    if file.is_file():
                        file.unlink()
                dir_path.rmdir()
            else:
                dir_path.unlink()
