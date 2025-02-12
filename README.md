# himatcal

![PyPI - version](https://img.shields.io/pypi/v/himatcal)
![supported python versions](https://img.shields.io/pypi/pyversions/himatcal)
![PyPI - Downloads](https://img.shields.io/pypi/dd/himatcal)

Some scripts to perform material simulation.

> [!WARNING]
> 🚧 This repository is still under construction. 🚧

## Installation

This project uses [Poetry](https://python-poetry.org/) for dependency management. To install the package and its dependencies:

1. First, ensure you have Poetry installed:
   ```bash
   curl -sSL https://install.python-poetry.org | python3 -
   ```

2. Clone the repository and install dependencies:
   ```bash
   git clone https://github.com/CCSun21/himatcal.git
   cd himatcal
   poetry install
   ```

3. To activate the virtual environment:
   ```bash
   poetry shell
   ```

Optional dependencies are organized into groups. To install specific groups:
- For molecule-related features: `poetry install --with molecule`
- For development tools: `poetry install --with dev`

