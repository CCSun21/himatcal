# himatcal

![PyPI - version](https://img.shields.io/pypi/v/himatcal)
![supported python versions](https://img.shields.io/pypi/pyversions/himatcal)
![PyPI - Downloads](https://img.shields.io/pypi/dd/himatcal)

Some scripts to perform material simulation.

> [!WARNING]
> 🚧 This repository is still under construction. 🚧

## Installation

This project uses [uv](https://github.com/astral-sh/uv) for dependency management. To install the package and its dependencies:

1. First, ensure you have uv installed:

   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. Clone the repository and install dependencies:

   ```bash
   git clone https://github.com/CCSun21/himatcal.git
   cd himatcal
   uv pip install -e .
   ```

3. To create and activate a virtual environment:
   ```bash
   uv venv
   source .venv/bin/activate  # On Unix/macOS
   # or
   .venv\Scripts\activate  # On Windows
   ```

Optional dependencies are organized into extras. To install specific extras:

- For molecule-related features: `uv pip install -e ".[molecule]"`
- For development tools: `uv pip install -e ".[dev]"`
