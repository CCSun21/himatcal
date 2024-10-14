# Installation

## using micromamba

!!! Note
    Installing with micromamba and uv is the recommended way to install himatcal.

```sh
micromamba env create -n himatcal python==3.10
micromamba activate himatcal

pip install uv
uv pip install himatcal
uv pip install git+https://gitlab.com/ase/ase.git
```
