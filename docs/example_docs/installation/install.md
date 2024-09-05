# Pip Installing

Now it's time to install your Python package! You will want to install your Python package in "editable" mode, which means you won't have to re-install your code every time you make updates to it. Additionally, you will want to install several optional dependencies (listed under the `[project.optional-dependencies]` header in `pyproject.toml`) to ensure that you can test and build the documentation for your code.

With all this in mind, you will want to run the following in the command line from the base of the package directory:

```bash
pip install -e .[dev,docs]
```

Here, the `-e` means editable mode, the `.` means the current directory, and the `[dev,docs]` means it will install the "dev" and "docs" optional dependency set listed in the `pyproject.toml` file.

!!! Tip

    You should generally start from a clean Python environment, such as a new Conda environment if you are using Anaconda or one of its variants.

To make sure you installed your package successfully, open a Python console and run `import <MyPackageName>`. It should return without any errors. If there are errors, it's likely because you forgot to replace a "template" placeholder with the name of your package.

# Prepare the Environment

## using micromamba

```sh
micromamba env create -n himatcal python==3.10
micromamba activate himatcal

pip install uv
uv pip install himatcal
uv pip install git+https://gitlab.com/ase/ase.git
```

## extras

```sh
micromamba install -c conda-forge xtb
uv pip install git+https://github.com/Quantum-Accelerators/xtb_ase.git
uv pip install git+https://github.com/CCSun21/pysisyphus.git@dev
```

## aimnet2

```sh
pip install torch-cluster
```