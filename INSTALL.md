# INSTALL

Instructions on how to install the `v2reduce` package

## Summary: Manual Install

- `conda env create -f environment.yml`
- `conda activate env_virus2_reduction`
- `pip install . --no-deps`

## Summary: Automated Install

- Install conda: [./install_conda.sh](install_conda.sh)
- Install `v2reduce` package: [./install_v2reduce.sh](install_v2reduce.sh)

## Details: Conda Environment Overview

- `conda` is required for this package installation. See [`./install_conda.sh`](`./install_conda.sh`)
- Run-time virtual environment is created and managed with `conda` under the users `$HOME` path.
- All external package dependencies and compatible versions are managed by `conda`
- The reduction packaging is as a `pip` package, and installed into the `conda` env with `pip install`

## Details: Conda Environment Creation

- Download conda installer from https://github.com/conda-forge/miniforgeInstall 
- Install conda under `$HOME/opt/conda` for user not system
- Create conda environment: `conda env create -f environment.yml`
- Activate conda env: `conda activate env_virus2_reduction`

## Details: Pip Install into Environment

- Confirm that `python` and `pip` live in the conda env
    - `echo $CONDA_PREFIX`
    - `which python pip`
- Install into conda env: 
    - `cd <repo-root-path>`
    - `pip install . --no-deps`
- Optional install methods:
    - editable install for dev testing: `pip install --editable . --no-deps`
    - optional: for testing dependencies after install: `pip install --no-build-isolation --no-deps .`

## Details: Pip Packaging, for Development and Deployment

Building Packages:
- activate the conda env: `conda activate env_virus2_reduction`
- install the pip package `build` tool: `conda install -c conda-forge build`
- build a src dist: 
    - build: `python -m build --sdist`
    - verify: `ls -l ./dist/v2reduce-2026.6.16.tar.gz`
    - install: `pip install ./dist/v2reduce-2026.6.16.tar.gz --no-deps`
- build a wheel: 
    - build: `python -m build --wheel`
    - verify: `ls -l ./dist/v2reduce-2026.6.16-py3-none-any.whl`
    - install: `pip install ./dist/v2reduce-2026.6.16-py3-none-any.whl --no-deps`
