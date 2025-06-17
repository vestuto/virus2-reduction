#!/bin/bash
#----------------------------------------
# Purpose: Install v2reduction into the conda environment
#----------------------------------------

# Is conda installed?
echo "Check if conda is installed and available on PATH ..."
if ! command -v conda >/dev/null 2>&1; then
    echo "Error: Conda is not installed or not on $PATH. Check your ~/.bashrc or install conda. Exiting ... " >&2
    exit 1
fi

# Deactivate env stack ...
echo "Notice: deactivate all envs in the activated env stack ..."
while [[ $CONDA_SHLVL -gt 0 ]]; do conda deactivate; done

# Create new conda env ...
echo "Notice: Creating new conda env from ./environment.yml ..."
conda env create -f ./environment.yml
conda clean

# Install v2reduction pkg into conda env ...
echo "Notice: Activating the conda env and pip installing the pkg into it ..."
conda activate env_virus2_reduction
pip install . --no-deps

# Finish ...
echo "Notice: v2reduction pkg has been installed into env_virus2_reduction"
echo "Notice: to use v2reduction, in every new shell run, $ conda activate env_virus2_reduction"
