# README

This repository contains a python package for VIRUS2 data reduction

## Provenance

- The VIRUS2 data reduction pipeline was derived from the VIRUS-P data reduction pipeline.
- This present repository is a "manual" fork from an [existing repo](https://github.com/maya-debski/Antigen) containing data reduction scripts for VIRUS-P and VIRUS2
- Original reduction pipeline `Antigen/reduct_virus2.py` was authored by [Maya Debski <maya.h.debski@gmail.com>](https://github.com/maya-debski) and [Greg Zeiman <grzeimann@gmail.com>](https://github.com/grzeimann)
- Original source code is available from the public GitHub repository [Antigen](https://github.com/maya-debski/Antigen): `git clone https://github.com/maya-debski/Antigen`
- Overhaul of that script into an installable python package was developed by [Jason Vestuto <vestuto@gmail.com>](https://github.com/vestuto)
- Other major changes to local/global scoping was made to help with anticipated updates and long term care.
- Didn't want to force the original devs into a PR, so all restructuring here was done in a stand-alone repo, in smaller commits to help illustrate changes, but still allow a rebase if wanted later after travel allowed a group discussion. 


## Installation

- see [./INSTALL.md](./INSTALL.md)

## Usage: option 1: Restructured package 

- `conda activate env_virus2_reduction`
- `v2reduce_process.py --help`

## Usage: option 2: original repo

- See [Antigen](https://github.com/maya-debski/Antigen) for help. Usage example from email communication below:
- Original repo does not structure the code as a python package, so there is no install; just copy/download and run using the explicit file path as a stand-alone script.
- Example use: [./Antigen/reduce_virus2.py](./Antigen/reduce_virus2.py)
```bash
(base) cns-r-pmaa65420:~ grz85$ python /Users/grz85/work/Antigen/reduce_virus2.py /Users/grz85/work/v2_data/ /Users/grz85/work/v2_data/reduc -ra
[INFO - 2025-06-13 14:25:26,064] Sorting Files
[INFO - 2025-06-13 14:25:26,067] Making master bias frames
[INFO - 2025-06-13 14:25:27,404] Making master flat frames
[INFO - 2025-06-13 14:25:28,002] Making master arc frames
[INFO - 2025-06-13 14:25:28,600] Getting trace for each master flat
[INFO - 2025-06-13 14:25:29,119] Getting wavelength for each master arc
[INFO - 2025-06-13 14:25:29,876] Getting fiber to fiber for each master domeFlat
```
- It will produce fiber extracted spectra in the `$CURDIR/reduc` folder for each (FITs file exposure) frame of type `science`.  
- To view the test image, use the `$CURDIR/IFUcen/IFUcen_VIRUS2_D3G.txt` file as follows:
