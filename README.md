# README

This repository contains a python package for VIRUS2 data reduction

## Provenance

- The VIRUS2 data reduction pipeline was derived from the VIRUS-P data reduction pipeline.
- This present repository is a fork from an existing repo containing data reduction scripts for VIRUS-P and VIRUS2
- Original reduction pipeline `Antigen/reduct_virus2.py` was authored by [Greg Zeiman <grzeimann@gmail.com>](https://github.com/grzeimann) and [Maya Debski <maya.h.debski@gmail.com>](https://github.com/maya-debski)
- Original source code is available from the public GitHub repository here: `git clone https://github.com/maya-debski/Antigen`
- Overhaul of that script into an installable python package was developed by [Jason Vestuto <vestuto@gmail.com>](https://github.com/vestuto)

## Installation

- see [./INSTALL.md](./INSTALL.md)

## Usage

- Example use: [./Antigen/reduce_virus2.py](./docs/Antigen/reduce_virus2.py)
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
