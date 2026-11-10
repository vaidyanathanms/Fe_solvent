# Simulation and Analysis Codes for Iron-TFSI systems

Repository contains all the GROMACS simulation and Python/FORTRAN analysis codes for understanding the structural and dynamic properties of Iron (II) and Iron (III) systems with TFSI anions in water.

## Features

- Code to analyze scattering intensities from radial distribution functions using Octave open-source software. Please note that using the code for other cations or anions require editing the `compute_Iofq.m` in `src_mat` to reflect the scattering factors. 
- Can be extended easily to other anions by just adding the corresponding `top/itp` files in GROMACS.
- Analysis codes are robust to different anions and solvents.
- Creates GROMACS-ready simulation codes using Python (see `src_py/genconf.py`).

## Directory Structure

- `src_gmx`: Contains all the GROMACS files (itp/top and initial structures).
- `src_pack`: Contains PACKMOL codes for generating initial structures.
- `src_mat`: Contains Octave open-source codes for generating scattering intensities from radial distribution functions.
- `src_f90`: FORTRAN codes for analyzing outputs from GROMACS trajectories to compute various properties.
- `src_py`: Python codes for analyzing GROMACS trajectories and plotting results.
- `src_sh`: Codes for running simulations in HPC (specifically the HPCs at ORNL).

## Contact

Author: Vaidyanathan M. Sethuraman
* **Email:** [mathamangal1@ornl.gov](mailto:mathamangal1@ornl.gov)
* **Email:** [msvaidyanathan1729[at]gmail[dot]com](mailto:msvaidyanathan1729[at]gmail[dot]com)


## License

Distributed under the MIT License. See `LICENSE` for more information.
