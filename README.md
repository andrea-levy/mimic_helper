# mimic_helper
Collection of python scripts to make working with MiMiC easier

## Run the scripts
All the scripts are coded in python and require standard packages.
 
Anyway, the required packages can be found in the `environment.yml` file and installed in a conda environment
```
conda env create -f environment.yml
conda activate mimic_helper
```

All the scripts are contained in the `script` folder and present a help fuction, providing a brief description and the input needed.
This can be accessed via
```
python3 script_name.pi -h
```
The scripts accept positional and optional arguments, to be specified via the command line. 

## Scripts
The scripts in `mimic_helper` are designed to automatize repetitive processes when working with MiMiC.
In particular, they are used for pre- or post-processing of GROMACS and/or CPMD files.

The scripts contained in `mimic_helper` include:
* CPDM scripts
	* `temp_check.py`	- Extract temperature for QM and MM regions along a trajectory
	* `geofile_extract.py`	- Extract the geometry from a trajectory at a selected step
	* `traj_xyz_convert.py`	- Convert a trajectory file into a .xyz file, easily visualizable with e.g. VMD
