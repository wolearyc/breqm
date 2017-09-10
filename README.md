README
======

BREQM - Boundary Region Embedding Quantum Mechanics/Molecular Mechanics

This software performs BREQM using VASP (for QM) and LAMMPS (for MM). The scripts contained within the main directory drive the calculation. Files used and produced by VASP and LAMMPS are contained within the vasp_run and lammps_run directories; these are called running directories, and must be set-up manually (continue reading!)

Requirements
------------
The code must be run with a python distribution with the numpy and ase packages installed. On WAGLAND, such a distribution can be found that /ul/wolearyc/.local/bin/python

Downloading
-----------
To download this code, simply type
```
git clone https://github.com/wolearyc/breqm
```
This will create a directory "breqm" containing this software.

Interfacing with VASP and LAMMPS
--------------------------------
First, create two directories, vasp_run and lammps_run, in the breqm directory. I call these "running directories", since it is where the VASP and LAMMPS calculations are run. The running directories must contain two scripts, lammps.run and vasp.run, which contain the commands to running the respective software. 

By default, this tool expects VASP to be run through PBS. Additionally, this tool expects LAMMPS to be run locally.

For LAMMPS, this tool only writes a data file ('lammps.data') and an input file ('lammps.in') within lammps_run. Therefore, all other files necessary for LAMMPS (e.g. force fields, control, etc.) must be  provided by the user within the lammps_run directory. The user can control the LAMMPS input file by changing a string in lammps.py.

For VASP, this tool only writes a POSCAR and an INCAR within vasp_run. Again, all other files required by VASP (KPOINTS, POTCAR, etc.) must be provided by the user within the vasp_run directory. The INCAR can be modified by changing a string in vasp.py. 


Running 
--------
An example run is included with the code. See example.py. It is also a good idea to read the SURF report (pdf) included with the code.

Overall, the total system should be initialized from a VASP file into an "Atoms" object. Then the BREQM regions should be established by instantiating a "Regions" object.  Finally, NVT or NVE dynamics can be run. 

Bonds should not be broken, so the user should specify a bond cutoff dictionary for relevant bonds. The cutoff dictionary (which can be empty) is passed when running dynamics. 

Proper order of atoms is essential to many parts of these scripts. The atom and element order specified in the input file is conserved. Because of this, the element order in the input VASP file should be identical to the desired atom order to be used in the LAMMPS calculation.


Output
------
This tool writes the trajectory to trj.xyz, which can be visualized in VMD. Remaining output is done to stdout. To redirect to a log file, try:
```
python -u example.py > log &
```

