README
======

BREQM - Boundary Region Embedding Quantum Mechanics/Molecular Mechanics

This software performs BREQM using VASP (for QM) and LAMMPS (for MM). The scripts contained within the main directory drive the calculation. Files used and produced by VASP and LAMMPS are contained within the vasp_run and lammps_run directories; these are called running directories.

Requirements
------------
The tool must be run with a python distribution with the numpy and ase packages installed. On WAGLAND, such a distribution can be found that /ul/wolearyc/.local/bin/python

Overview
--------
First, the total system should be initialized from a VASP file into an Atoms object. Then the BREQM regions should be established by instantiating a Regions class. Finally, NVT or NVE dynamics can be run. 

Atoms objects can contain bonding information, which is especially important whenever bonds should be preserved when splitting. The user should specify a bond cutoff dictionary for relevant bonds. The cutoff dictionary (which can be empty) is passed when running dynamics. 

Proper order of atoms is essential to many parts of these scripts. The atom and element order specified in the input file is conserved. Because of this, the element order in the input VASP file should be identical to the desired atom order to be used in the LAMMPS calculation.

Interfacing with VASP and LAMMPS
--------------------------------
Two scripts, vasp.run and lammps.run, should be present in the running directories. These scripts should run VASP and LAMMPS, respectively. Either software can be run through PBS, though I recommend only running VASP through PBS. If running with PBS, ensure that the block_pbs function is used before the atoms are updated; this code is contained in the calc_forces functions in lammps.py and vasp.py.

For LAMMPS, this tool only writes a data file ('lammps.data') and an input file ('lammps.in'). The data file comes from an Atoms object, but the input file is written by the user in string form (see lammps.py). This string should contain two format specifiers at two specific locations. The first specifier should be located when the atom order is specified in the pair_coeff command.
```
pair_coeff      * * lammps.ffield {0}
```
The second specifier should be located on some group that will be frozen later.
```
group           frozen id {1}
```

For VASP, this tool only writes a POSCAR and an INCAR. The POSCAR comes form an Atoms object, but the input file is written by the user in the string form (see vasp.py). This string should not contain any format specifiers, since all atom fixing is done using the 'Selective Dynamics' tag within the POSCAR.

Besides the files explicitly discussed above, the user is responsible for adding the remaining files required to perform the LAMMPS and VASP calculations within the running directories.

Output
------
This tool writes the trajectory to trj.xyz, which can be visualized in VMD. Remaining output is done to stdout. To redirect to a log file, try:
```
python -u example.py > log &
```

Example
-------
An example calculation is included with this distribution and is coded in example.py

