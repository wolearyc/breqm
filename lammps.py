# Author: Willis O'Leary

import os
import itertools
import numpy as np
from util import *

run_dir = 'lammps_run' # Running directory

def calc_forces(atoms, regions):
    """Sets up and runs LAMMPS calculation via a lammps.run script in the 
    run directory. Upon completion, updates atom forces.
    """
    # data file
    write_data(atoms, '{0}/lammps.data'.format(run_dir))

    cmd = 'cd {0} && ./lammps.run'.format(run_dir)
    p = Popen(cmd, stderr=STDOUT, stdout=PIPE, shell=True)
    jobID = p.communicate()[0]

    # For PBS: wait until job completes
    #block_pbs(jobID)
    #block = True
    #while block:
    #    time.sleep(1)
    #    block = os.path.isfile("vasp_run/vasp.flag") # flag does not exist when VASP finishes


    updt(atoms, run_dir, regions)

def write_data(atoms, file):
    """Writes data to a file in lammps data format."""
    f = open(file, 'w')
    # Header
    f.write('# Generated data file\n\n')
    f.write('{0} atoms\n'.format(sum(atoms.num_elements)))
    f.write('{0} atom types\n\n'.format(len(atoms.elements)))
    # Simulation box
    f.write('0.0 {0}  xlo xhi\n'.format(atoms.cell[0][0]))
    f.write('0.0 {0}  ylo yhi\n'.format(atoms.cell[1][1]))
    f.write('0.0 {0}  zlo zhi\n'.format(atoms.cell[2][2]))
    xy, xz, yz = atoms.cell[1][0], atoms.cell[2][0], atoms.cell[2][1]
    f.write('{0} {1} {2}  xy xz yz\n\n'.format(xy, xz, yz))
    # IDs and Masses
    f.write('Masses\n\n')
    atom_id = 1
    for element in atoms.elements:
        f.write('{0} {1} # {2}\n'.format(atom_id, masses[element], element))
        atom_id += 1
    f.write('\n')
    # Atoms
    f.write('Atoms\n\n')
    index = 0
    for atom_id in range(1, len(atoms.elements)+1):
        for i in range(atoms.num_elements[atom_id-1]):
            f.write('{0:>6} {1:>3} 0.0 {2} {3} {4}\n'.format(index+1, atom_id, 
                *atoms.positions[index]))
            index += 1

    f.close()

def updt(atoms, dir, regions):
    """Updates data with that read from a LAMMPS trajectory file. This function
    reads the first structure entry, and expects the following dump format:
    'id type x y z fx fy fz'. 
    """
    if atoms.forces is None:
        atoms.forces = [None] * len(atoms)  

    try:
        f = open('{0}/lammps.trj'.format(dir), 'r')
        for _ in xrange(9):
            next(f)
        for i in xrange(len(atoms)):
            nums = [float(n) for n in next(f).split()]
            atoms.positions[i] = np.array(nums[2:5])
            if not regions.fixed_mm(*atoms.positions[i]):       
                atoms.forces[i] = np.array(nums[5:8]) * 4.184e-4
        f.close()
    except:
        f = open('{0}/lammps.data'.format(dir), 'r')
        for _ in xrange(17):
            next(f)
        for i in xrange(len(atoms)):
            nums = [float(n) for n in next(f).split()]
            atoms.positions[i] = np.array(nums[3:])
        f.close()




