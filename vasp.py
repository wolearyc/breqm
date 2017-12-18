# Author: Willis O'Leary

import re
import os
import itertools
import numpy as np
from atoms import *

run_dir = 'vasp_run'  # Running directory

def calc_forces(atoms, regions):
    """Sets up and runs VASP calculation via a vasp.run script in the 
    run directory. Upon completion, updates atom forces. 
    """
    # POSCAR
    write_poscar(atoms,'{0}/POSCAR'.format(run_dir))

    #cmd = 'cd {0} && ./vasp.run'.format(run_dir)
    #p = Popen(cmd, stderr=STDOUT, stdout=PIPE,  shell=True)
    #jobID = p.communicate()[0]
    
    # A flag to indicate calc_forces was called to vasp.pbs // Collin
    os.system("echo 'flag' > vasp_run/vasp.flag")

    # For PBS: wait until job completes
    #block_pbs(jobID)
    block = True
    while block:
        time.sleep(1)
        block = os.path.isfile("vasp_run/vasp.flag") # flag does not exist when VASP finishes

    updt(atoms, run_dir, regions)

def write_poscar(atoms, file):
    """Writes data to a file in POSCAR format. Uses selective dynamics in 
    cartesian coordinate system. Does NOT write velocities."""
    f = open(file, 'w')
    # Comment and scaling factor
    f.write('Generated POSCAR\n1.0\n')
    # Cell matrix
    f.write('{0} {1} {2}\n'.format(*atoms.cell[0]))
    f.write('{0} {1} {2}\n'.format(*atoms.cell[1]))
    f.write('{0} {1} {2}\n'.format(*atoms.cell[2]))
    # Element list and number
    for element in atoms.elements:
        f.write(element + ' ')
    f.write('\n')
    for num in atoms.num_elements:
        f.write(str(num) + ' ')
    # Selective dynamics and cartesian
    f.write('\nSelective Dynamics\nCartesian\n')
    # Position entries
    for i, position in enumerate(atoms.positions):
        selection = 'T T T'
        if atoms.fixed[i]:
            selection = 'F F F'
        x,y,z = position
        f.write('  {0:.5f} {1:.5f} {2:.5f} {3}\n'.format(x,y,z, selection))
    f.write('\n')    
    f.close()

def read_poscar(file):
    """Returns Atoms object containing data read from a POSCAR file. Does not
    read in velocities, instead setting the velocity of each atom to None."""
    symbols       = []
    positions     = []
    velocities    = []
    cell          = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])


    f = open(file, 'r')
    # Comment and scaling factor
    next(f)
    scale = float(next(f))
    # Cell matrix
    cell[0] = np.array([float(n)*scale for n in next(f).split()])
    cell[1] = np.array([float(n)*scale for n in next(f).split()])
    cell[2] = np.array([float(n)*scale for n in next(f).split()])
    # Element list and number
    elements = next(f).split()
    num_elements = [int(n) for n in next(f).split()]
    # Selective dynamics and coordinate system
    selective = False
    direct = False
    line = next(f)
    if 'Selective' in line:
        selective = True
        line = next(f)
        if 'Direct' in line:
            direct = True
    elif 'Direct' in line:
        direct = True
    # Atom positions
    for _ in xrange(sum(num_elements)):
        entries = next(f).split()
        r = np.array(tuple(float(n) for n in entries[:3]))
        if direct:
            r = cell.dot(r)
        positions.append(r)
    
    f.close()

    velocities = [None] * sum(num_elements)
    return Atoms(elements, num_elements, positions, velocities, cell)

def updt(atoms, dir, regions):
    """Updates Atoms object using the CONTCAR and OUTCAR in a directory. Reads
    in forces and positions (in case of wrap around), but ignores velocity data.
    """
    f = open('{0}/CONTCAR'.format(dir), 'r')
    for _ in xrange(7):
        next(f)
    # Selective dynamics and coordinate system
    selective = False
    direct = False
    line = next(f)
    if 'Selective' in line:
        selective = True
        line = next(f)
        if 'Direct' in line:
            direct = True
    elif 'Direct' in line:
        direct = True
    # Atom positions
    for i in xrange(len(atoms)):
        entries = next(f).split()
        r = np.array(tuple(float(n) for n in entries[:3]))
        if direct:
            r = atoms.cell.dot(r)
        atoms.positions[i] = r
    f.close()

    atoms.forces = [None] * len(atoms)
    f = open('{0}/OUTCAR'.format(dir), 'r')
    line = next(f)
    while line[:9] != ' POSITION':
        line = next(f)
    next(f)
    for i in xrange(len(atoms)):
        if not regions.fixed_qm(*atoms.positions[i]):
            SplitLine = list(filter(None, re.split(r"\ |(-|[()])", next(f))))
            FixedSplitLine = []
            for j in range(len(SplitLine)):
                if SplitLine[j] == '-':
                    SplitLine[j+1] = -1 * float(SplitLine[j+1])
                else:
                    FixedSplitLine.append(float(SplitLine[j]))
            nums = FixedSplitLine #[float(s) for s in splitLine]
            atoms.forces[i] = np.array(nums[3:]) * 0.00964895
    f.close()

