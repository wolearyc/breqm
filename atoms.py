# Author: Willis O'Leary

import random
import math
import itertools
import ase
import numpy as np
from scipy.stats import maxwell
from util import *


class Atoms:
    """A periodic system of atoms.
    
    My Atom's class, inspired by ASE's Atoms class, is designed with simplicity 
    and function in mind. Atoms can be split up into multiple Atoms, which can
    be modified and merged back into the original Atoms object. To increase QM 
    performance, atoms are organized by element, and relative order is always 
    maintained during splits and merges. 
    """
    def __init__(self, 
        elements, 
        num_elements, 
        positions, 
        velocities,
        cell):
        """Constructor.
        elements: list of contained element symbols in order of specification
        num_elements: number of contained elements correponding to 'elements'
        positions: cartesian tuple locations
        velocities: cartesian tuple velocities
        cell: unit cell matrix
        """
        self.elements       = elements
        self.num_elements   = num_elements
        self.positions      = positions
        self.cell           = cell
        self.velocities     = velocities
        self.forces         = None
        self.prev_forces    = None
        self.fixed          = [False] * sum(num_elements)
        self.merge_indexes  = None
        self.bond_info      = None
        self.zeta           = 0

    def __len__(self):
        """Returns number of atoms."""
        return sum(self.num_elements)

    def add_force(self, index, force):
        """Adds a force to a specific atom."""
        if force is not None:
            if self.forces is None:
                self.forces = [None] * len(self)
            if not self.forces[index]:
                self.forces[index] = force
            else:
                self.forces[index] += force

    def enumerate(self):
        """Enumerates the object for each index, element pair."""
        iterators = []
        for element, num in zip(self.elements, self.num_elements):
            iterators.append(itertools.repeat(element, num)) 
        return enumerate(itertools.chain(*iterators))

    def matches(self, condition, ignore_bonds=True):
        """Returns whether each atom matches a certain condition, optionally including atoms bonded to matched atoms. 
        """
        marked = [False] * len(self)
        def mark(i):
            if not marked[i]:
                marked[i] = True
                if self.bond_info is not None:
                    for other_i in self.bond_info[i]:
                        mark(other_i)

        for i, position in enumerate(self.positions):
            if condition(*position):
                mark(i)
        return marked

    def split(self, condition):
        """Returns new Atoms with atoms of interest (that satisfy condition). 
        Atoms bonded to atoms of interest are also included. 
        """
        elements       = []
        num_elements   = []
        positions      = []
        velocities     = []
        forces         = None if self.forces is None else []
        prev_forces    = None if self.prev_forces is None else []
        merge_indexes  = []

        marked = self.matches(condition, ignore_bonds=False)

        for i, element in self.enumerate():
            if marked[i]:
                if not elements or element is not elements[-1]:
                    elements.append(element)
                    num_elements.append(0)
                num_elements[-1] += 1
                positions.append(self.positions[i])
                velocities.append(self.velocities[i])
                if self.forces is not None:
                    forces.append(self.forces[i])
                if self.prev_forces is not None:
                    prev_forces.append(self.prev_forces[i])
                merge_indexes.append(i)

        if sum(num_elements) == 0:
            warn('0 atoms split')
        result = Atoms(elements, num_elements, positions, velocities, self.cell)
        result.merge_indexes = merge_indexes
        result.forces = forces
        result.prev_forces = prev_forces
        return result

    def merge(self, atoms):
        """Merges with another object assumed to be previously split. Clears 
        bonding and unfixes everything.
        """
        self.bond_info = None
        self.unfix_all()

        for i in xrange(len(atoms)):
            merge_i = atoms.merge_indexes[i]
            self.positions[merge_i] = atoms.positions[i]
            self.velocities[merge_i] = atoms.velocities[i]
            if atoms.forces is not None:
                self.add_force(merge_i, atoms.forces[i])
            if atoms.prev_forces is not None:
                if self.prev_forces is None:
                    self.prev_forces = [None] * len(self)
                self.prev_forces[merge_i] = atoms.prev_forces[i]

    def init_velocities(self, temp):
        """Initializes velocies randomly according to maxwell boltzmann 
        distribution. 
        """
        random.seed(None)
        for i, element in self.enumerate():
            scale = math.sqrt(kb * temp / masses[element])
            speed = maxwell.rvs(scale=scale, size=1)[0]
            theta = random.uniform(0, 360)
            phi = random.uniform(0, 180)
            vx = speed * math.sin(phi) * math.cos(theta)
            vy = speed * math.sin(phi) * math.sin(theta)
            vz = speed * math.cos(phi)
            self.velocities[i] = np.array((vx, vy, vz))

    def rescale_velocities(self, temp):
        """Rescales atom velocities randomly according to maxwell boltzmann
        distribution."""
        random.seed(None)
        for i, element in self.enumerate():
            scale = math.sqrt(kb * temp / masses[element])
            speed = maxwell.rvs(scale=scale, size=1)[0]
            self.velocities[i] *= speed/np.linalg.norm(self.velocities[i])

    def temp(self):
        """Returns simulation temperature."""
        s = 0
        for i, element in self.enumerate():
            v = self.velocities[i]
            s += masses[element] * v.dot(v)
        return 1/(kb * 3 * (len(self) - 1)) * s

    def step_nve(self, dt):
        """Performs MD step in an NVE ensemble using velocity verlet."""
        if not self.forces:
            error('cannot step without forces')
        for i, element in self.enumerate():
            if not self.fixed[i] and self.forces[i] is not None:
                m = masses[element]
                if self.prev_forces is not None and self.prev_forces[i] is not None:
                    # Steps 4 and 1
                    p_a = self.prev_forces[i]/m
                    a = self.forces[i]/m
                    self.velocities[i] += (p_a+a)*dt/2
                    # Step 2
                    self.positions[i] += self.velocities[i]*dt + a/2 * dt**2
                    # ...step 3...
                else:
                    a = self.forces[i]/m
                    # Step 2
                    self.positions[i] += self.velocities[i]*dt + a/2 * dt**2
                    # ...step 3...
        self.prev_forces = self.forces
        self.forces = None

    def step_nvt(self, dt, Q, temp):
        """Performs MD step in an NVT ensemble using a Nose-Hoover thermostat"""

        def nose_sum(m):
            result = 0
            for i, element in self.enumerate():
                m = masses[element]
                result += m * self.velocities[i].dot(self.velocities[i]) / 2
            result += -len(self) * (3*len(self) + 1)/2 * kb * temp
            return  dt/(2*Q) * result

        if not self.forces:
            error('cannot step without forces')

        for i, element in self.enumerate():
            if not self.fixed[i] and self.forces[i] is not None:
                m = masses[element]
                z = self.zeta
                if self.prev_forces is not None:
                    # Step 6 (velocity is out of date)
                    v = self.velocities[i]
                    p_a = self.prev_forces[i]/m
                    self.velocities[i] = (v + dt/2*p_a)/(1+dt/2*z)
                # Step 1
                a = self.forces[i]/m
                v = self.velocities[i]
                self.positions[i] += v*dt + (a-z*v) * dt**2 / 2
                # Step 4
                self.zeta += nose_sum(m)
                # Step 2
                self.velocities[i] += dt/2 * (a-z*v)
        for i, element in self.enumerate():
            if not self.fixed[i] and self.forces[i] is not None:
                m = masses[element]
                # Step 5
                self.zeta += nose_sum(m)
        # ...Step 3...
        self.prev_forces = self.forces
        self.forces = None

    def fix(self, condition):
        """Fixes atoms meeting a condition, paying no attention to bonding."""
        for i, position in enumerate(self.positions):
            if condition(*position):
                self.fixed[i] = True

    def unfix_all(self):
        """Unfixes all atoms."""
        for i in xrange(len(self)):
            self.fixed[i] = False

    def update_bonding(self, cutoffs, condition):
        """Clears and updates bonds using cutoff dictionary.
        Condition narrows atoms to evaluate, increasing speed.
        """
        self.bond_info = [set() for _ in xrange(len(self))]

        rel_indexes = []
        rel_positions = []
        rel_elements = []
        index = 0
        for i, element in enumerate(self.elements):
            for _ in xrange(self.num_elements[i]):
                if condition(*self.positions[index]):
                    rel_indexes.append(index)
                    rel_positions.append(self.positions[index])
                    rel_elements.append(element)
                index += 1

        ase_atoms = ase.Atoms(positions=rel_positions, cell=self.cell, pbc=[1,1,0])
        distance_info = ase_atoms.get_all_distances(mic=True)

        for distances, i1, e1 in zip(distance_info, rel_indexes, rel_elements):
            for distance, i2, e2 in zip(distances, rel_indexes, rel_elements):
                if i1 == i2:
                    continue
                bond = '{0}-{1}'.format(e1, e2) 
                cutoff = 0
                if bond in cutoffs:
                    cutoff = cutoffs[bond]
                elif bond[::-1] in cutoffs:
                    cutoff = cutoffs[bond[::-1]]
                else:
                    continue
                if distance <= cutoff:
                    self.bond_info[i1].add(i2)
                    self.bond_info[i2].add(i1)


