# Author: Willis O'Leary

from threading import Thread
import datetime
from abc import ABCMeta, abstractmethod
from atoms import *
import lammps
import vasp

class Regions:
    """Abstract class that specifies BREQM regions. Regions are specified 
    by methods which return True or False, depending on whether the atom
    positions passed to the method are contained within a certain region.""" 
    __metaclass__ = ABCMeta
    def __init__(self, cutoffs):
        self.cutoffs = cutoffs

    @abstractmethod
    def qm(self, x,y,z):
        """Returns True if atom should be included in QM calculation."""
        pass
    @abstractmethod
    def mm(self, x,y,z):
        """Returns True if atom should be included in MM calculation."""
        pass
    @abstractmethod
    def fixed_qm(self, x,y,z):
        """Returns True if atom should be fixed in QM calculation."""
        pass
    @abstractmethod
    def fixed_mm(self, x,y,z):
        """Returns True if atom should be fixed in MM calculation."""
        pass
    @abstractmethod
    def mm_boundary(self, x,y,z):
        """Defines general mm calculation boundary to search for bonding."""
        pass
    @abstractmethod
    def qm_boundary(self, x,y,z):
        """defines general qm calculation boundary to search for bonding."""
        pass
    def boundary(self, x,y,z):
        """Returns True if atom is within the MM/QM boundary."""
        return self.qm(x,y,z) and self.mm(x,y,z)

class ZIntervalRegions(Regions):
    """BREQM regions specified by z-cutoffs. This region specification is
    useful for slab/solvent systems.
    """
    def __init__(self, qm_interval, mm_interval, midpoint, cutoffs):
        self.cutoffs = cutoffs
        self.qm_z1, self.qm_z2 = qm_interval
        self.mm_z1, self.mm_z2 = mm_interval
        self.midpoint          = midpoint

    def qm(self,x,y,z):
        return self.qm_z1 <= z <= self.qm_z2

    def mm(self,x,y,z):
        return self.mm_z1 <= z <= self.mm_z2

    def fixed_qm(self,x,y,z):
        if self.qm_z1 < self.mm_z1:
            return z >= self.midpoint
        else:
            return z <= self.midpoint

    def fixed_mm(self,x,y,z):
        if self.qm_z1 > self.mm_z1:
            return z >= self.midpoint
        else:
            return z <= self.midpoint

    def mm_boundary(self, x,y,z):
        return self.mm_z1-1.5<=z<=self.mm_z1+1.5 or self.mm_z2-1.5<=z<=self.mm_z2+1.5

    def qm_boundary(self, x,y,z):
        return self.qm_z1-1.5<=z<=self.qm_z1+1.5 or self.qm_z2-1.5<=z<=self.qm_z2+1.5

def append_trj(atoms):
    with open('trj.xyz', 'a') as f:
        f.write('{0}\ntrajectory\n'.format(len(atoms)))
        for i, element in atoms.enumerate():
            f.write('{0} {1} {2} {3}\n'.format(element, *atoms.positions[i])) 

def clear_trj():
    f = open('trj.xyz', 'w')
    f.close()

def run_heating_md(atoms, regions, dt, mm_ratio, steps, T_i, T_f, freq):
    """Performs a heating MD equilibration, using QM timestep dt and MM timestep
     dt*dt_ratio. System is heated from T_i to T_f in steps (QM). Velocity is 
     rescaled every freq steps. Between rescalings, an NVE ensemble is used 
     (velocity verlet). 
    """
    print """
========================================
============== HEATING MD ==============
========================================
Temperature:  {0:.1f} K --> {1:.1f} K
QM Timestep:  {2:.2f} fs
MM Timestep:  {3:.2f} fs

Running {4} steps.
Rescaling velocities every {5} steps.
Starting...""".format(T_i, T_f, dt, float(dt)/mm_ratio, steps, freq)
    
    atoms.init_velocities(T_i)
    clear_trj()
    
    atoms.update_bonding(regions.cutoffs, regions.qm_boundary)
    qm = atoms.split(regions.qm)
    qm.fix(regions.fixed_qm)

    for step in xrange(1, steps+1):
        print ' Step {0} '.format(step).center(40, '-')
        # Split systems and apply fixes
        atoms.update_bonding(regions.cutoffs, regions.mm_boundary)
        mm = atoms.split(regions.mm)
        mm.fix(regions.fixed_mm)
        atoms.update_bonding(regions.cutoffs, regions.qm_boundary)
        qm = atoms.split(regions.qm)
        qm.fix(regions.fixed_qm)
        # Determine shared atoms (to exclude from MM forces)
        overlap_indexes = []
        mm_i, qm_i = 0, 0
        while mm_i < len(mm) and qm_i < len(qm):
            mm_merge_i = mm.merge_indexes[mm_i]
            qm_merge_i = qm.merge_indexes[qm_i]
            if mm_merge_i == qm_merge_i:
                overlap_indexes.append(mm_i+1) # add 1 b/c LAMMPS
                qm_i += 1
                mm_i += 1
            elif mm_merge_i > qm_merge_i:
                qm_i += 1
            else:
                mm_i += 1

        qmThread = Thread(target=vasp.calc_forces, args=(qm,))
        if step % freq == 0:
            qm.rescale_velocities(T_i + (T_f-T_i) * float(step)/steps)
        print 'Starting QM...'
        qmThread.start()

        for i in xrange(mm_ratio):
            # Rescale velocities based on MM steps
            print 'Starting MM...',
            if (step * mm_ratio + i) % freq == 0:
                mm.rescale_velocities(T_i + (T_f-T_i) * float(step)/steps)
            lammps.calc_forces(mm, overlap_indexes)
            mm.step_nve(float(dt)/mm_ratio)
            print 'done.'

        qmThread.join()
        print 'QM done.'
        atoms.merge(mm)
        atoms.merge(qm)
        atoms.step_nve(float(dt))

        print '   QM Temperature: {0:.1f} K'.format(qm.temp())
        print '   MM Temperature: {0:.1f} K'.format(mm.temp())
        print 'Total Temperature: {0:.1f} K'.format(atoms.temp())
        append_trj(atoms)

def run_nvt_md(atoms, regions, dt, mm_ratio, steps, T, freq):
    """Performs a nvt, Nose-Hoover MD equilibration, using QM timestep dt and MM
     timestep dt/mm_ratio. System is kept at temperature T using a nose
     frequency freq. 
    """
    print """
========================================
================ NVT MD ================
========================================
Temperature:  {0:.1f} K
QM Timestep:  {1:.2f} fs
MM Timestep:  {2:.2f} fs

Running {3} steps.
Rescaling velocities every {4} steps.
Starting...""".format(T, dt, float(dt)/mm_ratio, steps, freq)
    
    atoms.init_velocities(T)
    clear_trj()
    
    atoms.update_bonding(regions.cutoffs, regions.qm_boundary)
    qm = atoms.split(regions.qm)
    qm.fix(regions.fixed_qm)

    for step in xrange(1, steps+1):
        print ' Step {0} '.format(step).center(40, '-')
        # Split systems and apply fixes
        atoms.update_bonding(regions.cutoffs, regions.mm_boundary)
        mm = atoms.split(regions.mm)
        mm.fix(regions.fixed_mm)
        atoms.update_bonding(regions.cutoffs, regions.qm_boundary)
        qm = atoms.split(regions.qm)
        qm.fix(regions.fixed_qm)
        # Determine shared atoms (to exclude from MM forces)
        overlap_indexes = []
        mm_i, qm_i = 0, 0
        while mm_i < len(mm) and qm_i < len(qm):
            mm_merge_i = mm.merge_indexes[mm_i]
            qm_merge_i = qm.merge_indexes[qm_i]
            if mm_merge_i == qm_merge_i:
                overlap_indexes.append(mm_i+1) # add 1 b/c LAMMPS
                qm_i += 1
                mm_i += 1
            elif mm_merge_i > qm_merge_i:
                qm_i += 1
            else:
                mm_i += 1

        qmThread = Thread(target=vasp.calc_forces, args=(qm,))
        print 'Starting QM...'
        qmThread.start()

        for i in xrange(mm_ratio):
            # Rescale velocities based on MM steps
            print 'Starting MM...',
            lammps.calc_forces(mm, overlap_indexes)
            mm.step_nvt(float(dt)/mm_ratio)
            print 'done.'

        qmThread.join()
        print 'QM done.'
        atoms.merge(mm)
        atoms.merge(qm)
        atoms.step_nve(float(dt))

        print '   QM Temperature: {0:.1f} K'.format(qm.temp())
        print '   MM Temperature: {0:.1f} K'.format(mm.temp())
        print 'Total Temperature: {0:.1f} K'.format(atoms.temp())
        append_trj(atoms)
