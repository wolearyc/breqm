# Author: Willis O'Leary

from atoms import *
from breqm import *

# Define relevant bond cutoff dictionary (this one's good for water)
cutoffs = {'H-O' : 1.5}
# Read in input file (vasp format)
total = vasp.read_poscar('example_input.vasp')
# Define BREQM regions. In this example, we separate regions using
# z intervals. Interval 0-13 is included in QM calculations. 
# Interval 9-40 is included in MM calcuations. 9-11 is the 
# QM boundary, while 11-13 is the MM boundary.
regions = ZIntervalRegions((0,13), (9,40), 11, cutoffs)

# Run 2 ps heating dynamics
run_heating_md(total, regions, 1.2, 5, 1666, 50, 300, 1.0)

# Run 2 ps NVT dynamics
# run_nvt_md(total, regions, 1.2, 5, 1666, 300.0, 1.0)
