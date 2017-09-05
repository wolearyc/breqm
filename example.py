# Author: Willis O'Leary

from atoms import *
from breqm import *

# Define relevant bond cutoff dictionary (this one's good for water)
cutoffs = {'H-O' : 1.5}
# Read in input file (vasp format)
total = vasp.read_poscar('example_input.vasp')
# Define BREQM regions
regions = ZIntervalRegions((0,13), (9,40), 11, cutoffs)

# Run 2 ps heating dynamics
# run_heating_md(total, regions, 1.2, 5, 1666, 50, 300, 20)

# Run 2 ps NVT dynamics
# run_nvt_md(total, regions, 1.2, 5, 1666, 300.0, 1.0)
