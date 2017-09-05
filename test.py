#!/home/wolearyc/.local/bin/python
# Author: Willis O'Leary

from atoms import *
from breqm import *

# define bond cutoff dictionary - this is used when getting bonding information
cutoffs = {'H-O' : 1.5}
total = vasp.read_poscar('input.vasp')
regions = ZIntervalRegions((0,13), (9,40), 11, cutoffs)
run_heating_md(total, regions, 1.2, 5, 1666, 50, 300, 20)

