import random
import numpy as np
from ase.lattice.cubic import SimpleCubic

ase_atoms = SimpleCubic(
        symbol='Te',
        size=(6, 6, 6),
        pbc=True,
        latticeconstant=3.28537)

old_p = ase_atoms.get_positions()
old_p += 1*np.random.random_sample(old_p.shape)-0.5
ase_atoms.set_positions(old_p)
ase_atoms.set_pbc((True,True,True))
ase_atoms.wrap()

chems = ase_atoms.get_chemical_symbols()

indexs = [x for x in range(len(chems))]
random.shuffle(indexs)

for x in range(32):
	chems[indexs[x]] = 'Ge'
ase_atoms.set_chemical_symbols(chems)
ase_atoms.write('data', format='lammps-data')
