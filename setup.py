from ase.io import read
from ase.calculators.espresso import Espresso

from scipy.spatial import distance_matrix

import numpy as np

import glob
import sys
import os

pseudopotentials = {'Sb': 'Sb.pbe-n-kjpaw_psl.1.0.0.UPF',
                    'Ge': 'Ge.pbe-n-kjpaw_psl.1.0.0.UPF'}

input_data = {
  'control': {
    'disk_io': 'none',
    'tprnfor': True,
    'pseudo_dir': '/zfshomes/tarbaugh/downloads/qe-7.2_intel/pseudo'
  },
  'system': {
    'ecutwfc': 30,
    'ecutrho': 150,
    'occupations': 'smearing',
    'smearing': 'gaussian',
    'degauss': 0.015
  }
}

atoms_list = read('quench.xyz', format='lammps-dump-text', index=':')

for i in range(100):
    try:
        os.makedirs(str(i))
    except:
        continue
    os.chdir(str(i))

    atoms = atoms_list[i*18]
    syms = atoms.get_chemical_symbols()
    newsyms = []
    for s in syms:
        if s == 'H':
            newsyms.append('Ge')
        if s == 'He':
            newsyms.append('Sb')
    
    atoms.set_chemical_symbols(newsyms)

    atoms.wrap()

    calc = Espresso(pseudopotentials=pseudopotentials, kspacing=0.03, input_data=input_data)

    atoms.calc = calc
    atoms.get_potential_energy()

    atoms.write('atoms_out.extxyz')
    os.chdir('..')
