from ase.io import read
from ase.calculators.espresso import Espresso

from scipy.spatial import distance_matrix

import numpy as np

import glob
import sys
import os

pseudopotentials = {'Sb': 'Sb.pbe-n-kjpaw_psl.1.0.0.UPF'}

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

NUM_CFGS = 100
SHEAR = 0.01
SCALE = 0.1

mp_cifs = glob.glob('*.cif')
if len(mp_cifs) > 1:
    print("Two cifs!")
    exit(1)

mp_atoms = read(mp_cifs[0])

for i in range(NUM_CFGS):
    try:
        os.makedirs(str(i))
    except:
        continue
    os.chdir(str(i))

    atoms = mp_atoms.repeat((2, 2, 3))

    old_p = atoms.get_positions()

    distances = distance_matrix(old_p, old_p)
    np.fill_diagonal(distances, np.inf)
    smallest_distance = np.min(np.min(distances))
    DISP = 0.125*smallest_distance

    old_p += 2*DISP*np.random.random_sample(old_p.shape)-DISP
    atoms.set_positions(old_p)

    atoms.set_pbc((True,True,True))
    atoms.wrap()

    shear_values = 2*SHEAR*np.random.random_sample(size=3)-SHEAR
    scale_values = 2*SCALE*np.random.random_sample(size=3)-SCALE + 1
    shear_matrix = np.array([[scale_values[0], shear_values[0], shear_values[1]],
                            [0, scale_values[1], shear_values[2]],
                            [0, 0, scale_values[2]]])

    atoms.set_cell(np.dot(atoms.get_cell(), shear_matrix))
    atoms.set_positions(np.dot(atoms.get_positions(), shear_matrix))

    atoms.set_pbc((True,True,True))
    atoms.wrap()

    calc = Espresso(pseudopotentials=pseudopotentials, kspacing=0.03, input_data=input_data)

    atoms.calc = calc
    atoms.get_potential_energy()

    atoms.write('atoms_out.extxyz')
    os.chdir('..')
