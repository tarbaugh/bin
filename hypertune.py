import flare

from flare import gp, struc, output, predict, md, otf, env,\
  otf_parser
from flare.kernels import mc_simple
from flare.utils import md_helper

import numpy as np
import time
import matplotlib.pyplot as plt
import pickle
from ase.visualize import view

from ase.io import read

dft_atoms = read('otf_dft.xyz', index=':')

# valid_species = valid.get_atomic_numbers()
# valid_positions = valid.get_positions()
# valid_cell = valid.get_cell()
# valid_forces = valid.get_forces()

# training_structure = struc.Structure(cell, species, positions)
# validation_structure = struc.Structure(cell, species, valid_positions)

kernels = ['twobody', 'threebody']
component = 'mc'
hyps = np.array([0.1, 0.1, 0.1, 2.0, 0.5])  # initial (bad) choice of hyps
cutoffs = {'twobody': 6.0, 'threebody': 3.5}  # cutoff radii in A
maxiter = 100  # max number of hyperparameter optimziation steps

gp_model = gp.GaussianProcess(
  kernels=kernels,
  component=component,
  hyps=hyps,
  cutoffs=cutoffs,
  maxiter=50
)

for snapshot in range(5):
    species = dft_atoms[snapshot].get_atomic_numbers()
    positions = dft_atoms[snapshot].get_positions()
    cell = dft_atoms[snapshot].get_cell()
    forces = dft_atoms[snapshot].get_forces()

    training_structure = struc.Structure(cell, species, positions)

    # add the structure to the training set of the GP
    gp_model.update_db(training_structure, forces)

gp_model.set_L_alpha()

print(gp_model.likelihood)

gp_model.train(print_progress=True)