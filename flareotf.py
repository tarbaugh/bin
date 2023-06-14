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

def get_LJ_forces(structure, lj_parameters):
  """Calculate multicomponent Lennard Jones forces on a structure of atoms.
  dft_kwargs is assumed to be a dictionary containing the cutoff, an ordered
  list of species, and arrays containing epsilon and sigma values."""

  cutoff = lj_parameters['cutoff']
  epsilons = lj_parameters['epsilons']
  sigmas = lj_parameters['sigmas']
  spec_list = lj_parameters['species']

  forces = np.zeros((structure.nat, 3))
  # Loop over atoms in the structure.
  for m in range(structure.nat):
      # Create atomic environment.
      environment = env.AtomicEnvironment(structure, m, np.array([cutoff]))
      ind1 = spec_list.index(environment.ctype)

      # Loop over atoms in the environment to compute the total force.
      for n in range(len(environment.etypes)):
          ind2 = spec_list.index(environment.etypes[n])
          eps = epsilons[ind1, ind2]
          sig = sigmas[ind1, ind2]

          # Compute LJ force.
          bond_vals = environment.bond_array_2[n]
          r = bond_vals[0]

          dE_dr = 4 * eps * (-12 * sig ** 12 / (r ** 13) +
                             6 * sig ** 6 / (r ** 7))

          forces[m, 0] += dE_dr * bond_vals[1]
          forces[m, 1] += dE_dr * bond_vals[2]
          forces[m, 2] += dE_dr * bond_vals[3]

  return forces

class lj_module:
  def parse_dft_input(file_name):
      """We assume the input is a pickled dictionary containing positions
      (in angstrom), species (as a list of strings), cell (as a 3x3 matrix of
      cell vectors), and masses (as a dictionary assigning each species a mass
      in AMU)."""

      input_file = open(file_name, 'rb')
      struc_dict = pickle.load(input_file)
      input_file.close()

      # Convert masses to MD units (energy = eV, length = A, )
      masses = struc_dict['masses']
      mass_convert = 0.000103642695727
      for species in masses:
          masses[species] *= mass_convert

      return struc_dict['positions'], struc_dict['species'], \
          struc_dict['cell'], masses

  def run_dft_par(dft_input=None, structure=None, dft_loc=None, n_cpus=None,
                npool=None, mpi=None, dft_kwargs=None, dft_out=None):
    return get_LJ_forces(structure, dft_kwargs)

cutoff = 5.0
epsilons = np.array([[3, 2.5], [2.5, 3.]])
sigmas = np.array([[2.0, 1.9], [1.9, 2.1]])
spec_list = [47, 53]  # silver and iodine atomic numbers

lj_params = {'cutoff': cutoff, 'epsilons': epsilons, 'sigmas': sigmas,
             'species': spec_list}

alat = 2.54
unit_cell = np.eye(3) * alat

# define bcc positions
unit_positions = np.array([[0, 0, 0],
                           [1/2, 1/2, 1/2]]) * alat

# make a supercell (54 atoms)
sc_size = 3
positions = md_helper.get_supercell_positions(sc_size, unit_cell, unit_positions)
cell = unit_cell * sc_size

#  positions to give nonzero force on first frame
for atom_pos in positions:
    for coord in range(3):
        atom_pos[coord] += (2*np.random.random()-1) * 0.05

# create initial structure
species = ['Ag', 'I'] * sc_size ** 3
struc_curr = struc.Structure(cell, species, positions)

# create pseudo input file
mass_dictionary = {'Ag': 108, 'I': 127}
input_dictionary = {'positions': positions, 'species': species, 'cell': cell,
                    'masses': mass_dictionary}
input_file_name = 'lj.in'
with open(input_file_name, 'wb') as input_file:
    pickle.dump(input_dictionary, input_file)

kernels = ['twobody']
component = 'mc'
hyps = np.array([0.1, 1., 0.06])
hyp_labels = ['Sigma', 'Length Scale', 'Noise']
cutoffs = {'twobody': 5.0}

gp_model = gp.GaussianProcess(
  kernels=kernels,
  component=component,
  hyps=hyps,
  hyp_labels=hyp_labels,
  cutoffs=cutoffs,
  maxiter=1000
)

dt = 0.001  # time step (ps)
number_of_steps = 1000
dft_loc = None  # path to dft executable would usually go here
std_tolerance_factor = -0.01  # 10 meV/A
init_atoms = [0, 25, 50]  # initial atoms added to the training set
max_atoms_added = 5  # number of atoms added when dft is called
freeze_hyps = 5  # no hyperparameter optimization after this many updates

# rescale the temperature halfway through the simulation
rescale_steps = [10]  # rescale at step 10
rescale_temps = [1000]

otf_model = otf.OTF(
    # MD arguments
    dt,
    number_of_steps,
    rescale_steps=rescale_steps,
    rescale_temps=rescale_temps,
    # FLARE arguments
    gp=gp_model,
    # OTF arguments
    std_tolerance_factor=std_tolerance_factor,
    init_atoms=init_atoms,
    max_atoms_added=max_atoms_added,
    freeze_hyps=freeze_hyps,
    # DFT arguments
    force_source=lj_module,
    dft_loc=dft_loc,
    dft_input=input_file_name,
    dft_kwargs=lj_params,
    n_cpus=8
    )

otf_model.run()

output_file = 'otf_run.out'
otf_trajectory = otf_parser.OtfAnalysis(output_file)

gp_model.write_model('model_out.gp')

# plot temperature vs. simulation time
times = otf_trajectory.times
temps = otf_trajectory.thermostat['temperature']
dft_times = otf_trajectory.dft_times[1:]  # exclude t = 0

for n, dft_time in enumerate(dft_times):
    otf_ind = times.index(dft_time)
    plt.plot(dft_times[n], temps[otf_ind], 'kx')

plt.plot(times, temps)
plt.xlabel('time (ps)')
plt.ylabel('temperature (K)')
plt.savefig('ok.png')

cell = np.eye(3) * 1000
seps = np.arange(2, 5, 0.01)
specs = [[47, 47], [47, 53], [53, 53]]
store_frcs = np.zeros((3, len(seps)))
gp_frcs = np.zeros((3, len(seps)))
gp_stds = np.zeros((3, len(seps)))
for m, sep in enumerate(seps):
    pos = np.array([[0, 0, 0], [sep, 0, 0]])
    for n, spec_curr in enumerate(specs):
        struc_curr = struc.Structure(cell, spec_curr, pos)
        env_curr = env.AtomicEnvironment(struc_curr, 1, np.array([cutoff]))
        frcs = \
            lj_module.run_dft_par(structure=struc_curr, dft_kwargs=lj_params)
        store_frcs[n, m] = frcs[1, 0]

        # predict the x component of the force
        pred, var = gp_model.predict(env_curr, 1)

        gp_frcs[n, m] = pred
        gp_stds[n, m] = np.sqrt(var)

# plot GP predictions vs ground truth
cols = ['b', 'r', 'g']
for n in range(3):
    plt.plot(seps, store_frcs[n], color=cols[n], linestyle='-')
    plt.plot(seps, gp_frcs[n], color=cols[n], linestyle='--')
    plt.fill_between(seps, gp_frcs[n] + 3 * gp_stds[n],
                     gp_frcs[n] - 3 * gp_stds[n],
                     color=cols[n], alpha=0.4)

plt.xlabel('separation ($\AA$)')
plt.ylabel('force (eV/$\AA$)')
plt.xlim(2, 3.4)
plt.ylim(-10, 15)
plt.savefig('nice.png')