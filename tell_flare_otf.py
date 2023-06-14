from asynchat import simple_producer
from selectors import EpollSelector
import numpy as np
from ase import units, Atoms
from ase.spacegroup import crystal
from ase.build import bulk


from flare import otf_parser, struc, env
from flare.gp import GaussianProcess
from flare.utils.parameter_helper import ParameterHelper
from flare.predict import predict_on_structure
import matplotlib.pyplot as plt

from flare.ase.otf import ASE_OTF

np.random.seed(12345)

a = 5.336
super_cell = bulk('Ge', 'diamond', a=a, cubic=True)
sc_pos = super_cell.get_positions()
print("Original pos:")
print(sc_pos)
for atom_pos in sc_pos:
    for coord in range(3):
        atom_pos[coord] += (2*np.random.random()-1) * 0.05
super_cell.set_positions(sc_pos)
super_cell.wrap()
print("Adujusted pos:")
print(super_cell.get_positions())
super_cell.set_pbc((True, True, True))

kernels = ['twobody']
hyps = np.array([0.92961609, 0.31637555, 0.05])
hyp_labels = ['Sigma', 'Length Scale', 'Noise']
cut = {'twobody': 7.0}

gp_model = GaussianProcess(
    kernels = kernels,
    component = 'mc', # If you are using ASE, please set to "mc" no matter for single-component or multi-component
    hyps = hyps,
    cutoffs = cut,
    hyp_labels = hyp_labels,
    opt_algorithm = 'L-BFGS-B',
    n_cpus = 12
)

from flare.ase.calculator import FLARE_Calculator

flare_calculator = FLARE_Calculator(gp_model,
                                    par = True,
                                    mgp_model = None,
                                    use_mapping = False)

super_cell.set_calculator(flare_calculator)

from ase.calculators.lj import LennardJones
lj_calc = LennardJones(rc=8.5526300, epsilon=2.7017100, sigma=2.1381600, smooth=True)

from ase import units
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution,
                                         Stationary, ZeroRotation)

MaxwellBoltzmannDistribution(super_cell, temperature_K=2000)
Stationary(super_cell)  # zero linear momentum
ZeroRotation(super_cell)  # zero angular momentum

md_engine = 'VelocityVerlet'
md_kwargs = {}

from flare.ase.otf import ASE_OTF

otf_params = {'init_atoms': [0, 1, 2, 3],
              'output_name': 'otf',
              'std_tolerance_factor': 2,
              'max_atoms_added' : 4,
              'freeze_hyps': 10,
              'write_model': 3} # If you will probably resume the training, please set to 3

test_otf = ASE_OTF(super_cell,
                   timestep = 1 * units.fs,
                   number_of_steps = 500,
                   dft_calc = lj_calc,
                   md_engine = md_engine,
                   md_kwargs = md_kwargs,
                   **otf_params,
                   n_cpus=12)

test_otf.run()

gp_model = test_otf.gp
gp_model.write_model('model_out.gp')

output_file = 'otf.out'
otf_trajectory = otf_parser.OtfAnalysis(output_file)

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
seps = np.arange(2, 7, 0.01)

store_frcs = np.zeros((1, len(seps)))
gp_frcs = np.zeros((1, len(seps)))
gp_stds = np.zeros((1, len(seps)))
cutoff = 7.0
n = 0
for m, sep in enumerate(seps):
    pos = np.array([[0, 0, 0], [sep, 0, 0]])
    atoms = Atoms(positions=pos, symbols=['Ge', 'Ge'], cell=cell, calculator=lj_calc)

    struc_curr = struc.Structure(cell, ['Ge', 'Ge'], pos)
    env_curr = env.AtomicEnvironment(struc_curr, 1, np.array([cutoff]))
    store_frcs[n, m] = atoms.get_forces()[1][0]
    print(store_frcs[n,m])

    # predict the x component of the force
    pred, var = gp_model.predict(env_curr, 1)
    gp_frcs[n, m] = pred
    print(gp_frcs[n,m])
    gp_stds[n, m] = np.sqrt(var)

# plot GP predictions vs ground truth
cols = ['b', 'r', 'g']

plt.plot(seps, store_frcs[n], color=cols[n], linestyle='-')
plt.plot(seps, gp_frcs[n], color=cols[n], linestyle='--')
plt.fill_between(seps, gp_frcs[n] + 3 * gp_stds[n],
                    gp_frcs[n] - 3 * gp_stds[n],
                    color=cols[n], alpha=0.4)

plt.xlabel('separation ($\AA$)')
plt.ylabel('force (eV/$\AA$)')
plt.xlim(2, 7.0)
plt.ylim(-10, 20)
plt.savefig('nice.png')