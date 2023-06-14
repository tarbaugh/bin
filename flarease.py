import numpy as np
from ase import units
from ase.spacegroup import crystal
from ase.build import bulk
from ase.io import read
from ase import Atoms
import dpdata
np.random.seed(12345)

file = 'POSCAR.mp-938_GeTe'
REPLICATION = 2

perturbed_system = dpdata.System(file, fmt='vasp/poscar').replicate((REPLICATION,REPLICATION,REPLICATION)).perturb(pert_num=1,
    cell_pert_fraction=0,
    atom_pert_distance=.25,
    atom_pert_style='normal')

masses = []
for i in perturbed_system['atom_types']:
    mass_dict = {0: 72.63, 1: 127.603}
    masses.append(mass_dict[int(i)])

symbols = []
for i in perturbed_system['atom_types']:
    symbols_dict = {0: 'Ge', 1: 'Te'}
    symbols.append(symbols_dict[int(i)])

super_cell = Atoms(
    positions=perturbed_system['coords'][0],
    numbers=perturbed_system['atom_types'],
    masses=masses,
    cell=perturbed_system['cells'][0],
    pbc=(True, True, True)
)

# a = 3.52678
# super_cell = bulk('Ge', 'diamond', a=a, cubic=True)
# super_cell = read('data_6', format='lammps-data')

# super_cell = read('data.inp', format='vasp')
# super_cell.set_pbc((True,True,True))
# print(super_cell.get_tags())

from flare.gp import GaussianProcess
from flare.utils.parameter_helper import ParameterHelper

# set up GP hyperparameters
kernels = ['twobody', 'threebody'] # use 2+3 body kernel
parameters = {'cutoff_twobody': 6.0,
              'cutoff_threebody': 3.5}
pm = ParameterHelper(
    kernels = kernels,
    random = True,
    parameters=parameters
)

hm = pm.as_dict()
hyps = hm['hyps']
cut = hm['cutoffs']
print('hyps', hyps)

gp_model = GaussianProcess(
    kernels = kernels,
    component = 'mc', # If you are using ASE, please set to "mc" no matter for single-component or multi-component
    hyps = hyps,
    cutoffs = cut,
    hyp_labels = ['sig2','ls2','sig3','ls3','noise'],
    opt_algorithm = 'L-BFGS-B',
    n_cpus = 1
)

# from flare.mgp import MappedGaussianProcess

# grid_params = {'twobody':   {'grid_num': [64]},
#                'threebody': {'grid_num': [16, 16, 16]}}

# mgp_model = MappedGaussianProcess(grid_params,
#                                   unique_species = [6],
#                                   n_cpus = 1,
#                                   var_map=None)

from flare.ase.calculator import FLARE_Calculator

flare_calculator = FLARE_Calculator(gp_model,
                                    par = True,
                                    mgp_model = None,
                                    use_mapping = False)

super_cell.set_calculator(flare_calculator)

from ase.calculators.lammpsrun import LAMMPS

lmp_params = {'pair_style': 'quip',
              'pair_coeff': ['* * /zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml "Potential xml_label=GAP_2017_5_7_60_20_48_12_148" 32 52']}
files = ['/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml', 
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1481',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1482',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1483',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1484',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1485',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1486',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1487',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1488',
    '/zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml.sparseX.GAP_2017_5_7_60_20_48_12_1489',]

lammps = LAMMPS(parameters=lmp_params, files=files)
lammps.set(**lmp_params)

# import os
# from ase.calculators.cp2k import CP2K

# # ---------------- set up executable ----------------
# label = 'C'
# input_file = label+'.inp'
# output_file = label+'.out'
# no_cpus = 4
# pw_loc = '/zfshomes/tarbaugh/downloads/n79cp2k/cp2k-9.1/exe/local/cp2k_shell.ssmp'

# # serial
# os.environ['ASE_CP2K_COMMAND'] = f'{pw_loc} < {input_file} > {output_file}'

# ## parallel qe using mpirun
# # os.environ['ASE_ESPRESSO_COMMAND'] = f'mpirun -np {no_cpus} {pw_loc} -npool {npool} < {input_file} > {output_file}'

# ## parallel qe using srun (for slurm system)
# # os.environ['ASE_ESPRESSO_COMMAND'] = 'srun -n {no_cpus} --mpi=pmi2 {pw_loc} -npool {npool} < {input_file} > {output_file}'


# # -------------- set up input parameters --------------
# input_data = {'control':   {'prefix': label,
#                             'pseudo_dir': './',
#                             'outdir': './out',
#                             'calculation': 'scf'},
#               'system':    {'ibrav': 0,
#                             'ecutwfc': 60,
#                             'ecutrho': 360},
#               'electrons': {'conv_thr': 1.0e-9,
#                             'electron_maxstep': 100,
#                             'mixing_beta': 0.7}}

# # ----------------  pseudo-potentials -----------------
# ion_pseudo = {'C': 'C.pz-rrkjus.UPF'}

# # -------------- create ASE calculator ----------------
# dft_calc = Espresso(pseudopotentials=ion_pseudo, label=label,
#                     tstress=True, tprnfor=True, nosym=True,
#                     input_data=input_data, kpts=(8, 8, 8))

from ase import units
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution,
                                         Stationary, ZeroRotation)

MaxwellBoltzmannDistribution(super_cell, temperature_K=943)
Stationary(super_cell)  # zero linear momentum
ZeroRotation(super_cell)  # zero angular momentum

md_engine = 'VelocityVerlet'
md_kwargs = {}

from flare.ase.otf import ASE_OTF

otf_params = {
    'output_name': 'otf',
    'std_tolerance_factor': -0.01,
    'max_atoms_added' : 10,
    'freeze_hyps': 10,
    'write_model': 3 # If you will probably resume the training, please set to 3
}

test_otf = ASE_OTF(super_cell,
                   timestep = 1 * units.fs,
                   number_of_steps = 50,
                   dft_calc = lammps,
                   md_engine = md_engine,
                   md_kwargs = md_kwargs,
                   **otf_params)

test_otf.run()