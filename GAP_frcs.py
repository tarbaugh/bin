import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.calculators.lammpsrun import LAMMPS

lmp_params = {'pair_style': 'quip',
              'pair_coeff': ['* * /zfshomes/tarbaugh/downloads/GAP_GST/gp_merged_2bdmbd.xml "Potential xml_label=GAP_2017_5_7_60_20_48_12_148" 32']}
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

cell = np.eye(3) * 1000
seps = np.arange(2, 6, 0.1)
# specs = [[0, 0], [0, 1], [1, 1]]
specs = [[0, 0], [1, 1]]
# masses = [[72.63, 72.63], [72.63, 127.603], [127.603, 127.603]]
masses = [[72.63, 72.63], [127.603, 127.603]]
frcs = np.zeros((3, len(seps)))

for m, sep in enumerate(seps):
    pos = np.array([[0, 0, 0], [sep, 0, 0]])
    print(sep)
    for n, spec_curr in enumerate(specs):
        atoms = Atoms(positions=pos, numbers=[0, 0], cell=cell, masses=masses[0], calculator=lammps)

        # predict the x component of the force

        frcs[n, m] = atoms.get_forces()[1][0]

# plot GP predictions vs ground truth
cols = ['b', 'r', 'g']
for n in range(2):
    plt.plot(seps, frcs[n], color=cols[n], linestyle='--')

plt.xlabel('separation ($\AA$)')
plt.ylabel('force (eV/$\AA$)')
plt.xlim(2, 6.0)
plt.ylim(-10, 15)
plt.savefig('fig.png')