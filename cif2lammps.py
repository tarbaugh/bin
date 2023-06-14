import argparse
import dpdata

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

REPLICATION = 3
SCALE = 0.95 #.989

for i in range(100):
    # d_lmp = dpdata.System(file, fmt = 'vasp/poscar').replicate((REPLICATION,REPLICATION,REPLICATION))
    # d_lmp['coords'][0] = (SCALE+i*0.002)*d_lmp['coords'][0]
    # d_lmp['cells'][0] = (SCALE+i*0.002)*d_lmp['cells'][0]
    perturbed_system = dpdata.System(file, fmt='vasp/poscar').replicate((REPLICATION,REPLICATION,REPLICATION)).perturb(pert_num=1,
        cell_pert_fraction=0,
        atom_pert_distance=.25,
        atom_pert_style='normal')
    perturbed_system['coords'][0] = (SCALE+i*0.002)*perturbed_system['coords'][0]
    perturbed_system['cells'][0] = (SCALE+i*0.002)*perturbed_system['cells'][0]
    perturbed_system.to('vasp/poscar', './938_data/data{}.inp'.format(i))

# Tellurium Notes:
# Scaling factor : density
# 0.95 = 7.05
# 0.985 = 6.3291
# 0.989 = 6.2526
# 0.99 = 6.2336
# 1.0 = 6.048
# 1.05 = 5.224
# 1.1 = 4.544
# Lower densitites show HCP structure