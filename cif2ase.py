import argparse
import dpdata

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

REPLICATION = 1
SCALE = 1.0 #.989


perturbed_system = dpdata.System(file, fmt='vasp/poscar').replicate((REPLICATION,REPLICATION,REPLICATION)).perturb(pert_num=1,
    cell_pert_fraction=0,
    atom_pert_distance=.25,
    atom_pert_style='normal')

print(perturbed_system['coords'][0])
print(perturbed_system['atom_types'])
# perturbed_system['coords'][0] = (SCALE)*perturbed_system['coords'][0]
# perturbed_system['cells'][0] = (SCALE)*perturbed_system['cells'][0]
# perturbed_system.to('vasp/poscar', './data.inp')

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