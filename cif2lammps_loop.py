import argparse
import dpdata

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

REPLICATION = 3
scale = 0.95 #.989
i = 0

while scale <= 1.25:
    d_lmp = dpdata.System(file, fmt = 'vasp/poscar').replicate((REPLICATION,REPLICATION,REPLICATION))
    d_lmp['coords'][0] = scale*d_lmp['coords'][0]
    d_lmp['cells'][0] = scale*d_lmp['cells'][0]
    print(d_lmp.data)
    d_lmp.to('lammps/lmp', 'conf{}.lmp'.format(i))
    i += 1
    scale += 0.0125

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