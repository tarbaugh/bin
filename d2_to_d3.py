from ase.calculators.dftd3 import DFTD3
from ase import Atoms
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

add = 0
 
with open(file, 'r') as f:
    tot_CFGS = 0
    for i, l in enumerate(f):
        if i == 2:
            num_A = int(l)
        if l.startswith('BEGIN_CFG'):
            tot_CFGS += 1

cells = np.zeros((tot_CFGS, 3,3)).astype(np.float64)
coords = np.zeros((tot_CFGS, num_A, 3)).astype(np.float64)
forces = np.zeros((tot_CFGS, num_A, 3)).astype(np.float64)
ids = np.zeros((tot_CFGS, num_A)).astype(int)
types = np.zeros((tot_CFGS, num_A)).astype(int)
energies = np.zeros(tot_CFGS).astype(np.float64)

with open(file, 'r') as f:
    curr_config = -1
    for _, l in enumerate(f):
        if l.startswith('BEGIN_CFG'):
            iter = 0
            curr_config += 1

        elif iter == 2:
            if int(l) != num_A:
                print("Error: number of atoms changed")
                exit(1)
        elif iter == 4:
            xyz = l.split()
            cells[curr_config][0][0] = float(xyz[0])
            cells[curr_config][0][1] = float(xyz[1])
            cells[curr_config][0][2] = float(xyz[2])
        
        elif iter == 5:
            xyz = l.split()
            cells[curr_config][1][0] = float(xyz[0])
            cells[curr_config][1][1] = float(xyz[1])
            cells[curr_config][1][2] = float(xyz[2])
        
        elif iter == 6:
            xyz = l.split()
            cells[curr_config][2][0] = float(xyz[0])
            cells[curr_config][2][1] = float(xyz[1])
            cells[curr_config][2][2] = float(xyz[2])
        
        elif iter >= 8 and iter < 8+num_A:
            xyz = l.split()

            ids[curr_config][iter-8] = int(xyz[0])
            types[curr_config][iter-8] = int(xyz[1])+1

            coords[curr_config][iter-8][0] = float(xyz[2])
            coords[curr_config][iter-8][1] = float(xyz[3])
            coords[curr_config][iter-8][2] = float(xyz[4])

            forces[curr_config][iter-8][0] = float(xyz[5])
            forces[curr_config][iter-8][1] = float(xyz[6])
            forces[curr_config][iter-8][2] = float(xyz[7])
        
        elif iter == 9+num_A:

            energies[curr_config] = float(l)
        iter += 1

for i in range(tot_CFGS):
    d2 = DFTD3(old=True, xc='b-lyp')
    d3 = DFTD3(old=False, xc='b-lyp')

    c = [cells[i][0][0], cells[i][1][1], cells[i][2][2]]
    atoms = Atoms(symbols=types[i], positions=coords[i], cell=c, pbc=[1,1,1])
    
    atoms.calc = d2
    es = atoms.get_potential_energy()
    fs = atoms.get_forces()[0][0]

    energies[i] -= es
    forces[i] -= fs

    atoms.calc = d3
    es = atoms.get_potential_energy()
    fs = atoms.get_forces()[0][0]

    energies[i] += es
    forces[i] += fs

    with open('converted.cfg', 'a') as file:
        file.write("BEGIN_CFG\n")
        file.write(" Size\n")
        file.write("    {}\n".format(num_A))
        file.write(" Supercell\n")
        file.write("    {}    0.0    0.0\n".format(c[0]))
        file.write("    0.0    {}    0.0\n".format(c[1]))
        file.write("    0.0    0.0    {}\n".format(c[2]))
        file.write(" AtomData:  id type    cartes_x    cartes_y    cartes_z    fx    fy    fz\n")
        for n in range(num_A):
            t = types[i][n]-1
            file.write("    {}    {}    {}    {}    {}    {}    {}    {}\n".format(ids[i][n], t, coords[i][n][0], coords[i][n][1], coords[i][n][2], 
                forces[i][n][0], forces[i][n][1], forces[i][n][2]))
        file.write(" Energy\n")
        file.write("    {}\n".format(energies[i]))
        file.write("END_CFG\n")
        file.write("\n")

# write_force(ids=ids, forces=forces, symbols=types, positions=coords, cells=cells, iter=i)