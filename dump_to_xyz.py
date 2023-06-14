from ase import Atoms
from ase.io import write
import numpy as np
import glob, os

def pyread(file, types):
    '''Unwraps coordinates, optionally storing to new dump
    Input: [np.ndarray, np.ndarray, [optional: float], [optional: bool]]
        (# Configurations, # Atoms, 3), (# Configurations, # Atoms, 3), [optional: cubic box length], [optional: whether to write to knew dump]
    Output: np.ndarray (# Configurations, # Atoms, 3)
    '''
    with open(file) as f:
        atoms = None
        tot_pos = []
        tot_cells = []
        tot_types = []
        tot_times = []
        for i, l in enumerate(f):
            if l.startswith("ITEM: TIMESTEP"):
                timestep_i = i + 1
                if atoms != None:
                    if types != None:
                        atoms_iter = 0
                        while atoms_iter < len(atoms):
                            if atoms[atoms_iter][1] not in types:
                                atoms.pop(atoms_iter)
                            else:
                                atoms_iter += 1
                    atoms = sorted(atoms, key=lambda x: x[0])
                    npatoms = np.array(atoms)
                    tot_pos.append(npatoms[:,2:5])
                    tot_cells.append([x_hi, y_low, z_low, x_low, y_hi, z_low, x_low, y_low, z_hi])
                    tot_types.append(npatoms[:,1])
                atoms = []
            elif l.startswith("ITEM: NUMBER OF ATOMS"):
                numA_i = i + 1
            elif l.startswith("ITEM: BOX BOUNDS pp pp pp"):
                cells_i = i + 1
            elif l.startswith("ITEM: ATOMS id type x y z fx fy fz"):
                atoms_i = i + 1
            elif i == timestep_i:
                tot_times.append(int(l))
            elif i == numA_i:
                numA = int(l)
            elif i == cells_i:
                tmp = list(map(float, l.split()))
                x_low, x_hi = tmp[0], tmp[1]
            elif i == cells_i+1:
                tmp = list(map(float, l.split()))
                y_low, y_hi = tmp[0], tmp[1]
            elif i == cells_i+2:
                tmp = list(map(float, l.split()))
                z_low, z_hi = tmp[0], tmp[1]
            elif i >= atoms_i and i < atoms_i+numA:
                lis = l.split()
                new_lis = []
                for i, v in enumerate(lis):
                    if i <= 1:
                        new_lis.append(int(v))
                    else:
                        new_lis.append(float(v))
                atoms.append(new_lis)

    return np.array(tot_pos), tot_cells, tot_types, tot_times

for f in glob.glob("*.xyz"):
    wrapped_pos, cells, types, times = pyread(f, None)
    
    # ASE Solution
    # tot_cfgs = []
    # for i in range(len(wrapped_pos)):
    #     wrapped_pos[i][:,0] += cells[i][3]
    #     wrapped_pos[i][:,1] += cells[i][1]
    #     wrapped_pos[i][:,2] += cells[i][2]
    #     tot_cfgs.append(Atoms(positions=wrapped_pos[i], cell=[cells[i][0]+cells[i][3], cells[i][4]+cells[i][1], cells[i][8]+cells[i][2]], numbers=types[i]))
    # write(filename=f+".xyz", images=tot_cfgs, format="xyz")

    with open(f+".starr", mode="w") as fout:
        for i in range(len(wrapped_pos)):
            wrapped_pos[i][:,0] -= cells[i][3]
            wrapped_pos[i][:,1] -= cells[i][1]
            wrapped_pos[i][:,2] -= cells[i][2]
            fout.write("{}\n".format(len(wrapped_pos[i][:,0])))
            fout.write("Atoms. Timestep: {}\n".format(int(times[i])))
            for j in range(len(wrapped_pos[i][:,0])):
                fout.write("{} {} {} {}\n".format(int(types[i][j]), wrapped_pos[i][j][0], wrapped_pos[i][j][1], wrapped_pos[i][j][2]))
