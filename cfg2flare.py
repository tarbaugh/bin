from multiprocessing import cpu_count
from turtle import position
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

write_file = "out.train"

tot_pos = []
tot_forces = []
tot_ids = []
tot_types = []
tot_cells = []
add = 0
 
with open(file, 'r') as f:
    for i, l in enumerate(f):
        if i == 2:
            num_A = int(l)
            break

with open(file, 'r') as f:
    i = 0
    for _, l in enumerate(f):
        if i == 2:
            num_A = int(l)
            cell = np.zeros((3,3)).astype(np.float64)
            coords = np.zeros((num_A, 3)).astype(np.float64)
            forces = np.zeros((num_A, 3)).astype(np.float64)
            ids = np.zeros((num_A,1)).astype(int)
            types = np.zeros((num_A,1)).astype(int)
        if i % (num_A + 12) == 4:
            xyz = l.split()
            cell[0][0] = float(xyz[0])
            cell[0][1] = float(xyz[1])
            cell[0][2] = float(xyz[2])
        
        if i % (num_A + 12) == 5:
            xyz = l.split()
            cell[1][0] = float(xyz[0])
            cell[1][1] = float(xyz[1])
            cell[1][2] = float(xyz[2])
        
        if i % (num_A + 12) == 6:
            xyz = l.split()
            cell[2][0] = float(xyz[0])
            cell[2][1] = float(xyz[1])
            cell[2][2] = float(xyz[2])
        
        if i % (num_A + 12) >= 8 and i % (num_A + 12) < 8+num_A:
            xyz = l.split()

            ids[i%(num_A + 12)-8][0] = int(xyz[0])
            types[i%(num_A + 12)-8][0] = int(xyz[1])+1

            coords[i%(num_A + 12)-8][0] = float(xyz[2])
            coords[i%(num_A + 12)-8][1] = float(xyz[3])
            coords[i%(num_A + 12)-8][2] = float(xyz[4])

            forces[i%(num_A + 12)-8][0] = float(xyz[5])
            forces[i%(num_A + 12)-8][1] = float(xyz[6])
            forces[i%(num_A + 12)-8][2] = float(xyz[7])
        
        if i == 10+num_A:
            if 'Feature' in l.split():
                add = 1
            else:
                add = 0
                

        if i == add+12+num_A:
            add = 0
            i = 0
            tot_pos.append(coords)
            tot_forces.append(forces)
            tot_ids.append(ids)
            tot_types.append(types)
            tot_cells.append(cell)
        
        i += 1
    
    # tot_pos.append(coords)
    # tot_forces.append(forces)
    # tot_ids.append(ids)
    # tot_types.append(types)
    # tot_cells.append(cell)

from flare.gp import GaussianProcess
from flare.struc import Structure

strucs = []
num_strucs = len(tot_pos)
i = 0

kernels = ['twobody', 'threebody']
component = 'mc'
hyps = np.array([0.1, 0.1, 0.1, 2.0, 0.5])  # initial (bad) choice of hyps
cutoffs = {'twobody': 6.0, 'threebody': 3.5}  # cutoff radii in A
maxiter = 100  # max number of hyperparameter optimziation steps

gp_model = GaussianProcess(
    kernels=kernels,
    component=component,
    hyps=hyps,
    cutoffs=cutoffs,
    maxiter=50
)

for snapshot in range(len(tot_pos)//2):
    # create flare structure
    training_positions = tot_pos[snapshot]
    training_forces = tot_forces[snapshot]
    species = np.reshape(tot_types[i], len(tot_types[i]))
    training_structure = Structure(cell, species, training_positions)
    # strucs.append(Structure(cell=tot_cells[i], 
    #     species=['Ge', 'Te'], 
    #     positions=tot_pos[i], 
    #     mass_dict={1:72.64, 2:127.603}, 
    #     species_labels=np.reshape(tot_types[i], len(tot_types[i]))))

    # add the structure to the training set of the GP
    gp_model.update_db(training_structure, np.array(tot_forces[snapshot]))

gp_model.set_L_alpha()

