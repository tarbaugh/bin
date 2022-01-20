import sys
import numpy as np
import argparse
import random

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

read_file = file
write_file = "input.pos"

lammps_dump_file = open(read_file, "r")

line1 = ''
lattice = [0,0,0,0,0,0]
pos_x = []
pos_y = []
pos_z = []
f_x = []
f_y = []
f_z = []
epa = []
energy = 0
idx = 0

valid_count = 0

for line in lammps_dump_file:
    if 'ITEM: NUMBER' in line:
        line1 = 'ITEM: NUMBER'
        continue
    if 'ITEM: NUMBER' in line1:
        cfg_size = int(line)
        epa = []
        energy = 0
        line1 = ''
        arr = np.zeros((cfg_size, 9)).astype(np.float64)
        lattice = [0,0,0,0,0,0]
    if 'ITEM: BOX' in line:
        line1 = 'ITEM: BOX'
        continue
    if 'ITEM: BOX' in line1:
        idx1 = 0
        for x in line.split():
            lattice[2*idx+idx1] = float(x)
            idx1 += 1
        idx += 1
        if (idx == 3): 
            line1 = ''
            idx = 0
        continue
    if 'ITEM: ATOMS' in line:
        line1 = 'ITEM: ATOMS'
        continue
    if 'ITEM: ATOMS' in line1:
        idx1 = 0
        for x in line.split():
            arr[idx][idx1] = float(x)
            if (idx1 == 2): 
                if float(x) < lattice[0]:
                    new_x = float(x) + (lattice[1]-lattice[0])
                elif float(x) > lattice[1]:
                    new_x = float(x) - (lattice[1]-lattice[0])
                else:
                    new_x = float(x)
                # pos_x.append(new_x)
                arr[idx][idx1] = new_x
            if (idx1 == 3):
                if float(x) < lattice[2]:
                    new_y = float(x) + (lattice[3]-lattice[2])
                elif float(x) > lattice[3]:
                    new_y = float(x) - (lattice[3]-lattice[2])
                else:
                    new_y = float(x)
                # pos_y.append(new_y)
                arr[idx][idx1] = new_y
            if (idx1 == 4):
                if float(x) < lattice[4]:
                    new_z = float(x) + (lattice[5]-lattice[4])
                elif float(x) > lattice[5]:
                    new_z = float(x) - (lattice[5]-lattice[4])
                else:
                    new_z = float(x)
                # pos_z.append(new_z)
                arr[idx][idx1] = new_z
            # if (idx1 == 5): 
            #     f_x.append(float(x))
            # if (idx1 == 6): 
            #     f_y.append(float(x))
            # if (idx1 == 7): 
            #     f_z.append(float(x))
            if (idx1 == 8):
                epa.append(float(x))
            idx1 += 1
        idx += 1
        if (idx == cfg_size):
            if (valid_count % 1) == 0:
                arr = arr[arr[:, 0].argsort()]
                line1 = ''
                idx = 0
                mlip_configuration_efs = open(write_file, "w")
    
                mlip_configuration_efs.write("Input File\n")
                mlip_configuration_efs.write('\t{}\tatoms\n'.format(cfg_size))
                mlip_configuration_efs.write('\t3\tatom types\n\n')
                mlip_configuration_efs.write('{}\t{}\txlo xhi\n'.format(lattice[0], lattice[1]))
                mlip_configuration_efs.write('{}\t{}\tylo yhi\n'.format(lattice[2], lattice[3]))
                mlip_configuration_efs.write('{}\t{}\tzlo zhi\n\n'.format(lattice[4], lattice[5]))

                mlip_configuration_efs.write('Atoms\n\n')
                for i in range (0,cfg_size):
                    mlip_configuration_efs.write('{}\t{}\t{}\t{}\t{}\n'.format(int(arr[i][0]), int(arr[i][1]), arr[i][2], arr[i][3], arr[i][4]))

                mlip_configuration_efs.close()
                valid_count += 1
            else:
                line1 = ''
                idx = 0
                valid_count += 1
                continue

lammps_dump_file.close()
