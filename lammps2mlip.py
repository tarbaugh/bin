import sys
import numpy as np
import argparse
import random

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

read_file = file
write_file = "train.cfg"

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
            if (valid_count > 10) and (valid_count % 5) == 0:
                arr = arr[arr[:, 0].argsort()]
                line1 = ''
                idx = 0
                if ((valid_count) % 1) == -1:
                    mlip_configuration_efs = open("valid.cfg", "a")
                else:
                    mlip_configuration_efs = open(write_file, "a")
    
                mlip_configuration_efs.write("BEGIN_CFG\n")
                mlip_configuration_efs.write(" Size\n")
                mlip_configuration_efs.write("    "+str(cfg_size)+"\n")
                mlip_configuration_efs.write(" Supercell\n")
                mlip_configuration_efs.write("    "+str(lattice[1])+"    "+str(lattice[2])+"    "+str(lattice[4])+"\n")
                mlip_configuration_efs.write("    "+str(lattice[0])+"    "+str(lattice[3])+"    "+str(lattice[4])+"\n")
                mlip_configuration_efs.write("    "+str(lattice[0])+"    "+str(lattice[2])+"    "+str(lattice[5])+"\n")
                mlip_configuration_efs.write(" AtomData:  id type    cartes_x    cartes_y    cartes_z    fx    fy    fz\n")
                for i in range (0,cfg_size):
                    mlip_configuration_efs.write("    "+str(int(arr[i][0]))+"    "+str(int(arr[i][1]-1))+"    "+str(arr[i][2])+"    "+str(arr[i][3])+"    "+str(arr[i][4])+"    "+str(arr[i][5])+"    "+str(arr[i][6])+"    "+str(arr[i][7])+"\n")
                    energy += epa[i]
                mlip_configuration_efs.write(" Energy\n")
                mlip_configuration_efs.write("    "+str(energy)+"\n")
                mlip_configuration_efs.write("END_CFG\n")
                mlip_configuration_efs.write("\n")

                mlip_configuration_efs.close()
                valid_count += 1
            else:
                line1 = ''
                idx = 0
                valid_count += 1
                continue

lammps_dump_file.close()
