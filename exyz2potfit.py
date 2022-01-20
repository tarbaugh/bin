import sys
import numpy as np
import argparse
import random

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

lammps_dump_file = open(file, "r")

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
    if len(line.split() == 1) in line:
        line1 = 'ITEM: NUMBER'
        continue
    if 'ITEM: NUMBER' in line1:
        cfg_size = int(line)
        line1 = ''
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
        if idx == 0:
            arr = np.zeros((cfg_size, 9)).astype(np.float64)
        for x in line.split():
            arr[idx][idx1] = float(x)
            if (idx1 == 2): 
                pos_x.append(float(x))
            if (idx1 == 3): 
                pos_y.append(float(x))
            if (idx1 == 4): 
                pos_z.append(float(x))
            if (idx1 == 5): 
                f_x.append(float(x))
            if (idx1 == 6): 
                f_y.append(float(x))
            if (idx1 == 7): 
                f_z.append(float(x))
            if (idx1 == 8):
                epa.append(float(x))
            idx1 += 1
        idx += 1
        if (idx == cfg_size):
            if (valid_count % 13) == 0:
                arr = arr[arr[:, 0].argsort()]
                line1 = ''
                idx = 0
                mlip_configuration_efs = open("train.cfg", "a")
                for i in range (0,cfg_size):
                    energy += epa[i]
                mlip_configuration_efs.write("#N "+str(cfg_size)+" 1\n")
                mlip_configuration_efs.write("#C Ge Te\n")
                mlip_configuration_efs.write("#X "+str(lattice[1])+" "+str(lattice[2])+" "+str(lattice[4])+"\n")
                mlip_configuration_efs.write("#Y "+str(lattice[0])+" "+str(lattice[3])+" "+str(lattice[4])+"\n")
                mlip_configuration_efs.write("#Z "+str(lattice[0])+" "+str(lattice[2])+" "+str(lattice[5])+"\n")
                mlip_configuration_efs.write("#E "+str(energy)+"\n")
                mlip_configuration_efs.write("#F\n")
                for i in range (0,cfg_size):
                    mlip_configuration_efs.write(str(int(arr[i][1])-1)+" "+str(arr[i][2])+" "+str(arr[i][3])+" "+str(arr[i][4])+" "+str(arr[i][5])+" "+str(arr[i][6])+" "+str(arr[i][7])+"\n")
                mlip_configuration_efs.write("\n")

                mlip_configuration_efs.close()
                valid_count += 1
            else:
                line1 = ''
                idx = 0
                valid_count += 1
                continue

lammps_dump_file.close()
