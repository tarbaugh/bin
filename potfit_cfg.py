import sys
import numpy as np

file = "train_GeTe.xyz"

def w(cfg_size, lattice, energy, posf, types):
    with open('potfit.cfg', "a") as mlip_configuration_efs:    
        mlip_configuration_efs.write("#N "+str(cfg_size)+" 1\n")
        mlip_configuration_efs.write("#C Ge Te\n")
        mlip_configuration_efs.write("#X "+str(lattice[0])+" "+str(lattice[3])+" "+str(lattice[6])+"\n")
        mlip_configuration_efs.write("#Y "+str(lattice[1])+" "+str(lattice[4])+" "+str(lattice[7])+"\n")
        mlip_configuration_efs.write("#Z "+str(lattice[2])+" "+str(lattice[5])+" "+str(lattice[8])+"\n")
        mlip_configuration_efs.write("#E "+str(energy)+"\n")
        mlip_configuration_efs.write("#F\n")
        for i in range (0,cfg_size):
            mlip_configuration_efs.write(str(types[i])+"    "+str(posf[i][0])+"    "+str(posf[i][1])+"    "+str(posf[i][2])+"    "+str(posf[i][3])+"    "+str(posf[i][4])+"    "+str(posf[i][5])+"\n")
        mlip_configuration_efs.write("\n")

with open(file, 'r') as f:
    for i, l in enumerate(f):
        if len(l.split()) == 1:
            if i != 0:
                w(cfg_size=numAtoms, lattice=lat, energy=nrgy, posf=xyzf, types=t)
            numAtoms = int(l.split()[0])
            xyzf = np.zeros((numAtoms, 6))
            count = 0
            t = [-1]*numAtoms
        elif len(l.split()) > 15:
            params = l.split()
            lat = [0.0]*9
            for iter in range(9):
                lat[iter] = params[iter+1]
            lat[-1] = lat[-1][:-1]
            nrgy = params[21][7:]
        else:
            l = l.split()
            if str(l[0]) == 'Ge':
                t[count] = 0
            if str(l[0]) == 'Te':
                t[count] = 1
            xyzf[count][0] = float(l[1])
            xyzf[count][1] = float(l[2])
            xyzf[count][2] = float(l[3])
            xyzf[count][3] = float(l[4])
            xyzf[count][4] = float(l[5])
            xyzf[count][5] = float(l[6])
            count += 1
