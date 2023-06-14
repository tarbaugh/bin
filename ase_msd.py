from ase.io import read
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('f')
args = parser.parse_args()
file = args.f

dump = read(file, index=":", format="lammps-dump-text")

windows = []

for j in range(1, len(dump)):

    atoms = dump[0]
    positions = atoms.get_positions()
    x_dim, y_dim, z_dim = atoms.get_cell()
    x_dim = x_dim[0]
    y_dim = y_dim[1]
    z_dim = z_dim[2]
    assert float(x_dim) == float(y_dim)
    assert float(y_dim) == float(z_dim)
    prev_positions = positions.copy()

    count = 0
    tmp = []

    for i in range(j, len(dump), j):
        atoms = dump[i]
        positions = atoms.get_positions()
        x_dim, y_dim, z_dim = atoms.get_cell()
        x_dim = x_dim[0]
        y_dim = y_dim[1]
        z_dim = z_dim[2]
        assert float(x_dim) == float(y_dim)
        assert float(y_dim) == float(z_dim)

        disp = positions - prev_positions

        for iiter in range(len(disp)):
            for jiter in range(3):
                if disp[iiter][jiter] >= x_dim/2:
                    disp[iiter][jiter] = disp[iiter][jiter] - x_dim
                if disp[iiter][jiter] <= -x_dim/2:
                    disp[iiter][jiter] = disp[iiter][jiter] + x_dim
        
        disp = disp*disp
        av = np.average(np.sqrt(np.sum(disp, axis=1)))
        tmp.append(av)
        count += 1

        prev_positions = positions.copy()
    
    windows.append(sum(tmp)/count)

print(windows)
import matplotlib.pyplot as plt
plt.plot(windows)
plt.show()