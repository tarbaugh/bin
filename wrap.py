from ase.io import read
import numpy as np
import os

c = 18.6078

a = read('GeTeMD-pos-1.xyz', index=':')
f = read('GeTeMD-frc-1.xyz', index=':')

if os.path.isfile('md.nb'):
    os.remove('md.nb')
for i, b in enumerate(a):
    b.set_cell(c * np.identity(3))
    b.set_pbc((1,1,1))
    b.wrap()
    with open('md.nb', 'a') as file:
        file.write("ITEM: TIMESTEP\n")
        file.write("{}\n".format(i))
        file.write("ITEM: NUMBER OF ATOMS\n")
        file.write("{}\n".format(len(b.get_atomic_numbers())))
        file.write("ITEM: BOX BOUNDS pp pp pp\n")
        file.write("{} {}\n".format(0.0, c))
        file.write("{} {}\n".format(0.0, c))
        file.write("{} {}\n".format(0.0, c))
        file.write("ITEM: ATOMS id type x y z fx fy fz\n")
        pos = b.get_positions()
        frcs = f[i].get_positions()
        chems = b.get_chemical_symbols()
        assert chems == a[i].get_chemical_symbols()
        for n in range(len(pos)):
            if chems[n] == 'Te':
                t = 2
            else:
                t = 1
            file.write("{} {} {} {} {} {} {} {}\n".format(n, t, pos[n][0], pos[n][1], pos[n][2], frcs[n][0]*51.42208619083232, frcs[n][1]*51.42208619083232, frcs[n][2]*51.42208619083232))