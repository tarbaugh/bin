from ase.io import read, write
from ase import Atoms
import numpy as np
import argparse
import cmath
import math

def unwrap(dump, unwrapped, cell=None, write_xyz=False):
    atoms = dump[0]
    prev_positions = atoms.get_positions()
    tallies = np.zeros(prev_positions.shape)

    for i in range(1, len(dump)):
        atoms = dump[i]
        positions = atoms.get_positions()

        if cell != None:
            assert type(cell) == type(1.0) or type(cell) == type(1)
            x_dim = float(cell)

        else:
            x_dim, y_dim, z_dim = atoms.get_cell()
            x_dim = x_dim[0]
            y_dim = y_dim[1]
            z_dim = z_dim[2]
            ## Assert square cell and isotropic expansion
            assert float(x_dim) == float(y_dim)
            assert float(y_dim) == float(z_dim)

        disp = positions - prev_positions

        for iiter in range(disp.shape[0]):
            for jiter in range(disp.shape[1]):
                if disp[iiter][jiter] >= x_dim/2:
                    tallies[iiter][jiter] -= x_dim
                if disp[iiter][jiter] <= -x_dim/2:
                    tallies[iiter][jiter] += x_dim
        
        unwrapped[i] +=  tallies
        
        prev_positions = positions.copy()
    
    if write:
        new_atoms = []
        for p in unwrapped:
            new_atoms.append(Atoms(positions=p))
        write("unwrapped.xyz", new_atoms, format = "xyz")

def fs(r, dirs):
    t = []
    c = 0
    while c < len(r):
        sum = 0
        j = 0
        av = 0
        while j < len(r[0]):
            disp = r[c][j]-r[0][j]
            for dir in dirs:
                dot = np.dot(dir, disp)
                sum += cmath.exp(complex(0,dot))
                av += 1
            j += 1
        
        sum /= (j+av)
        t.append(sum)
        c += 1
    return t

def autocorrelation(x):
    N = x.shape[0]
    F = np.fft.fft(x, n=2*N, axis=0)
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD, axis=0)
    res = (res[:N]).real
    n = np.arange(1, N+1)[::-1]  # N to 1

    return res/n[:, np.newaxis]

def msd_fft(r):
    N = r.shape[0]
    D = np.square(r).sum(axis=2) 
    D = np.append(D,np.zeros(r.shape[:2]), axis=0)
    Q = 2*D.sum(axis=0)
    S1 = np.zeros(r.shape[:2])
    for m in range(N):
        Q -= (D[m-1, :] + D[N-m, :])
        S1[m, :] = Q/(N-m)

    # The second term can be computed via autocorrelation
    corrs = []
    for i in range(r.shape[2]):
        corrs.append(autocorrelation(r[:, :, i]))
    S2 = np.sum(corrs, axis=0)

    return (S1 - 2*S2).mean(axis=-1)

def fs_msd(r, dirs):
    msd = msd_fft(r)
    return np.array([cmath.exp((-1/6)*dirs[0][0]*dirs[0][0]*m) for m in msd])

def main():
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
    parser.add_argument('-k', nargs=1)
    parser.add_argument('-f', nargs='+')
    parser.add_argument('--format', dest='ASE_FORMAT', default="lammps-dump-text")
    args = parser.parse_args()
    ASE_FORMAT = args.ASE_FORMAT
    k = float(args.k[0])
    print("K chosen as {}".format(k))
    dirs = np.array([[k, 0, 0], [0, k, 0], [0, 0, k], [k/math.sqrt(2), k/math.sqrt(2), 0], \
        [k/math.sqrt(2), 0, k/math.sqrt(2)], [0, k/math.sqrt(2), k/math.sqrt(2)], [k/math.sqrt(3), \
        k/math.sqrt(3), k/math.sqrt(3)], [-k, 0, 0], [0, -k, 0], [0, 0, -k], [-k/math.sqrt(2), -k/math.sqrt(2), 0], \
        [-k/math.sqrt(2), 0, -k/math.sqrt(2)], [0, -k/math.sqrt(2), -k/math.sqrt(2)], [-k/math.sqrt(3), -k/math.sqrt(3), -k/math.sqrt(3)]])
    files = args.f

    colors = ["red", "orange", "gold", "green", "blue"]
    for i, file in enumerate(files):
        dump = read(file, index=":", format=ASE_FORMAT)
        unwrapped_positions = np.array(list(map(lambda x: x.get_positions(), dump)))
        unwrap(dump, unwrapped_positions, write_xyz = False)
        t = fs_msd(unwrapped_positions, dirs)
        xs = np.linspace(0, len(t), num=len(t))
        plt.plot(xs, t, color=colors[i%len(colors)], label=file[:3]+"K")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()