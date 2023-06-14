from ase.io import read, write
from ase import Atoms
import numpy as np
import argparse

kB = 1.380649*10e-23

def unwrap(dump, unwrapped, write_xyz=False):
    atoms = dump[0]
    prev_positions = atoms.get_positions()
    tallies = np.zeros(prev_positions.shape)

    for i in range(1, len(dump)):
        atoms = dump[i]
        positions = atoms.get_positions()
        x_dim, y_dim, z_dim = atoms.get_cell()
        x_dim = x_dim[0]
        y_dim = y_dim[1]
        z_dim = z_dim[2]
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

def msd(r):
    msd = np.zeros(len(r))
    for i in range(len(r)):
        disp = r[i] - r[0]
        disp = np.square(disp)
        disp = disp.sum(axis=1)
        msd[i] = disp.mean()

    return msd

def msd_straight_forward(r):
    shifts = np.arange(len(r))
    msd = np.zeros(shifts.size)    

    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs).sum(axis=2)
        msd[i] = sqdist.mean()

    return msd

def main():
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
    parser.add_argument('f', nargs='+')
    args = parser.parse_args()
    files = args.f

    colors = ["red", "orange", "gold", "green", "blue"]
    viscs = []
    temps = []
    for i, file in enumerate(files):
        dump = read(file, index=":", format="lammps-dump-text")
        unwrapped_positions = np.array(list(map(lambda x: x.get_positions(), dump)))
        unwrap(dump, unwrapped_positions, write_xyz = False)
        fft = msd_fft(unwrapped_positions)
        xs = np.linspace(0, len(fft), num=len(fft))
        a, _ = np.polyfit(xs[len(xs//2):]*10e-12, fft[len(fft//2):], 1)
        D = 10e-10*a/6
        visc = kB*int(file[:3])/(6*np.pi*D)
        viscs.append(D)
        temps.append(int(file[:3]))
        # plt.plot(xs, fft, color = colors[i%len(colors)], label = file[:3]+"K; visc. = {}".format(visc))
        plt.plot(xs, fft, color=colors[i%len(colors)], label=file[:3]+"K; D = {}".format(D))
    plt.legend()
    plt.show()
    plt.plot(temps, viscs, marker='o')
    plt.show()

if __name__ == "__main__":
    main()