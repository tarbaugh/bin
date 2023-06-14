from ase.io import read, write
from ase import Atoms
import numpy as np
import argparse
import cmath

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
                    tot_cells.append([x_low, x_hi, y_low, y_hi, z_low, z_hi])
                    tot_types.append(npatoms[:,2])
                atoms = []
            elif l.startswith("ITEM: NUMBER OF ATOMS"):
                numA_i = i + 1
            elif l.startswith("ITEM: BOX BOUNDS pp pp pp"):
                cells_i = i + 1
            elif l.startswith("ITEM: ATOMS id type x y z fx fy fz"):
                atoms_i = i + 1
            elif i == timestep_i:
                timestep = int(l)
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

    return np.array(tot_pos), tot_cells

def unwrap(wrapped, unwrapped, cells):
    '''Unwraps coordinates, optionally storing to new dump
    Input: [np.ndarray, np.ndarray, [optional: float], [optional: bool]]
        (# Configurations, # Atoms, 3), (# Configurations, # Atoms, 3), [optional: cubic box length], [optional: whether to write to knew dump]
    Output: np.ndarray (# Configurations, # Atoms, 3)
    '''
    prev_positions = wrapped[0]
    tallies = np.zeros(prev_positions.shape)

    for i in range(1, len(wrapped)):
        positions = wrapped[i]
        cell = cells[i]
        disp = positions - prev_positions

        x_dim = abs(cell[1] - cell[0])
        y_dim = abs(cell[3] - cell[2])
        z_dim = abs(cell[5] - cell[4])

        for iiter in range(disp.shape[0]):
            for jiter in range(disp.shape[1]):
                if disp[iiter][jiter] >= x_dim/2:
                    tallies[iiter][jiter] -= x_dim
                if disp[iiter][jiter] <= -x_dim/2:
                    tallies[iiter][jiter] += x_dim
        
        unwrapped[i] +=  tallies
        
        prev_positions = positions.copy()

def unwrap_ase(dump, unwrapped, cell=None, write_xyz=False):
    '''Unwraps coordinates, optionally storing to new dump
    Input: [np.ndarray, np.ndarray, [optional: float], [optional: bool]]
        (# Configurations, # Atoms, 3), (# Configurations, # Atoms, 3), [optional: cubic box length], [optional: whether to write to knew dump]
    Output: np.ndarray (# Configurations, # Atoms, 3)
    '''
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
        write('unwrapped.xyz', new_atoms, format = 'xyz')

def autocorrelation(x):
    '''Calculates autocorrelation of cartesian coordinates in one plane
    Input: np.ndarray (# Configurations, # Atoms, 1)
    Output: np.ndarray ()
    '''
    N = x.shape[0]
    F = np.fft.fft(x, n=2*N, axis=0)
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD, axis=0)
    res = (res[:N]).real
    n = np.arange(1, N+1)[::-1]

    return res/n[:, np.newaxis]

def msd_fft(r):
    '''Calculates MSD of trajectory; credit to Freud analysis module and subsequent contributors
    Input: np.ndarray (# Configurations, # Atoms, 3)
    Output: np.ndarray (# Configurations)
    '''
    N = r.shape[0]
    D = np.square(r).sum(axis=2) 
    D = np.append(D,np.zeros(r.shape[:2]), axis=0)
    Q = 2*D.sum(axis=0)
    S1 = np.zeros(r.shape[:2])
    for m in range(N):
        Q -= (D[m-1, :] + D[N-m, :])
        S1[m, :] = Q/(N-m)

    corrs = []
    for i in range(r.shape[2]):
        corrs.append(autocorrelation(r[:, :, i]))
    S2 = np.sum(corrs, axis=0)

    return (S1 - 2*S2).mean(axis=-1)

def fs(r, k):
    t = []
    c = 0
    dirs = np.array([[k, 0, 0], [0, k, 0], [0, 0, k]])
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

def fs_msd(r, k):
    '''Calculates intermediate scattering function Fs using FFT MSD
    Input: [np.ndarray, int] [(# Configurations, # Atoms, 1), Wave vector magnitude]
    Output: np.ndarray ()
    '''
    msd = msd_fft(r)
    return msd, np.array([cmath.exp((-1/6)*k*k*m) for m in msd])

def validate_args(args):
    if args.t != None:
        types = list(map(int, args.t))
    else:
        types = None
    if args.CALC_FS:
        assert args.k != None
        print('K chosen as '+str(args.k[0]))
        return float(args.k[0]), types
    return None, types

def main():
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
    parser.add_argument('--msd', dest='CALC_MSD', default=False, action='store_true')
    parser.add_argument('--fs', dest='CALC_FS', default=False, action='store_true')
    parser.add_argument('--avg', dest='AVG', default=False, action='store_true')
    parser.add_argument('-f', nargs='+')
    parser.add_argument('-t', nargs='+')
    parser.add_argument('-k', nargs=1, default=None)
    parser.add_argument('--format', dest='ASE_FORMAT', default='lammps-dump-text')
    args = parser.parse_args()
    k, types = validate_args(args)
    files = args.f

    # types = dump[0].get_atomic_numbers()
    # types = list(set(types))
    # use_types = []
    # for t in types:
    #     inp = input("Use type {}?\n".format(t))
    #     if inp.lower() != "n" or inp.lower() != "no":
    #         use_types.append(t)
    # types = dump[0].get_atomic_numbers()

    # del_rows = []
    # for i, t in enumerate(types):
    #     if t not in use_types:
    #         del_rows.append(t)

    colors = ['red', 'orange', 'gold', 'green', 'blue']
    if args.AVG:
        tot_xs = []
        tot_ys = []

    curr_temps = []
    for file in files:
        if file.split('.')[0] not in curr_temps:
            curr_temps.append(file.split('.')[0])

    if args.AVG:
        for ct in curr_temps:
            tot_xs = []
            tot_ys = []
            for file in files:
                dump = read(file, index=':', format=args.ASE_FORMAT)
                unwrapped_positions = np.array(list(map(lambda x: x.get_positions(), dump)))
                unwrap(dump, unwrapped_positions, write_xyz=False)
                if file.split('.')[0] != ct:
                    continue
                if args.CALC_MSD:
                    ys = msd_fft(unwrapped_positions)
                    xs = np.linspace(0, len(ys), num=len(ys))
                    if args.AVG:
                        tot_xs.append(xs)
                        tot_ys.append(ys)

            xs = np.zeros(shape=tot_xs[0].shape)
            for x in tot_xs:
                xs += x
            xs /= len(tot_xs)
            ys = np.zeros(shape=tot_ys[0].shape)
            for y in tot_ys:
                ys += y
            ys /= len(tot_ys)
            plt.plot(xs, ys, color=colors[i%len(colors)], label=file.split('.')[0]+'K')

    else:
        for i, file in enumerate(files):
            wrapped_pos, cells = pyread(file, types)
            unwrapped_pos = wrapped_pos.copy()
            # unwrapped_positions = np.array(list(map(lambda x: x.get_positions(), dump)))
            unwrap(wrapped_pos, unwrapped_pos, cells)
            if args.CALC_MSD and args.CALC_FS:
                msd_ys, fs_ys = fs_msd(unwrapped_pos, k)
                xs = np.linspace(0, len(msd_ys), num=len(msd_ys))
                
                plt.plot(xs, msd_ys, color=colors[i%len(colors)], label=file.split('.')[0]+'K')
                # plt.legend()
                # plt.show()
                plt.plot(xs, fs_ys, color=colors[i%len(colors)], label=file.split('.')[0]+'K')
                # plt.legend()
                # plt.show()
            elif args.CALC_MSD:
                ys = msd_fft(unwrapped_pos)
                # np.savetxt(file+"msd", ys)
                xs = np.linspace(0, len(ys), num=len(ys))
                plt.plot(xs, ys, color=colors[i%len(colors)], label=file.split('.')[0]+'K')
                # plt.legend()
                # plt.show()
            elif args.CALC_FS:
                # _, ys = fs_msd(unwrapped_positions, k)
                ys = fs(unwrapped_pos, k)
                xs = np.linspace(0, len(ys), num=len(ys))
                plt.plot(xs, ys, color=colors[i%len(colors)], label=file.split('.')[0]+'K')
                # plt.legend()
                # plt.show()
        plt.legend()
        plt.show()
if __name__ == '__main__':
    main()