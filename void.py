import numpy as np
from scipy.spatial import KDTree
import argparse

parser = argparse.ArgumentParser(description='LAMMPS .xyz MSD code')
parser.add_argument('f', nargs='+')
parser.add_argument('r', nargs=1)

def fileLen(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getAtomsTimes(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            if i == 0:
                return(int(l))
            # if i == 1:
            #     fTime = int(l.split()[-1])
            # if i != 1 and l.split()[0] == 'Atoms.':
            #     dTime = int(l.split()[-1]) - fTime


def build_Arr(filename, nA, npArr):
    with open(filename) as f:
        c = -1
        for i, l in enumerate(f):
            if i % (nA+2) == 0:
                c += 1
            elif i % (nA+2) == 1:
                pass
            else:
                txyz = l.split()
                for p in range(3):
                    npArr[c][(i%(nA+2))-2][p] = float(txyz[p+1])

def void(c, nAtoms, npArr, r):
    points = []
    dim = 100
    x_max = y_max = z_max = 0.0
    for a in range(nAtoms):
        if npArr[c][a][0] > x_max:
            x_max = npArr[c][a][0]
        if npArr[c][a][1] > y_max:
            y_max = npArr[c][a][1]
        if npArr[c][a][2] > z_max:
            z_max = npArr[c][a][2]

    tree = KDTree(npArr[c])
    zeds = np.zeros((dim, dim, dim))
    x_mult = x_max/dim
    y_mult = y_max/dim
    z_mult = z_max/dim
    for i in range(10, dim-10):
        for j in range(10, dim-10):
            for k in range(10, dim-10):
                q = tree.query_ball_point((i*x_mult, j*y_mult, k*z_mult), r, return_length=True)
                if q == 0:
                    zeds[i][j][k] = 1
                    points.append((i*x_mult, j*y_mult, k*z_mult))

    return 100*np.sum(zeds)/(dim*dim*dim), points




def main():
    args = parser.parse_args()
    fname = str(args.f[0])
    rad = float(args.r[0])
    print("Computing {} with radius {}".format(fname, rad))

    fLength = fileLen(filename=fname)
    nAtoms = getAtomsTimes(filename=fname)
    cycles = int(fLength/(nAtoms+2))

    trajs = np.zeros((cycles, nAtoms, 3))
    build_Arr(filename=fname, nA=nAtoms, npArr=trajs)

    # for c in range(cycles):
    #     vol = void(c, nAtoms, npArr=trajs)
    #     print(c, vol)

    vol, points = void(0, nAtoms, npArr=trajs, r=rad)
    print(vol)
    with open("1500KMSD copy.xyz", 'w') as f:
        for x in points:
            f.write("3 {} {} {}\n".format(x[0], x[1], x[2]))
    print(len(points))
        
if __name__ == "__main__":
    main()