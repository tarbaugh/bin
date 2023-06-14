import numpy as np
import os
import argparse

def fileLen(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getAtomsTimes(filename):
    bool2 = False
    with open(filename) as f:
        for i, l in enumerate(f):
            if i == 3:
                nA = int(l)
            if i == 1:
                fTime = int(l)
            if bool2 == True:
                dTime = int(l.split()[-1]) - fTime
                return nA, fTime, dTime 
            if i != 0 and len(l.split()) > 1:
                if l.split()[1] == 'TIMESTEP':
                    bool2 = True

def autocorrFFT(x):
  N = len(x)
  F = np.fft.fft(x, n=2*N)
  PSD = F * F.conjugate()
  res = np.fft.ifft(PSD)
  res = (res[:N]).real
  n = N*np.ones(N)-np.arange(0,N)
  return res/n 

def msd_fft(r):
  N=len(r)
  D=np.square(r).sum(axis=1) 
  D=np.append(D,0) 
  S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
  Q=2*D.sum()
  S1=np.zeros(N)
  for m in range(N):
      Q=Q-D[m-1]-D[N-m]
      S1[m]=Q/(N-m)
  return S1-2*S2

def build_Arr(filename, nA, npArr, npDelts):
    minmax = np.zeros(6)
    with open(filename) as f:
        c = -1
        for i, l in enumerate(f):
            if i % (nA+9) == 0:
                c += 1
            elif i % (nA+9) < 9:
                pass
            else:
                txyz = l.split()
                if c > 0:
                    for p in range(3):
                        tmp = float(txyz[p+1])
                        if tmp < minmax[2*p]:
                            minmax[2*p] = tmp
                        if tmp > minmax[2*p+1]:
                            minmax[2*p+1] = tmp
    
    xdiff = minmax[1] - minmax[0]
    ydiff = minmax[3] - minmax[2]
    zdiff = minmax[5] - minmax[4]
    boxS = max([xdiff, ydiff, zdiff]) 
    
    # Boxsize is determined as the cube of largest xyz length - this will lead to larger diffusions in uneven boxes
    # This could be corrected if which side of the box crossed is detected

    with open(filename) as f:
        c = -1
        for i, l in enumerate(f):
            if i % (nA+9) == 0:
                c += 1
            elif i % (nA+9) < 9:
                pass
            else:
                txyz = l.split()
                if c > 0:
                    for p in range(3):
                        dr = (float(txyz[p+1]) + npDelts[(i%(nA+2))-2][p]) - npArr[(i%(nA+2))-2][c-1][p]
                        if dr >= boxS/2:
                            npDelts[(i%(nA+2))-2][p] -= boxS
                        if dr <= -boxS/2:
                            npDelts[(i%(nA+2))-2][p] += boxS
                for p in range(3):
                    npArr[(i%(nA+2))-2][c][p] = float(txyz[p+1]) + npDelts[(i%(nA+2))-2][p]

def getMSD(fname):
    fLength = fileLen(filename=fname)
    nAtoms, _, deltaTime = getAtomsTimes(filename=fname)
    cycles = int(fLength/(nAtoms+2))

    trajs = np.zeros((nAtoms, cycles, 3))
    deltas = np.zeros((nAtoms, 3))
    build_Arr(filename=fname, nA=nAtoms, npArr=trajs, npDelts=deltas)

    result = np.zeros_like(msd_fft(trajs[0]))
    for xyz in trajs:
        result = np.add(result, msd_fft(xyz))
    result /= nAtoms
    np.savetxt(fname+'.msddat', result)

def main():
    path = os.getcwd()
    # for _,dirs,_ in os.walk(os.getcwd()):
    #     for dir in dirs:
    #         for root,_,files in os.walk(dir):
    #             for file in files:
    #                 if file.endswith(".xyz"):
    #                     f = os.path.join(path, root)
    #                     f = os.path.join(f, file)
                        # getMSD(fname=f)
    parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
    parser.add_argument('f')
    args = parser.parse_args()
    file = args.f

    getMSD(fname=file)
     

if __name__ == "__main__":
    main()