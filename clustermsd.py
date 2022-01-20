import numpy as np
import argparse

parser = argparse.ArgumentParser(description='LAMMPS .xyz MSD code')
parser.add_argument('f', nargs='+')

def fileLen(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getAtomsTimes(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            if i == 0:
                nA = int(l)
            if i == 1:
                fTime = int(l.split()[-1])
            if i != 1 and l.split()[0] == 'Atoms.':
                dTime = int(l.split()[-1]) - fTime
                return nA, fTime, dTime

def getArrays(filename, length, nA, npArr):
    s_header = 2
    # l = []
    i = 0
    while (s_header + nA) <= length+2:
        df = np.genfromtxt(filename, skip_header=s_header, max_rows=nA)
        iter = 0
        for x in df:
            npArr[iter][i] = x[1:]
            iter += 1
        s_header += (nA+2)
        i += 1

def autocorrFFT(x):
  N=len(x)
  F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
  PSD = F * F.conjugate()
  res = np.fft.ifft(PSD)
  res= (res[:N]).real   #now we have the autocorrelation in convention B
  n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
  return res/n #this is the autocorrelation in convention A

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

def build_Arr(filename, nA, npArr, npDelts, boxS):
    with open(filename) as f:
        c = -1
        for i, l in enumerate(f):
            if i % (nA+2) == 0:
                c += 1
            elif i % (nA+2) == 1:
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



def main():  
    # fname = "1000K100psNVT_3.xyz"
    args = parser.parse_args()
    fname = str(args.f[0])

    boxSize = float(input("Box size: \n"))

    fLength = fileLen(filename=fname)
    nAtoms, _, deltaTime = getAtomsTimes(filename=fname)
    cycles = int(fLength/(nAtoms+2))

    # Revert to old array shape -> use MDAnalysis
    trajs = np.zeros((nAtoms, cycles, 3))
    deltas = np.zeros((nAtoms, 3))
    # getArrays(filename=fname, length=fLength, nA=nAtoms, npArr=trajs)
    build_Arr(filename=fname, nA=nAtoms, npArr=trajs, npDelts=deltas, boxS=boxSize)
    # print(trajs)

    result = np.zeros_like(msd_fft(trajs[0]))
    for xyz in trajs:
        result = np.add(result, msd_fft(xyz))
    result /= nAtoms
    np.savetxt(fname+'.msddat', result)
    

    

if __name__ == "__main__":
    main()