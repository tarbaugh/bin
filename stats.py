import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('-t', type=int, nargs='?', default=0)
parser.add_argument('log')

args = parser.parse_args()

logfile = args.log
transientNum = args.t

def fileLen(filename=logfile):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getStart(s='Per MPI rank memory allocation', filename=logfile):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                return (i+1)
    return -1

def getEnd(s='Loop time of', filename=logfile):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                return (i)
    return -1

def main():
    filelength = fileLen()
    startline = (getStart() + transientNum)
    endline = getEnd()

    if startline == -1:
        return 'Not Log File'

    if endline == -1:
        return 'Not Log File'

    srows = [i for i in range(0,startline)]+[i for i in range(endline,filelength)]

    df = pd.read_table(logfile, delim_whitespace=True, skiprows=srows)

    print(df.describe())

if __name__ == "__main__":
    main()

