import pandas as pd
import argparse
import os

def fileLen(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getStart(s, filename):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                return (i+1)
    return -1

def getEnd(s, filename):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                return (i)
    return -1

def getFiles(ending):
    fs = []
    for _,dirs,_ in os.walk(os.getcwd()):
        for dir in dirs:
            for root,_,files in os.walk(dir):
                for file in files:
                    if file.endswith(ending):
                        fs.append(os.path.join(root, file))
    return fs

def main():
    parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
    parser.add_argument('-t', type=int, nargs='?', default=0)
    parser.add_argument('endings')
    parser.add_argument('x')
    parser.add_argument('y')

    args = parser.parse_args()

    ends = args.endings
    x_label = args.x
    y_label = args.y
    transientNum = args.t


    files = getFiles(ending=ends)
    avgs = []
    for f in files:
        filelength = fileLen(filename=f)
        startline = (getStart(s='Per MPI rank memory allocation', filename=f) + transientNum)
        endline = getEnd(s='Loop time of', filename=f)

        if startline == -1:
            return 'Not Log File'

        if endline == -1:
            return 'Not Log File'

        srows = [i for i in range(0,startline)]+[i for i in range(endline,filelength)]

        df = pd.read_table(f, delim_whitespace=True, skiprows=srows)

        avgs.append((df.loc[:,x_label].mean(), df.loc[:,y_label].mean()))

    avgs.sort(key=lambda x: x[0])

    for x in avgs:
        s = ""
        for y in x:
            s += str(y)+"\t"
        print(s)

if __name__ == "__main__":
    main()

