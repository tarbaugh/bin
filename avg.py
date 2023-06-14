import pandas as pd
import argparse
import os
import sys

def fileLen(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getStart(s, filename):
    starts = []
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                starts.append(i+1)
    return starts

def getEnd(s, filename):
    ends = []
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                ends.append(i)
    return ends

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
    parser.add_argument('-f', type=str, nargs='?', default="")
    parser.add_argument('-e', type=str, nargs='?', default="")
    parser.add_argument('--labels', action="extend", nargs='+', type=str)

    args = parser.parse_args()

    ends = args.e
    data_file = args.f
    labels = args.labels
    transientNum = args.t

    if ends == "" and data_file == "":
        sys.exit("Must pass a file through `-f file` or list of endings in subdirectories through `-e ending`")
    
    if labels == "":
        sys.exit("Must pass a list of labels through `--labels label_1 label_2`")

    avgs = []

    if ends != "":
        files = getFiles(ending=ends)
        for f in files:
            filelength = fileLen(filename=f)
            tmpS = getStart(s='Per MPI rank memory allocation', filename=f)
            startlines = [x+transientNum for x in tmpS]

            endlines = getEnd(s='Loop time of', filename=f)

            if len(startlines) == 0 or len(endlines) == 0 or len(startlines) != len(endlines):
                return 'Invalid or not log file'

            srows = []
            for iter in range(len(startlines)):
                if iter == 0:
                    srows += [i for i in range(0,startlines[iter])]
                else:
                    srows += [i for i in range(endlines[iter-1],startlines[iter])]
                if iter == len(startlines)-1:    
                    srows += [i for i in range(endlines[iter],filelength)]
                else:
                    srows += [i for i in range(endlines[iter],startlines[iter+1])]

            if transientNum > 0:
                srows.remove(tmpS[0])

            df = pd.read_table(f, delim_whitespace=True, skiprows=srows)

            tmp_avgs = []
            for l in labels:
                tmp_avgs.append(df.loc[:,l].mean())
            avgs.append(tmp_avgs)
    
    else:
        filelength = fileLen(filename=data_file)
        tmpS = getStart(s='Per MPI rank memory allocation', filename=data_file)
        startline = (tmpS + transientNum)
        endline = getEnd(s='Loop time of', filename=data_file)

        if startline == -1:
            return 'Not Log File'

        if endline == -1:
            return 'Not Log File'

        srows = [i for i in range(0,startline)]+[i for i in range(endline,filelength)]
        if transientNum > 0:
            srows.remove(tmpS)

        df = pd.read_table(data_file, delim_whitespace=True, skiprows=srows)

        tmp_avgs = []
        for l in labels:
            tmp_avgs.append(df.loc[:,l].mean())
        avgs.append(tmp_avgs)

    avgs.sort(key=lambda x: x[0])

    s = ""
    for l in labels:
        s += str(l)+"\t"
    print(s)
    
    for x in avgs:
        s = ""
        for y in x:
            s += str(y)+"\t"
        print(s)

if __name__ == "__main__":
    main()