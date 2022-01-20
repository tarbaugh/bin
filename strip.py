import numpy as np
import argparse

parser = argparse.ArgumentParser(description='LAMMPS data extraction tool')
parser.add_argument('log', nargs='+')

def fileLen(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return (i+1)

def getStarts(s='Step', filename=''):
    starts = []
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if s in l:
                starts.append((i+1, getEnd(filename=filename, startAt=i), l))
    return starts

def getEnd(s='Loop time of', filename='', startAt=0):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if i >= startAt:
                if s in l:
                    return i
    return -1

def main():
    args = parser.parse_args()
    logfiles = args.log

    logfile = logfiles[0]
    startEnds = getStarts(filename=logfile) # + transientNum

    if len(startEnds) == 0:
        return 'Not Log File'

    for _, i_check, _ in startEnds:
        if i_check == -1:
            return 'Log File Error'

    # srows = [i for i in range(0,startline)]+[i for i in range(endline,filelength)]
    # df = pd.read_table(logfile, delim_whitespace=True, skiprows=srows)

    s, f, first_h = startEnds[0]
    print("Collected: ", first_h[:-2], "...")
    df = np.genfromtxt(logfile, skip_header=s, max_rows=(f-s))

    i = 1
    while i < len(startEnds):
        s, f, h = startEnds[i]
        if h != first_h:
            print("Appended data has different labels:\n", h)
            break
        s += 1
        df = np.append(df, np.genfromtxt(logfile, skip_header=s, max_rows=(f-s)), axis=0)
        i += 1

    file_iter = 1
    while file_iter < len(logfiles):

        logfile = logfiles[file_iter]
        startEnds = getStarts(filename=logfile) # + transientNum

        if len(startEnds) == 0:
            return 'Not Log File'

        for _, i_check, _ in startEnds:
            if i_check == -1:
                return 'Log File Error'

        # srows = [i for i in range(0,startline)]+[i for i in range(endline,filelength)]
        # df = pd.read_table(logfile, delim_whitespace=True, skiprows=srows)

        i = 0
        while i < len(startEnds):
            s, f, h = startEnds[i]
            if h != first_h:
                print("Appended data has different labels:\n", h)
                break
            s += 1
            df = np.append(df, np.genfromtxt(logfile, skip_header=s, max_rows=(f-s)), axis=0)
            i += 1

        file_iter += 1
        
    if len(logfiles) > 1:
        np.savetxt('aggregate_logs.dat', df, fmt='%lf', header=startEnds[0][2][:-2])
        print("File saved as ", 'aggregate_logs.dat')
    else:
        np.savetxt(logfile+'.dat', df, fmt='%lf', header=startEnds[0][2][:-2])
        print("File saved as ", logfile+'.dat')

if __name__ == "__main__":
    main()
