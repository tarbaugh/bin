import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('-t', type=int, nargs='?', default=0)
parser.add_argument('log')

args = parser.parse_args()

logfile = args.log
transientNum = args.t


def main():

    with open(logfile) as f:
        for i, l in enumerate(f):
            if i > 0:
                break
            s = str(l)
    s = s.split()        
    df = pd.read_table(logfile, delim_whitespace=True, names=s[1:], skiprows=1)

    print(df.describe())

if __name__ == "__main__":
    main()