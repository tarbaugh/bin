def readCFGForces(file, mtp=False):
    with open(file) as f:
        mark = 99999999999
        tot_forces = []
        forces = []
        count = 0

        for i, l in enumerate(f):
            if l.startswith(" AtomData:"):
                mark = int(i)
                count += 1
            if i == mark+201:
                mark = 99999999999
                tot_forces.append(forces)
                forces = []
            if i >= mark+1:
                tmp_forces = l.split()[-3:]
                forces.append([float(j) for j in tmp_forces])

    return tot_forces

def readCFGCells(file, mtp=False):
    with open(file) as f:
        mark = 99999999999
        tot_cells = []
        cells = []
        count = 0

        for i, l in enumerate(f):
            if l.startswith(" Supercell"):
                mark = int(i)
                count += 1
            if i == mark+4:
                mark = 99999999999
                tot_cells.append(cells)
                cells = []
            if i >= mark+1:
                tmp_cells = l.split()
                cells.append([float(j) for j in tmp_cells])

    return tot_cells

def readCFGPosTypes(file, mtp=False):
    with open(file) as f:
        mark = 99999999999
        tot_pos = []
        pos = []
        types = []
        tot_types = []
        count = 0

        for i, l in enumerate(f):
            if l.startswith(" AtomData:"):
                mark = int(i)
                count += 1
            if i == mark+201:
                mark = 99999999999
                tot_pos.append(pos)
                pos = []
                tot_types.append(types)
                types = []
            if i >= mark+1:
                tmp_pos = l.split()[2:5]
                tmp_type = l.split()[1]
                types.append(int(tmp_type)+1)
                pos.append([float(j) for j in tmp_pos])

    return tot_pos, tot_types
def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f')
    args = parser.parse_args()

    pos, types = readCFGPosTypes(args.f)
    # forces = readCFGForces(args.f)
    cells = readCFGCells(args.f)

    from ase import Atoms
    from ase.io import write
    atoms = []
    for i in range(len(types)):
        a = Atoms(positions=pos[i],
            numbers=types[i],
            cell=cells[i],
            pbc=[1,1,1])
        a.wrap()
        atoms.append(a)

    write(filename="cfgs.xyz", images=atoms, format="extxyz")

if __name__ == "__main__":
    main()