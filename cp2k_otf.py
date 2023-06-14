import os
import argparse
import subprocess
import numpy as np
from ase.io import read

def run_cp2k(positions, types, cells, num_cores):
    with open('cp2k.inp', 'w') as f:
        f.write(
"""&GLOBAL
PROJECT GeTeMD
RUN_TYPE ENERGY_FORCE
PRINT_LEVEL LOW
&TIMINGS
    THRESHOLD 0.000001
&END TIMINGS
&END GLOBAL

&FORCE_EVAL
   METHOD Quickstep
   STRESS_TENSOR ANALYTICAL
   &DFT
      BASIS_SET_FILE_NAME BASIS_MOLOPT
      POTENTIAL_FILE_NAME GTH_POTENTIALS
      &SCF
         MAX_SCF    3
         SCF_GUESS  ATOMIC
         CHOLESKY   INVERSE
         ADDED_MOS  -1
         &SMEAR ON
            METHOD FERMI_DIRAC
            ELECTRONIC_TEMPERATURE 300
         &END SMEAR
         &MIXING
            METHOD BROYDEN_MIXING
         &END MIXING
         &PRINT
            &RESTART OFF
            &END RESTART
         &END PRINT
      &END SCF
      &MGRID
         NGRIDS        4
         CUTOFF        40
         REL_CUTOFF    60
      &END MGRID
      &XC
         &XC_FUNCTIONAL BLYP
         &END XC_FUNCTIONAL
      &END XC
   &END DFT
   &SUBSYS
      &KIND Te
         BASIS_SET DZVP-MOLOPT-SR-GTH-q6
         POTENTIAL GTH-BLYP-q6
      &END KIND
      &KIND Ge
         BASIS_SET DZVP-MOLOPT-SR-GTH-q4
         POTENTIAL GTH-BLYP-q4
      &END KIND
    """)  
        f.write("      &COORD\n")
        for i, a in enumerate(positions):
            f.write("         {} {} {} {}\n".format(types[i], a[0], a[1], a[2]))
        f.write("  &END COORD\n")
        f.write("      &CELL\n")
        f.write("         PERIODIC XYZ\n")
        f.write("         A {} {} {}\n".format(cells[0][0], cells[0][1], cells[0][2]))
        f.write("         B {} {} {}\n".format(cells[1][0], cells[1][1], cells[1][2]))
        f.write("         C {} {} {}\n".format(cells[2][0], cells[2][1], cells[2][2]))
        f.write("      &END CELL\n")
        f.write(
"""     &END SUBSYS
    &PRINT
        &STRESS_TENSOR ON
        &END STRESS_TENSOR
        &FORCES ON
        &END FORCES
    &END PRINT
&END FORCE_EVAL
    """)

    cmd = 'mpirun -np {} /zfshomes/tarbaugh/downloads/new_cp2k/cp2k-9.1/exe/local/cp2k.popt -i cp2k.inp -o cp2k.out'.format(num_cores)
    process = subprocess.Popen(
        cmd, 
        shell=True, 
        stdout=subprocess.PIPE)
    process.wait()

def gen_data_from_lammps(filename):
    atoms = read(filename, format="lammps-data", style="atomic")
    atoms.set_pbc((1, 1, 1))
    atoms.wrap()
    cells = atoms.get_cell()
    positions = atoms.get_positions()
    types = []
    for c in atoms.get_chemical_symbols():
        if c == 'H':
            types.append('Ge')
        else:
            types.append('Te')
    return positions, types, cells

def gather_strss_frcs(length):
    frcs = np.zeros((length, 3))
    with open('cp2k.out', 'r') as f:
        pre_frcs = 9999999
        forces_done = False
        j = 0
        for i, l in enumerate(f):
            if l.startswith("  Total energy:"):
                tot_eng = float(l.split()[-1])
            
            if l.startswith(" # Atom   Kind   Element"):
                pre_frcs = i
            
            if i > pre_frcs and not forces_done:
                if l.startswith(" SUM OF ATOMIC FORCES"):
                    forces_done = True
                else:
                    spl = l.split()
                    frcs[j][0] = float(spl[3])
                    frcs[j][1] = float(spl[4])
                    frcs[j][2] = float(spl[5])
                    j += 1
            
            if l.startswith(" STRESS|      x"):
                spl = l.split()
                stress_xx = float(spl[2])
                stress_xy = float(spl[3])
                stress_xz = float(spl[4])
            
            if l.startswith(" STRESS|      y"):
                spl = l.split()
                stress_yy = float(spl[3]) 
                stress_yz = float(spl[4])
            
            if l.startswith(" STRESS|      z"):
                spl = l.split()
                stress_zz = float(spl[4])
                break

        stresses = [stress_xx, stress_yy, stress_zz, stress_yz, stress_xz, stress_xy]
        return tot_eng, frcs, stresses

def write_mlip(positions, types, cells, forces, energy, stresses):
    with open('new_cfgs.cfg', 'a') as file:
        file.write("BEGIN_CFG\n")
        file.write(" Size\n")
        file.write("    {}\n".format(len(types)))
        file.write(" Supercell\n")
        file.write("    {}    {}    {}\n".format(cells[0][0], cells[0][1], cells[0][2]))
        file.write("    {}    {}    {}\n".format(cells[1][0], cells[1][1], cells[1][2]))
        file.write("    {}    {}    {}\n".format(cells[2][0], cells[2][1], cells[2][2]))
        file.write(" AtomData:  id type    cartes_x    cartes_y    cartes_z    fx    fy    fz\n")
        for n, curr_type in enumerate(types):
            if curr_type == 'Te':
                t = 1
            else:
                t = 0
            file.write("    {}    {}    {}    {}    {}    {}    {}    {}\n".format(n+1, t, positions[n][0], positions[n][1], positions[n][2], 
                forces[n][0]*51.42208619083232, forces[n][1]*51.42208619083232, forces[n][2]*51.42208619083232))
        file.write(" Energy\n")
        file.write("    {}\n".format(energy))
        file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        file.write("     {} {} {} {} {} {}\n".format(stresses[0], stresses[1], stresses[2], stresses[3], stresses[4], stresses[5]))
        file.write("END_CFG\n")
        file.write("\n")

def clean_up():
    cmd = 'cat cp2k.out >> tot.out && rm cp2k.out'
    process = subprocess.Popen(
        cmd, 
        shell=True, 
        stdout=subprocess.PIPE)
    process.wait()

def main():
    parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
    parser.add_argument('files')
    parser.add_argument('cores')
    args = parser.parse_args()
    num_files = int(args.files)
    num_cores = int(args.cores)

    if num_files > 1:
        for n in range(num_files):
            pos, t, cell = gen_data_from_lammps("lammps_input{}".format(n))
            run_cp2k(pos, t, cell, num_cores)
            enrg, frcs, strss = gather_strss_frcs(len(t))
            write_mlip(pos, t, cell, frcs, enrg, strss)
            clean_up()
    else:
        pos, t, cell = gen_data_from_lammps("lammps_input")
        run_cp2k(pos, t, cell, num_cores)
        enrg, frcs, strss = gather_strss_frcs(len(t))
        write_mlip(pos, t, cell, frcs, enrg, strss)
        clean_up()

if __name__ == "__main__":
    main()