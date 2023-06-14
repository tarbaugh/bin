import argparse

parser = argparse.ArgumentParser(description='LAMMPS basic data tool')
parser.add_argument('files')
parser.add_argument('cores')
args = parser.parse_args()
num_files = args.files
num_cores = args.cores

import os
from ase.io import read
from ase.calculators.cp2k import CP2K

CP2K.command="mpirun -np {} /zfshomes/tarbaugh/downloads/new_cp2k/cp2k-9.1/exe/local/cp2k_shell.psmp".format(num_cores)

inp_str = """
&FORCE_EVAL
   &DFT
      BASIS_SET_FILE_NAME BASIS_MOLOPT
      POTENTIAL_FILE_NAME GTH_POTENTIALS
      &SCF
         MAX_SCF    500
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
    &END SUBSYS
    &PRINT
      &STRESS_TENSOR ON
      &END STRESS_TENSOR
   &END PRINT
&END FORCE_EVAL"""

calc = CP2K(basis_set=None,basis_set_file=None,potential_file=None,
		max_scf=None,cutoff=None,pseudo_potential=None,xc=None,inp=inp_str)

if num_files > 1:
    for n in range(num_files):
        atoms = read("lammps_input{}".format(n), format="lammps-data", style="atomic")
        atoms.set_calculator(calc)
        atoms.set_pbc((1, 1, 1))
        atoms.wrap()
        atoms.get_forces()

        lattice = atoms.get_lattice()

        mlip_write = open("output_ef.cfg", "w")
        mlip_write.write("BEGIN_CFG\n")
        mlip_write.write(" Size\n")
        mlip_write.write("    "+str(len(atoms.get_chemical_symbols()))+"\n")
        mlip_write.write(" Supercell\n")
        mlip_write.write("    "+str(lattice[0])+"    0    0\n")
        mlip_write.write("    0    "+str(lattice[1])+"    0\n")
        mlip_write.write("    0    0    "+str(lattice[2])+"\n")
        mlip_write.write(" AtomData:  id type    cartes_x    cartes_y    cartes_z    fx    fy    fz\n")
        for i in range (0,cfg_size):
            mlip_write.write("    "+str(i+1)+"    "+str(t[i]-1)+"    "+str(pos_x[i])+"    "+str(pos_y[i])+"    "+str(pos_z[i])+"    "+str(f_x[i])+"    "+str(f_y[i])+"    "+str(f_z[i])+"\n")
            energy += epa[i]
        mlip_write.write(" Energy\n")
        mlip_write.write("    "+str(energy)+"\n")
        mlip_write.write("END_CFG\n")
        mlip_write.write("\n")

        mlip_write.close()  
        os.system("cat output_ef.dump >> ../train.cfg")

else:
    atoms = read("lammps_input")
    atoms.set_calculator(calc)
    atoms.get_forces()
    atoms.write("output_ef.dump", format="lammps-data")
    os.system("cat output_ef.dump >> ../train.cfg")