import random
import numpy as np
from ase.lattice.cubic import SimpleCubic
from ase.io import read


ase_atoms = read('vdw_796.xyz', index=-1)
# ase_atoms = SimpleCubic(
      #   symbol='Te',
      #   size=(6, 6, 6),
      #   pbc=True,
      #   latticeconstant=3.28537)

# old_p = ase_atoms.get_positions()
# old_p += 0.25*np.random.random_sample(old_p.shape)-0.125
# ase_atoms.set_positions(old_p)
ase_atoms.set_pbc((True,True,True))
ase_atoms.wrap()

pos = ase_atoms.get_positions()
cell = ase_atoms.get_cell()

indexs = [x for x in range(len(pos))]
random.shuffle(indexs)

atoms = []
for c in range(31):
    atoms.append(("Ge", pos[indexs[c]][0], pos[indexs[c]][1], pos[indexs[c]][2]))

for c in range(31, 216):
    atoms.append(("Te", pos[indexs[c]][0], pos[indexs[c]][1], pos[indexs[c]][2]))

with open('cp2k.inp', 'w') as f:
    f.write("""
&GLOBAL
   PROJECT GeTeMD
   RUN_TYPE MD
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
      &QS
         EPS_DEFAULT 1.0E-12
      &END QS
      &SCF
         MAX_SCF    500
         EPS_SCF    1.0E-06
         SCF_GUESS  ATOMIC
         &MIXING
            METHOD BROYDEN_MIXING
         &END MIXING
         &SMEAR ON
            METHOD FERMI_DIRAC
            ELECTRONIC_TEMPERATURE 300
         &END SMEAR
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
    for a in atoms:
        f.write("         {} {} {} {}\n".format(a[0], a[1], a[2], a[3]))
    f.write("      &END COORD\n")
    f.write("      &CELL\n")
    f.write("         PERIODIC XYZ\n")
    f.write("         A {} 0.000000000000000000e+00 0.000000000000000000e+00\n".format(cell[0][0]))
    f.write("         B 0.000000000000000000e+00 {} 0.000000000000000000e+00\n".format(cell[1][1]))
    f.write("         C 0.000000000000000000e+00 0.000000000000000000e+00 {}\n".format(cell[2][2]))
    f.write("      &END CELL\n")
    f.write(
"""     &END SUBSYS
   &PRINT
      &STRESS_TENSOR ON
      &END STRESS_TENSOR
   &END PRINT
&END FORCE_EVAL

&MOTION
   &MD
      ENSEMBLE NVT
      STEPS 7500
      TIMESTEP 2.0
      TEMPERATURE 1150.0
      &THERMOSTAT
         &NOSE                    #Uses the Nose-Hoover thermostat
            TIMECON 200           #timeconstant of the thermostat chain, how often does thermostat adjust your system 
         &END NOSE
      &END
   &END MD
   &PRINT
      &TRAJECTORY
         &EACH
            MD 1
         &END EACH
      &END TRAJECTORY
      &FORCES ON
      &END
      &RESTART
         BACKUP_COPIES 1
      &END RESTART
   &END PRINT
&END MOTION
""")