from ase.calculators.cp2k import CP2K
from ase.build import molecule

CP2K.command="env OMP_NUM_THREADS=2 mpiexec -np 16 /home/tphyswes/cp2k-9.1/exe/local/cp2k_shell.psmp"

inp_str = '''
&FORCE_EVAL
METHOD Quickstep
    &DFT
        &QS
            EPS_DEFAULT 1.0E-10
        &END QS
        &KPOINTS
            SCHEME GAMMA 12 12 8
        &END KPOINTS
        &SCF
            ADDED_MOS 10
            &SMEAR
                METHOD FERMI_DIRAC
                ELECTRONIC_TEMPERATURE [K] 500.0
            &END SMEAR
        &END SCF
        &XC
            &XC_FUNCTIONAL
                &GGA_X_RPW86
                &END GGA_X_RPW86
                &PW92
                &END PW92
            &END XC_FUNCTIONAL
            &vdW_POTENTIAL
                DISPERSION_FUNCTIONAL NON_LOCAL
                &NON_LOCAL
                    TYPE LMKLL
                    VERBOSE_OUTPUT
                    KERNEL_FILE_NAME vdW_kernel_table.dat
                    CUTOFF 300
                &END NON_LOCAL
            &END vdW_POTENTIAL
        &END XC
    &END DFT
    &PRINT
        &FORCES ON
        &END FORCES
  &END PRINT
&END FORCE_EVAL
'''

calc = CP2K(inp=inp_str)
tot_trajs = []
atoms = molecule('H2O', calculator=calc)
atoms.center(vacuum=2.0)
print(atoms.get_potential_energy())