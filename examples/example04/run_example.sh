#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use pw.x to simulate a charged two-dimensional "
$ECHO "system immersed in continuum solvent described by the SCCS model "
$ECHO
$ECHO "   O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, 064102 (2012) "
$ECHO
$ECHO "and with open boundary conditions along the axis perpendicular "
$ECHO "to the slab plane, as described in "
$ECHO
$ECHO "   O. Andreussi and N. Marzari, Phys. Rev. B 90, 245101 (2014) "

# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x"
PSEUDO_LIST="O.pbe-rrkjus.UPF C.pbe-rrkjus.UPF Pt.pbe-nd-rrkjus.UPF"

$ECHO
$ECHO "  executables directory: $BIN_DIR"
$ECHO "  pseudo directory:      $PSEUDO_DIR"
$ECHO "  temporary directory:   $TMP_DIR"
$ECHO "  checking that needed directories and files exist...\c"

# check for directories
for DIR in "$BIN_DIR" "$PSEUDO_DIR" ; do
    if test ! -d $DIR ; then
        $ECHO
        $ECHO "ERROR: $DIR not existent or not a directory"
        $ECHO "Aborting"
        exit 1
    fi
done
for DIR in "$TMP_DIR" "$EXAMPLE_DIR/results" ; do
    if test ! -d $DIR ; then
        mkdir $DIR
    fi
done
cd $EXAMPLE_DIR/results

# check for executables
for FILE in $BIN_LIST ; do
    if test ! -x $BIN_DIR/$FILE ; then
        $ECHO
        $ECHO "ERROR: $BIN_DIR/$FILE not existent or not executable"
        $ECHO "Aborting"
        exit 1
    fi
done

# check for pseudopotentials
for FILE in $PSEUDO_LIST ; do
    if test ! -r $PSEUDO_DIR/$FILE ; then
       $ECHO
       $ECHO "Downloading $FILE to $PSEUDO_DIR...\c"
            $WGET $PSEUDO_DIR/$FILE \
                http://www.quantum-espresso.org/upf_files/UPF/$FILE 2> /dev/null 
    fi
    if test $? != 0; then
        $ECHO
        $ECHO "ERROR: $PSEUDO_DIR/$FILE not existent or not readable"
        $ECHO "Aborting"
        exit 1
    fi
done
$ECHO " done"

# how to run executables
PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX --environ"
$ECHO
$ECHO "  running pw.x as: $PW_COMMAND"
$ECHO

### ELECTROSTATIC EMBEDDING PARAMETERS ##############################
verbose=0             # if GE 1 prints debug informations
                      # if GE 2 prints out gaussian cube files with 
                      # dielectric function, polarization charges, etc
                      # WARNING: if GE 2 lot of I/O, much slower
environ_thr='1.d0'    # electronic convergence threshold for the onset  
                      # of solvation correction
environ_type='input'  # type of environment
                      # input: read parameters from input
                      # vacuum: all flags off, no environ 
                      # water: parameters from experimental values 
                      #        and specifically tuned
solvent_mode='full'   # specify the charge density that is used to 
                      # build the dielectric cavity:
                      # electronic: use the electronic density (default)
                      # ionic: use a fictitious charge density calculated
                      #        for atomic-centered interlocking spheres, 
                      #        whose analytical expression is based on the
                      #        error function (see example02)
                      # full: same as electronic, but a small extra 
                      #       density is added on the nuclei to avoid
                      #       spurious holes in the density due to 
                      #       pseudopotentials. Only needed with some
                      #       atomic types and pp, in particular 
                      #       alogens and transition metals
### PERIODIC BOUNDARY CONDITIONS ####################################
env_electrostatic='.true.' # modify electrostatic embedding (required to
                           #   switch on PBC corrections in vacuum)
pbc_correction='parabolic' # correction scheme to remove PBC 
                           # none: periodic calculation, no correction 
                           # parabolic: quadratic real-space correction
pbc_dim=2                  # select the desired system dimensionality
                           #   0, 1 or 2: isolated, 1D or 2D system
pbc_axis=3                 # set the axis along the 1D direction or 
                           #   normal to the 2D plane (pbc_axis = 1, 2 
                           #   or 3 for x, y or z axis)
#####################################################################

for epsilon in 1 80 ; do 

    # clean TMP_DIR
    $ECHO "  cleaning $TMP_DIR...\c"
    rm -rf $TMP_DIR/*
    $ECHO " done"
    
    if [ $epsilon = "1" ]; then
      label="vacuum"
    else
      label="water"
    fi

    $ECHO "  running the scf calculation in $label"
    $ECHO "  with $pbc_correction periodic boundary correction"

  prefix=PtCO_${label}
  input=${prefix}'.in'
  output=${prefix}'.out'
  cat > $input << EOF 
 &CONTROL
   !
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR/'
   outdir = '$TMP_DIR/'
   prefix = '$prefix'
   tprnfor = .TRUE.
   verbosity = 'high'
   !
 /
 &SYSTEM
   !
   ibrav = 8
   celldm(1) = 10.6881
   celldm(2) = 0.866025
   celldm(3) = 3.95422
   nat = 10
   ntyp = 3
   ecutwfc = 35
   ecutrho = 280
   occupations = 'smearing'
   degauss = 0.03
   smearing = 'mv'
   nbnd = 80
   tot_charge = 1
   !
/
 &ELECTRONS
   !
   conv_thr = 5.D-9
   diagonalization = 'davidson'
   mixing_beta = 0.2
   electron_maxstep = 200
   !
 /
K_POINTS (automatic)
 1 1 1 0 0 0
ATOMIC_SPECIES
C 1 C.pbe-rrkjus.UPF
O 1 O.pbe-rrkjus.UPF
Pt 1 Pt.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS (bohr)
C        5.335084148   4.646723426  12.901029877
O        5.335009643   4.619623254  15.079854269
Pt       8.061327071   0.098057998   8.992142901
Pt       2.608989366   0.098058283   8.992140585
Pt       0.000036609   4.720846294   8.968756935
Pt       5.335159557   4.721612729   9.380196435
Pt       0.000041121   7.802951963   4.604626508
Pt       5.335161233   7.697749113   4.753489408
Pt       2.697860636   3.152173889   4.688412329
Pt       7.972463687   3.152174491   4.688415209
EOF
  cat > environ_${label}.in << EOF
 &ENVIRON
   !
   verbose = $verbose
   environ_thr = $environ_thr
   environ_type = '$environ_type'
   env_static_permittivity = $epsilon
   env_surface_tension = 0.D0
   env_pressure = 0.D0
   env_electrostatic = $env_electrostatic
   !
 /
 &BOUNDARY
   !
   solvent_mode = '$solvent_mode'
   !
 /
 &ELECTROSTATIC
   !
   pbc_correction = '$pbc_correction'
   pbc_dim = $pbc_dim                 
   pbc_axis = $pbc_axis   
   tol = 5.D-13
   !
 /
EOF
 
  cp environ_${label}.in environ.in   
  $PW_COMMAND < $input > $output 
  check_failure $?
  $ECHO " done"

done

evac=$(awk '/^!/ {en=$5}; END {print en}' PtCO_vacuum.out)
esol=$(awk '/^!/ {en=$5}; END {print en}' PtCO_water.out)
dgsol=$($ECHO "($esol+(-1)*$evac)*313.68" | bc -l) 

fermi_vac=$(awk '/the Fermi energy is/ {en=$5}; END {print en}'     PtCO_vacuum.out)
delta_vac=$(awk '/the Fermi energy shift/ {en=$11}; END {print en}' PtCO_vacuum.out)
fermi_vac_new=$($ECHO "$fermi_vac+$delta_vac" | bc -l)

fermi_sol=$(awk '/the Fermi energy is/ {en=$5}; END {print en}'     PtCO_water.out)
delta_sol=$(awk '/the Fermi energy shift/ {en=$11}; END {print en}' PtCO_water.out)
fermi_sol_new=$($ECHO "$fermi_sol+$delta_sol" | bc -l)

$ECHO "  Energy in vacuum                    = $evac  Ry        " >> results.txt
$ECHO "  Energy in solution                  = $esol  Ry        " >> results.txt
$ECHO "  Electrostatic Solvation Energy      = $dgsol Kcal/mol  " >> results.txt

$ECHO >> results.txt
$ECHO "  Fermi energy in vacuum   = $fermi_vac eV + $delta_vac eV = $fermi_vac_new eV " >> results.txt
$ECHO "  Fermi energy in solution = $fermi_sol eV + $delta_sol eV = $fermi_sol_new eV " >> results.txt

$ECHO
$ECHO "$EXAMPLE_DIR : done"
