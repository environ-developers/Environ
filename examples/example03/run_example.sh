#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use pw.x to calculate the electrostatic "
$ECHO "contribution to the solvation energy of a charged isolated system "
$ECHO "(pyridine cation in water) using the SCCS model "
$ECHO
$ECHO "   O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, 064102 (2012) "
$ECHO
$ECHO "coupled with different periodic-boundary correction schemes "
$ECHO
$ECHO "   O. Andreussi and N. Marzari, Phys. Rev. B 90, 245101 (2014) " 

# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x"
PSEUDO_LIST="N.pbe-rrkjus.UPF C.pbe-rrkjus.UPF H.pbe-rrkjus.UPF"

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
                http://www.quantum-espresso.org/pseudo/1.3/UPF/$FILE 2> /dev/null 
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
environ_thr='1.d-1'   # electronic convergence threshold for the onset  
                      # of solvation correction
environ_type='input'  # type of environment
                      # input: read parameters from input
                      # vacuum: all flags off, no environ 
                      # water: parameters from experimental values 
                      #        and specifically tuned
### PERIODIC BOUNDARY CONDITIONS ####################################
assume_isolated='none' # correction scheme to remove pbc 
                       # none: periodic calculation, no correction
                       # martyna-tuckerman: on-the-fly correction of
                       #   the potential in G-space, fast and no error
                       #   but requires a cell at least 2xsize of sytem
                       # pcc: on-the-fly correction of the potential
                       #   in real space, error on energy decays as
                       #   (cell size)^-5
                       # NOTE: the pcc correction can be activated by
                       # setting pbc_correction = 'parabolic' and 
                       # pbc_dim = 0 in the Environ input file
                       # in the ELECTROSTATIC namelist (recommended). 
                       # For vacuum calculations, one needs to 
                       # additionally set env_electrostatic = .true. 
                       # in the ENVIRON namelist (see e.g. example01).
#####################################################################

for assume_isolated in none martyna-tuckerman pcc ; do 

for epsilon in 01 80 ; do 

    # clean TMP_DIR
    $ECHO "  cleaning $TMP_DIR...\c"
    rm -rf $TMP_DIR/*
    $ECHO " done"
    
    if [ $epsilon = "01" ]; then
      label="vacuum"
    else
      label="water"
    fi

    $ECHO "  running the scf calculation in $label"
    $ECHO "  with $assume_isolated periodic boundary correction"

  prefix=pyridine_${label}_${assume_isolated}
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
   ecutrho = 300
   ecutwfc = 30
   ibrav = 1
   celldm(1) = 30
   nat = 12
   ntyp = 3
   tot_charge = 1 
   assume_isolated = '$assume_isolated'
   !
/
 &ELECTRONS
   !
   conv_thr = 5.D-9
   diagonalization = 'davidson'
   mixing_beta = 0.4
   electron_maxstep = 200
   !
 /
 &IONS
   ion_dynamics    = 'bfgs'
 /
K_POINTS (gamma)
ATOMIC_SPECIES  
N   14      N.pbe-rrkjus.UPF
C   12      C.pbe-rrkjus.UPF
H    1      H.pbe-rrkjus.UPF
ATOMIC_POSITIONS (angstrom)
H       -2.418655660   0.000000226   0.001997570
N       -1.294053698   0.000000088   0.001601933
C       -0.652551721   1.189250860   0.001097976
C        0.729271251   1.209271597   0.000774029
C        1.425569237   0.000000099   0.000697925
C        0.729271382  -1.209271463   0.001403411
C       -0.652551525  -1.189250323   0.001876697
H       -1.276762352   2.083974595   0.000904571
H        1.248379957   2.168950315   0.000425216
H        2.518135313  -0.000000088   0.000010538
H        1.248380238  -2.168950127   0.001697259
H       -1.276761908  -2.083974094   0.002582665
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
   !
 /
 &BOUNDARY
 /
 &ELECTROSTATIC
   tol = 1.D-12
 /
EOF

  cp environ_${label}.in environ.in   
  $PW_COMMAND < $input > $output 
  check_failure $?
  $ECHO " done"

done

evac=$(awk '/^!/ {en=$5}; END {print en}' pyridine_vacuum_${assume_isolated}.out)
esol=$(awk '/^!/ {en=$5}; END {print en}' pyridine_water_${assume_isolated}.out)
dgsol=$($ECHO "($esol+(-1)*$evac)*313.68" | bc -l) 

$ECHO "  Periodic boundary correction scheme = $assume_isolated " >> results.txt
$ECHO "  Electrostatic Energy in vacuum      = $evac  Ry        " >> results.txt
$ECHO "  Electrostatic Energy in solution    = $esol  Ry        " >> results.txt
$ECHO "  Electrostatic Solvation Energy      = $dgsol Kcal/mol  " >> results.txt
$ECHO "                                                         " >> results.txt

done

$ECHO
$ECHO "$EXAMPLE_DIR : done"
