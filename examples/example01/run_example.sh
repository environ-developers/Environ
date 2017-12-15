#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use pw.x to calculate the solvation energy "
$ECHO "and other solvent related quantites for a water molecule in water"
$ECHO "using a fully self consistent dielectric defined on the electronic"
$ECHO "density according to "
$ECHO "   O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, 064102 (2012) "

# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x"
PSEUDO_LIST="O.pbe-rrkjus.UPF H.pbe-rrkjus.UPF"

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
environ_type='vacuum' # type of environment
                      # input: read parameters from input
                      # vacuum: all flags off, no environ 
                      # water: parameters from experimental values 
                      #        and specifically tuned
### ITERATIVE SOLVER PARAMETERS #####################################
solver='iterative' # iterative solver
auxiliary='full'   # solve for the polarization charge density
tol='1.d-11'       # tolerance of the iterative solver
mix='0.6'          # polarization mixing for the iterative solver
#####################################################################

for environ_type in vacuum water ; do 

    # clean TMP_DIR
    $ECHO "  cleaning $TMP_DIR...\c"
    rm -rf $TMP_DIR/*
    $ECHO " done"
    
    $ECHO "  running the relax calculation in $environ_type"

  prefix=h2o_$environ_type
  input=${prefix}'.in'
  output=${prefix}'.out'
  cat > $input << EOF 
 &CONTROL
   !
   calculation = 'relax'
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
   celldm(1) = 20
   nat = 3
   ntyp = 2
   assume_isolated = 'pcc'
   !
/
 &ELECTRONS
   !
   conv_thr = 5.D-9
   diagonalization = 'cg'
   mixing_beta = 0.4
   electron_maxstep = 200
   !
 /
 &IONS
   ion_dynamics    = 'bfgs'
 /
K_POINTS (automatic)
 1 1 1 0 0 0
ATOMIC_SPECIES  
 H   1  H.pbe-rrkjus.UPF
 O  16  O.pbe-rrkjus.UPF
ATOMIC_POSITIONS (bohr)
O   11.79  12.05  11.50
H   13.45  11.22  11.50
H   10.56  10.66  11.50
EOF
  cat > environ_${environ_type}.in << EOF
 &ENVIRON
   !
   verbose = $verbose
   environ_thr = $environ_thr
   environ_type = '$environ_type'
   !
 /
 &BOUNDARY
 /
 &ELECTROSTATIC
   !
   tol = $tol
   mix = $mix
   solver = '$solver'
   auxiliary = '$auxiliary'
   !  
 /
EOF
   
  cp environ_${environ_type}.in environ.in
  $PW_COMMAND < $input > $output 
  check_failure $?
  $ECHO " done"

done

evac=$(awk '/^!/ {en=$5}; END {print en}' h2o_vacuum.out)
esol=$(awk '/^!/ {en=$5}; END {print en}' h2o_water.out)
dgsol=$($ECHO "($esol+(-1)*$evac)*313.68" | bc -l) 
ecav=$(awk 'BEGIN {en=0}; /cavitation energy/ {en=$4}; END {print en}' h2o_water.out) 
epres=$(awk 'BEGIN {en=0}; /PV energy/ {en=$4}; END {print en}' h2o_water.out)

$ECHO "  Solvation Energy     = $dgsol Kcal/mol" > results.txt
iprint=0
dgelec=$dgsol
if [ $ecav != 0 ]; then 
  iprint=1
  dgcav=$($ECHO "$ecav*313.68" | bc -l)
  $ECHO "  Cavitation Energy    =  $dgcav Kcal/mol" >> results.txt
  dgelec=$($ECHO "$dgelec+(-1)*$dgcav" | bc -l)
fi
if [ $epres != 0 ]; then
  iprint=1
  dgpres=$($ECHO "$epres*313.68" | bc -l)
  $ECHO "  PV Energy            = $dgpres Kcal/mol" >> results.txt
  dgelec=$($ECHO "$dgelec+(-1)*$dgpres" | bc -l)
fi
if [ $iprint != 0 ]; then
  $ECHO "  Electrostatic Energy = $dgelec Kcal/mol" >> results.txt
fi

$ECHO
$ECHO "$EXAMPLE_DIR : done"
