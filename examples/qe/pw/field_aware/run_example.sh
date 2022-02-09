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
$ECHO "solvation energy for a water molecule in water using a field-aware"
$ECHO " soft sphere interface defined on the electrostatic field at the "
$ECHO "interface normal and the atomic positions, according to "
$ECHO "M Truscott, O Andreussi, J. Phys. Chem. B 2019, 123, 16, 3513–3524"

# set the needed environment variables
x=$EXAMPLE_DIR
while test "$x" != "/"; do
    x=$(dirname "$x")
    if test -f "$x/environment"; then
      . "$x/environment"; break
    fi
done

# compatibility with QE for versions prior to 6.4
if [ -z $NETWORK_PSEUDO ]; then
    NETWORK_PSEUDO=http://www.quantum-espresso.org/wp-content/uploads/upf_files/
fi

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
                $NETWORK_PSEUDO/$FILE 2> /dev/null
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

### SYSTEM PARAMETERS ###############################################
ecutrho=300
ecutwfc=30
### ELECTROSTATIC EMBEDDING PARAMETERS ##############################
verbose=0             # if GE 1 prints debug informations
                      # if GE 2 prints out gaussian cube files with
                      # dielectric function, polarization charges, etc
                      # WARNING: if GE 2 lot of I/O, much slower
### ITERATIVE SOLVER PARAMETERS #####################################
tol='1.d-10'  # tolerance of the iterative solver
#####################################################################

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/*
$ECHO " done"

$ECHO "  running the relax calculations with field aware"

cat > h2o.in << EOF
&CONTROL
   calculation = 'relax'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR'
   outdir = '$TMP_DIR'
   prefix = 'h2o'
   tprnfor = .TRUE.
   verbosity = 'high'
   forc_conv_thr = 5.d-3
/
&SYSTEM
   ecutrho = $ecutrho
   ecutwfc = $ecutwfc
   ibrav = 1
   celldm(1) = 20
   nat = 3
   ntyp = 2
/
&ELECTRONS
   electron_maxstep = 200
   conv_thr = 5.d-9
/
&IONS
   ion_dynamics = 'bfgs'
/
&CELL
/
K_POINTS (automatic)
1 1 1 0 0 0
ATOMIC_SPECIES
H   1.008  H.pbe-rrkjus.UPF
O  15.999  O.pbe-rrkjus.UPF
ATOMIC_POSITIONS (bohr)
O   10.28  10.51  10.00
H   11.87   9.57  10.00
H    8.99   9.20  10.00
EOF

cat > h3o.in << EOF
&CONTROL
   calculation = 'relax'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR'
   outdir = '$TMP_DIR'
   prefix = 'h3o'
   tprnfor = .TRUE.
   verbosity = 'high'
   forc_conv_thr = 5.d-2
/
&SYSTEM
   ecutrho = $ecutrho
   ecutwfc = $ecutwfc
   ibrav = 1
   smearing = "mv"
   degauss = 0.01
   celldm(1) = 20
   nat = 4
   ntyp = 2
   tot_charge = 1
/
&ELECTRONS
   electron_maxstep = 200
   conv_thr = 5.d-9
/
&IONS
   ion_dynamics = 'bfgs'
/
&CELL
/
K_POINTS (automatic)
1 1 1 0 0 0
ATOMIC_SPECIES
H   1.008  H.pbe-rrkjus.UPF
O  15.999  O.pbe-rrkjus.UPF
ATOMIC_POSITIONS (bohr)
H  12.23 11.50 10.00
H  10.56  8.75 10.00
H   9.02 11.57 10.00
O  10.61 10.61 10.00
EOF

cat > ho.in << EOF
&CONTROL
   calculation = 'relax'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR'
   outdir = '$TMP_DIR'
   prefix = 'ho'
   tprnfor = .TRUE.
   verbosity = 'high'
   forc_conv_thr = 5.d-3
/
&SYSTEM
   ecutrho = $ecutrho
   ecutwfc = $ecutwfc
   ibrav = 1
   smearing = "mv"
   degauss = 0.01
   celldm(1) = 20
   nat = 2
   ntyp = 2
   tot_charge = -1
/
&ELECTRONS
   electron_maxstep = 200
   conv_thr = 5.d-9
/
&IONS
   ion_dynamics = 'bfgs'
/
&CELL
/
K_POINTS (automatic)
1 1 1 0 0 0
ATOMIC_SPECIES
H   1.008  H.pbe-rrkjus.UPF
O  15.999  O.pbe-rrkjus.UPF
ATOMIC_POSITIONS (bohr)
O   10.00   10.00   10.00
H   11.84   10.00   10.00
EOF

cat > vacuum.in << EOF
&ENVIRON
  env_electrostatic = .true.
  env_pressure = 0.d0
  env_static_permittivity = 1.d0
  env_surface_tension = 0.d0
  environ_restart = .false.
  environ_thr = 1.d-1
  environ_type = 'input'
  verbose = $verbose
/
&BOUNDARY
  alpha = 1.12d0
  solvent_mode = 'ionic'
  radius_mode = 'muff'
/
&ELECTROSTATIC
  auxiliary = 'none'
  pbc_correction = 'parabolic'
  pbc_dim = 0
  solver = 'direct'
  tol = $tol
/
EOF

cat > solution.in << EOF
&ENVIRON
  env_electrostatic = .true.
  env_pressure = -3.5d-1
  env_static_permittivity = 7.83d1
  env_surface_tension = 5.d1
  environ_restart = .false.
  environ_thr = 1.d-1
  environ_type = 'input'
  verbose = $verbose
/
&BOUNDARY
  alpha = 1.12d0
  radius_mode = 'muff'
  solvent_mode = 'ionic'
  field_aware = .TRUE.
/
&ELECTROSTATIC
  auxiliary = 'none'
  pbc_correction = 'parabolic'
  pbc_dim = 0
  solver = 'cg'
  tol = $tol
/
EOF

# VACUUM
cp vacuum.in environ.in
$ECHO "performing calculations in vacuum..."
for input in ho h2o h3o; do
    $ECHO "calculating energy for $input"
    $PW_COMMAND < $input.in > ${input}_vacuum.out
    check_failure $?
done
cp solution.in environ.in
$ECHO "performing calculations in solution..."
for input in ho h2o h3o; do
    $ECHO "calculating energy for $input"
    $PW_COMMAND < $input.in > ${input}_solution.out
    check_failure $?
done

for input in ho h2o h3o; do
    evac=$(awk '/^!/ {en=$5}; END {print en}' ${input}_vacuum.out)
    esol=$(awk '/^!/ {en=$5}; END {print en}' ${input}_solution.out)
    dgsol=$($ECHO "($esol+(-1)*$evac)*313.68" | bc -l)
    $ECHO "  Input                = $input "         >> results.txt
    $ECHO "  Solvation Energy     = $dgsol Kcal/mol" >> results.txt
    $ECHO
done

$ECHO
$ECHO "$EXAMPLE_DIR : done"
