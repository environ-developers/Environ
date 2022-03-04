#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use neb.x to calculate the"
$ECHO "minimum energy path (MEP) of the collinear proton transfer reaction:"
$ECHO "  H2+H <==> H+H2, within the Born-Oppenheimer approximation."
$ECHO
$ECHO "!!! Beware: neb.x DOES NOT READ FROM STANDARD INPUT"
$ECHO "!!! run as 'neb.x -inp input_file_name > output_file_name'"
$ECHO
# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="neb.x"
PSEUDO_LIST="HUSPBE.RRKJ3"

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
            $WGET $PSEUDO_DIR/$FILE $NETWORK_PSEUDO/$FILE 2> /dev/null
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
NEB_COMMAND="$PARA_PREFIX $BIN_DIR/neb.x $PARA_POSTFIX"
$ECHO
$ECHO "  running Born-Oppenheimer NEB as: $NEB_COMMAND"
$ECHO

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/pwscf*
$ECHO " done"

# NEB calculation. Automatic choice of the climbing image.
cat > H2+H.in << EOF
BEGIN
BEGIN_PATH_INPUT
&PATH
  restart_mode      = 'from_scratch'
  string_method     = 'neb',
  nstep_path        = 20,
  ds                = 2.D0,
  opt_scheme        = "broyden",
  num_of_images     = 7,
  k_max             = 0.3D0,
  k_min             = 0.2D0,
  CI_scheme         = "auto",
  path_thr          = 0.1D0,
/
END_PATH_INPUT
BEGIN_ENGINE_INPUT
&CONTROL
  prefix         = "H2+H"
  outdir         = "$TMP_DIR",
  pseudo_dir     = "$PSEUDO_DIR",
/
&SYSTEM
  ibrav                  = 0,
  nat                    = 3,
  ntyp                   = 1,
  ecutwfc                = 20.0D0,
  ecutrho                = 100.0D0,
  nspin                  = 2,
  starting_magnetization = 0.5D0,
  occupations            = "smearing",
  degauss                = 0.03D0,
/
&ELECTRONS
  conv_thr    = 1.D-8,
  mixing_beta = 0.3D0,
/
&IONS
/
ATOMIC_SPECIES
H  1.00794  HUSPBE.RRKJ3
BEGIN_POSITIONS
FIRST_IMAGE
ATOMIC_POSITIONS { bohr }
H   -4.56670009      0.00000000      0.00000000  1  0  0
H    0.00000000      0.00000000      0.00000000  0  0  0
H    1.55776676      0.00000000      0.00000000  1  0  0
LAST_IMAGE
ATOMIC_POSITIONS { bohr }
H   -1.55776676      0.00000000      0.00000000
H    0.00000000      0.00000000      0.00000000
H    4.56670009      0.00000000      0.00000000
END_POSITIONS
K_POINTS { gamma }
CELL_PARAMETERS { bohr }
  12.00000  0.00000  0.00000
   0.00000  5.00000  0.00000
   0.00000  0.00000  5.00000
END_ENGINE_INPUT
END
EOF

$ECHO

# NEB
$ECHO "  running Born-Oppenheimer NEB calculation for H2+H => H+H2... \c"
$NEB_COMMAND -inp H2+H.in  > H2+H.out
check_failure $?
$ECHO "done"

$ECHO

# Environ input
cat > environ.in << EOF
 &ENVIRON
   environ_thr = 1.d2
   env_pressure = 10. 
 /
 &BOUNDARY
   stype = 1
 /
EOF

# NEB with Environ
$ECHO "  running Born-Oppenheimer NEB calculation for H2+H => H+H2 under pressure... \c"
$NEB_COMMAND -inp H2+H.in --environ > H2+H_pressure.out
check_failure $?
$ECHO "done"

$ECHO

# Environ input
cat > environ.in << EOF
 &ENVIRON
   environ_thr = 1.d2
   env_static_permittivity = 1.d2
 /
 &BOUNDARY
 /
 &ELECTROSTATIC
 /
EOF

# NEB with Environ
$ECHO "  running Born-Oppenheimer NEB calculation for H2+H => H+H2 in dielectric... \c"
$NEB_COMMAND -inp H2+H.in --environ > H2+H_dielectric.out
check_failure $?
$ECHO "done"

$ECHO

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR... \c"
rm -rf $TMP_DIR/pwscf*
$ECHO "done"
$ECHO
$ECHO "$EXAMPLE_DIR: done"
