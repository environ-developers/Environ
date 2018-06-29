#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use pw.x, turbo_lanczos.x, and turbo_spectrum.x"
$ECHO "to calculate the absorption spectrum of the CH4 molecule in water"
$ECHO "using a fully self-consistent dielectric defined on the electronic density"
$ECHO "and using the Liouville-Lanczos approach within the time-dependent density" 
$ECHO "functional perturbation theory according to I. Timrov, O. Andreussi, A. Biancardi,"
$ECHO "N. Marzari, and S. Baroni, J. Chem. Phys. 142, 034111 (2015)."

# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x turbo_lanczos.x turbo_spectrum.x"
PSEUDO_LIST="H.pz-vbc.UPF C.pz-vbc.UPF"

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
                http://www.quantum-espresso.org/upf_files/$FILE 2> /dev/null 
    fi
    if test $? != 0; then
        $ECHO
        $ECHO "ERROR: $PSEUDO_DIR/$FILE not existent or not readable"
        $ECHO "Aborting"
        exit 1
    fi
done
$ECHO " done"


### ELECTROSTATIC EMBEDDING PARAMETERS ##############################
verbose=0                      # if GE 1 prints debug informations
                               # if GE 2 prints out gaussian cube files with 
                               # dielectric function, polarization charges, etc
                               # WARNING: if GE 2 lot of I/O, much slower
environ_thr='1.d-1'            # electronic convergence threshold for the onset  
                               # of solvation correction
environ_type="vacuum"          # type of environment
                               # input: read parameters from input
                               # vacuum: all flags off, no environ 
env_static_permittivity=78.5   # static permittivity of water
env_optical_permittivity=1.776 # optical permittivity of water
env_surface_tension=0.0        # surface tension (not supported by TDDFPT)
env_pressure=0.0               # pressure (not supported by TDDFPT)
### ITERATIVE SOLVER PARAMETERS #####################################
tol='1.d-12' # tolerance of the iterative solver
#####################################################################

for environ_type in vacuum input ; do 

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/*
$ECHO " done"

if [ "${environ_type}" = "input" ] ; then
   prefix_environ="solvent"
   flag="-environ"
else
   prefix_environ="vacuum"
   flag=" "
fi

# how to run executables
PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x $PARA_POSTFIX $flag"
TURBO_LANCZOS_COMMAND="$PARA_PREFIX $BIN_DIR/turbo_lanczos.x $PARA_POSTFIX $flag"
TURBO_SPECTRUM_COMMAND="$PARA_PREFIX $BIN_DIR/turbo_spectrum.x $PARA_POSTFIX $flag"
$ECHO
$ECHO "  running pw.x as: $PW_COMMAND"
$ECHO "  running turbo_lanczos.x as: $TURBO_LANCZOS_COMMAND"
$ECHO "  running turbo_spectrum.x as: $TURBO_SPECTRUM_COMMAND"
$ECHO


# Input for Environ
if [ "${environ_type}" = "input" ] ; then
cat > environ.in << EOF
 &ENVIRON
   !
   verbose = $verbose
   environ_thr = $environ_thr
   environ_type = '$environ_type'
   env_static_permittivity = $env_static_permittivity
   env_optical_permittivity = $env_optical_permittivity
   env_surface_tension = $env_surface_tension
   env_pressure = $env_pressure
   !
 /
 &BOUNDARY
 /
 &ELECTROSTATIC
   !
   tol = $tol
   !
 /
EOF
fi
 
prefix=CH4_$prefix_environ

# Ground-state SCF calculation
  cat > ${prefix}.scf.in << EOF
 &CONTROL
   !
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR/'
   outdir = '$TMP_DIR/'
   prefix = '$prefix'
   verbosity = 'high'
   !
 /
 &SYSTEM
   !
   ibrav = 1
   celldm(1) = 20
   nat = 5
   ntyp = 2
   ecutwfc = 25
   !
/
 &ELECTRONS
   !
   diagonalization = 'davidson'
   mixing_mode = 'plain'
   mixing_beta = 0.7
   conv_thr = 1.D-8
   electron_maxstep = 100
   !
 /
ATOMIC_SPECIES  
 H   1  H.pz-vbc.UPF
 C  12  C.pz-vbc.UPF
ATOMIC_POSITIONS {Angstrom}
C        0.000000000   0.000000000   0.000000000
H        0.642814093   0.642814093   0.642814093
H       -0.642814093  -0.642814093   0.642814093
H        0.642814093  -0.642814093  -0.642814093
H       -0.642814093   0.642814093  -0.642814093
K_POINTS {gamma}
EOF
$ECHO "  running the scf calculation...\c"
$PW_COMMAND < ${prefix}.scf.in > ${prefix}.scf.out
check_failure $?
$ECHO " done"



# TDDFPT calculation using the Lanczos algorithm
cat > ${prefix}.turbo-lanczos.in << EOF
&lr_input
    prefix='$prefix'
    outdir='$TMP_DIR/'
    restart_step = 250
    restart = .false.
/
&lr_control
    itermax = 3000
    ipol = 4
    n_ipol = 3
/
EOF
$ECHO "  running the TDDFPT calculation using the Lanczos algorithm...\c"
$TURBO_LANCZOS_COMMAND < ${prefix}.turbo-lanczos.in > ${prefix}.turbo-lanczos.out
check_failure $?
$ECHO " done"



# Post-processing spectrum calculation 
cat > ${prefix}.turbo-spectrum.in << EOF
&lr_input
  prefix = '$prefix',
  outdir = '$TMP_DIR/',
  itermax0 = 3000,
  itermax =  20000,
  extrapolation="osc"
  epsil = 0.01
  start = 0.0d0
  increment = 0.001d0
  end = 3.5d0
  ipol = 4
/
EOF
$ECHO "  running the post-processing spectrum calculation...\c"
$TURBO_SPECTRUM_COMMAND < ${prefix}.turbo-spectrum.in > ${prefix}.turbo-spectrum.out
check_failure $?
$ECHO " done"


done


$ECHO
$ECHO "$EXAMPLE_DIR : done"
