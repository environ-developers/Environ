#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use pw.x to simulate a 10xH20 water cluster "
$ECHO "immersed in continuum solvent. Turning on the solvent_radius inside the "
$ECHO "BOUNDARY namelist prevents dielectric from getting inside the cluster."
$ECHO

# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x"
PSEUDO_LIST="H.pbe-rrkjus.UPF O.pbe-rrkjus.UPF"

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
verbose=0                      # if GE 1 prints debug informations
                               # if GE 2 prints out gaussian cube files with 
                               # dielectric function, polarization charges, etc
                               # WARNING: if GE 2 lot of I/O, much slower
environ_thr='1.d0'             # electronic convergence threshold for the onset  
                               # of solvation correction
environ_type="water"           # type of environment
                               # input: read parameters from input
                               # vacuum: all flags off, no environ 
### BOUNDARY PARAMETERS #####################################
filling_threshold=8.25d-01     # if filling function GE threshold, 
                               # the dielectric is set to 1, so dielectric pockets
                               # of the solvent size will be prevented
                               # to prevent cylindrical pockets: decrease to approx. 0.65.
                               # For below 0.5, the dielectric boundary might be modified strongly
                               # will be treated as quantum mechanical
solvent_radius=3               # size of solvent molecules (radius(H20) approx 3 bohr)
### ITERATIVE SOLVER PARAMETERS #####################################
tol='1.d-12'  # tolerance of the iterative solver
#####################################################################

for solvent_awareness in standard solvent_aware ; do 

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/*
$ECHO " done"

if [ "${solvent_awareness}" = "standard" ] ; then
   label="standard"
else
   label="fill_pockets"
fi

prefix=H2OCluster_$label
input=${prefix}'.in'
output=${prefix}'.out'

# Input for Environ standard
if [ "${solvent_awareness}" = "standard" ] ; then
cat > environ.in << EOF
 &ENVIRON
   !
   verbose = $verbose
   environ_thr = $environ_thr
   environ_type = '$environ_type'
   !
 /
 &BOUNDARY
   solvent_mode = 'electronic'
 /
 &ELECTROSTATIC
   !
   tol = $tol
 /
EOF
fi

if [ "${solvent_awareness}" = "solvent_aware" ] ; then
cat > environ.in << EOF
 &ENVIRON
   !
   verbose = $verbose
   environ_thr = $environ_thr
   environ_type = '$environ_type'
   !
 /
 &BOUNDARY
   solvent_mode = 'electronic'
   filling_threshold = $filling_threshold
   solvent_radius = $solvent_radius
 /
 &ELECTROSTATIC
   !
   tol = $tol
   !
 /
EOF
fi

cp environ.in environ_${prefix}.in   


 

# Ground-state SCF calculations
  cat > $input << EOF
 &CONTROL
   !
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR/'
   outdir = '$TMP_DIR/'
   prefix = '$prefix'
   verbosity = 'high'
   tprnfor = .true.
   !
 /
 &SYSTEM
   !
   assume_isolated = 'mt'
   ecutrho = 280
   ecutwfc = 35
   ibrav = 0
   nat = 30
   ntyp = 2
   occupations = 'fixed'
   !
/
 &ELECTRONS
   !
   conv_thr =   2.0d-08
   electron_maxstep = 100
   mixing_beta =   3.0d-01
   !
 /
ATOMIC_SPECIES
H      1.00794 H.pbe-rrkjus.UPF 
O      15.9994 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS angstrom
O            9.7407800000       9.5971600000       6.5018633333 
H           10.2290800000      10.1831600000       5.9108633333 
H            8.8803800000      10.0472600000       6.6731633333 
O            9.3199800000       6.8106600000       6.0826633333 
H            9.4998800000       7.7702600000       6.1095633333 
H            9.6905800000       6.4740600000       6.9229633333 
O            7.3101800000      10.5947600000       7.2777633333 
H            6.8562800000      11.3708600000       6.9263633333 
H            6.6440800000       9.8498600000       7.2166633333 
O            6.7248800000       6.2260600000       6.3047633333 
H            7.6774800000       6.4711600000       6.1439633333 
H            6.4600800000       5.6954600000       5.5435633333 
O            9.8923800000       5.9950600000       8.7285633333 
H            8.9392800000       5.7891600000       8.9336633333 
H           10.4030800000       5.2569600000       9.0824633333 
O            7.2755800000       5.5673600000       9.0328633333 
H            6.8204800000       6.3086600000       9.4818633333 
H            6.9531800000       5.6153600000       8.1129633333 
O            7.8983800000       9.8013600000      10.0250633333 
H            7.7074800000      10.3022600000       9.2121633333 
H            8.7987800000       9.4352600000       9.8566633333 
O            6.0458800000       7.9070600000       9.9967633333 
H            5.5153800000       7.9994600000      10.7974633333 
H            6.7709800000       8.5922600000      10.0746633333 
O            5.5767800000       8.5916600000       7.2964633333 
H            5.8959800000       7.7908600000       6.8329633333 
H            5.5353800000       8.3223600000       8.2357633333 
O           10.3283800000       8.7649600000       9.2297633333 
H           10.2693800000       7.7933600000       9.1483633333 
H           10.3395800000       9.0758600000       8.3073633333 
K_POINTS automatic
1 1 1 0 0 0
CELL_PARAMETERS angstrom
     14.0000000000       0.0000000000       0.0000000000
      0.0000000000      14.0000000000       0.0000000000
      0.0000000000       0.0000000000      14.0000000000      
EOF

$PW_COMMAND < $input > $output
check_failure $?
rm environ.in

$ECHO " done" ${prefix}

done
