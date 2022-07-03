#!/bin/sh

# run from directory where this script is
cd $(echo $0 | sed 's/\(.*\)\/.*/\1/') # extract pathname
EXAMPLE_DIR=$(pwd)

# check whether echo has the -e option
if test "$(echo -e)" = "-e"; then ECHO=echo; else ECHO="echo -e"; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example tests K-edge X-ray absorption spectra calculation"
$ECHO

# set the needed environment variables
x=$EXAMPLE_DIR
while test "$x" != "/"; do
   x=$(dirname "$x")
   if test -f "$x/environment"; then
      . "$x/environment"
      break
   fi
done

# required executables and pseudopotentials
BIN_LIST="pw.x xspectra.x "
PSEUDO_LIST="Si_PBE_USPP.UPF O_PBE_USPP.UPF"
PSEUDO_DIR="$EXAMPLE_DIR/pseudo"
XSPECTRA_DIR="$PREFIX/XSpectra"
PLOT="$XSPECTRA_DIR/tools/upf2plotcore.sh"

$ECHO
$ECHO "  executables directory: $BIN_DIR"
$ECHO "  pseudo directory:      $PSEUDO_DIR"
$ECHO "  temporary directory:   $TMP_DIR"
$ECHO
$ECHO "  checking that needed directories and files exist...\c"

# check for directories
for DIR in "$BIN_DIR" "$PSEUDO_DIR"; do
   if test ! -d $DIR; then
      $ECHO
      $ECHO "ERROR: $DIR not existent or not a directory"
      $ECHO "Aborting"
      exit 1
   fi
done
for DIR in "$TMP_DIR" "$EXAMPLE_DIR/results"; do
   if test ! -d $DIR; then
      mkdir $DIR
   fi
done
cd $EXAMPLE_DIR/results

# check for executables
for FILE in $BIN_LIST; do
   if test ! -x $BIN_DIR/$FILE; then
      $ECHO
      $ECHO "ERROR: $BIN_DIR/$FILE not existent or not executable"
      $ECHO "Aborting"
      exit 1
   fi
done

# check for pseudopotentials
for FILE in $PSEUDO_LIST; do
   if test ! -r $PSEUDO_DIR/$FILE; then
      $ECHO
      $ECHO "ERROR: $PSEUDO_DIR/$FILE not existent or not readable"
      $ECHO "Aborting"
      exit 1
   fi
done
$ECHO " done"

# how to run executables
PW_COMMAND="$PARA_PREFIX $BIN_DIR/pw.x --environ $PARA_POSTFIX"
X_COMMAND="$PARA_PREFIX $BIN_DIR/xspectra.x --environ $PARA_POSTFIX"

$ECHO
$ECHO "  running pw.x as: $PW_COMMAND"
$ECHO "  running xspectra.x as: $X_COMMAND"

$ECHO
$ECHO "WARNING : All these calculations are underconverged"
$ECHO "    (These are simple quick tests)  "
$ECHO

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/*
$ECHO " done"
# extracting core wavefunction
$ECHO "  extracting core wavefunction from pseudo...\c"
$PLOT $PSEUDO_DIR/Si_PBE_USPP.UPF >./Si.wfc
$ECHO " done"

$ECHO
$ECHO "  self-consistent calculation for SiO2"
$ECHO "  Ultrasoft pseudopotentials are being used "
$ECHO

cat >SiO2.scf.in <<EOF
 &control
    calculation='scf',
    restart_mode='from_scratch',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/',
    prefix='SiO2',
 &end
 &system
    ibrav = 4 ,
    celldm(1)=9.28630318961333,
    celldm(3)=2.20020,
    nat = 9 ,
    ntyp = 3 ,
    nspin=1,
    ecutwfc = 20.0,
    ecutrho = 150.0,
    nbnd=30,
    smearing='mp',
    occupations='smearing',
    degauss=0.03,
 &end
 &electrons
    diagonalization='david',
    conv_thr = 1.d-9,
    mixing_mode = 'plain',
    mixing_beta = 0.3,
 &end
ATOMIC_SPECIES
Sih   28.086 Si_PBE_USPP.UPF
Si   28.086 Si_PBE_USPP.UPF
O  15.9994   O_PBE_USPP.UPF
ATOMIC_POSITIONS bohr
Sih  4.364561 0 0
Si  -2.1822805 3.779820333 6.81057242
Si  -2.1822805 -3.779820333 3.40528621
O    2.593199275 2.152889361 1.214665684
O    0.567857245 3.322221185 5.595906736
O    -3.16105652 1.169331825 8.025238104
O    -3.16105652 -1.169331825 2.190620526
O    0.567857245 -3.322221185 4.619951894
O    2.593199275 -2.152889361 -1.214665684
K_POINTS automatic
2 2 1 0 0 0
EOF
cat >environ.in <<EOF
 &ENVIRON
   !
   verbose = 1
   environ_thr = 1.d2
   env_static_permittivity = 100. 
   !
 /
 &BOUNDARY
  solvent_mode = 'ionic'
 /
 &ELECTROSTATIC
 /
EOF

$ECHO "  running pw.x for SiO2 ...\c"
$PW_COMMAND <SiO2.scf.in >SiO2.scf.out
check_failure $?
$ECHO " done"

$ECHO
$ECHO "  xanes calculations : Si K-edge"
$ECHO "  ultrasoft pseudopotentials are being used"
$ECHO "  cutting occupied states "
$ECHO
$ECHO "  x-ray absorption spectrum calculation"
$ECHO "  dipolar part"
$ECHO
$ECHO "  epsilon in the plane direction"
$ECHO

cat >environ.in <<EOF
 &ENVIRON
   !
   verbose = 1
   environ_restart = .true.
   environ_thr = 1.d2
   env_static_permittivity = 100. 
   !
 /
 &BOUNDARY
  solvent_mode = 'ionic'
 /
 &ELECTROSTATIC
 /
EOF

cat >SiO2.xspectra_dip_plane.in <<EOF
 &input_xspectra
    calculation='xanes_dipole'
    prefix='SiO2',
    outdir='$TMP_DIR/',
    xonly_plot=.false.,
    xniter=2000,
    xcheck_conv=50,
    xepsilon(1)=1.0,
    xepsilon(2)=1.0,
    xepsilon(3)=0.0,
    xiabs=1,
    x_save_file='SiO2.xspectra_dip_plane.sav',
    xerror=0.001,
 / 
 &plot
    xnepoint=1000,
    xgamma=0.8,
    xemin=-10.0,
    xemax=100.0,
    terminator=.true.,
    cut_occ_states=.true.,
 &end
 &pseudos
    filecore='Si.wfc',
    r_paw(1)=2.4,
 &end
 &cut_occ
    cut_desmooth=0.1,
 &end
2 2 1 0 0 0
EOF

$ECHO "  running xspectra.x ...\c"
$X_COMMAND <SiO2.xspectra_dip_plane.in >SiO2.xspectra_dip_plane.out
check_failure $?
mv xanes.dat SiO2.xspectra_dip_plane.dat
$ECHO " done"

$ECHO
$ECHO "  epsilon in the c axis direction"
$ECHO

cat >SiO2.xspectra_dip_c.in <<EOF
 &input_xspectra
    calculation='xanes_dipole'
    prefix='SiO2',
    outdir='$TMP_DIR/',
    xonly_plot=.false.,
    xniter=2000,
    xcheck_conv=50,
    xepsilon(1)=0.0,
    xepsilon(2)=0.0,
    xepsilon(3)=1.0,
    xiabs=1,
    x_save_file='SiO2.xspectra_dip_c.sav',
    xerror=0.001,
 /
 &plot
    xnepoint=1000,
    xgamma=0.8,
    xemin=-10.0,
    xemax=100.0,
    terminator=.true.,
    cut_occ_states=.true.,
 &end
 &pseudos
    filecore='Si.wfc',
    r_paw(1)=2.4,
 &end
 &cut_occ
    cut_desmooth=0.1,
 &end
2 2 1 0 0 0
EOF

$ECHO "  running xspectra.x ...\c"
$X_COMMAND <SiO2.xspectra_dip_c.in >SiO2.xspectra_dip_c.out
check_failure $?
mv xanes.dat SiO2.xspectra_dip_c.dat
$ECHO " done"
