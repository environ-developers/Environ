#!/bin/sh

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use pw.x to simulate a charged semiconductor- "
$ECHO "solution system, determining the voltage drop of the ultimate system "
$ECHO "compared to the flatband potential as described in "
$ECHO
$ECHO " Q. Campbell, D. Fisher, I. Dabo, Phys. Rev. Mater. 3 (1), 015404 (2019). " 
$ECHO
$ECHO "The electrochemical diffuse layer is modeled by solving parabolic boundary "
$ECHO "conditions with two helmholtz planes in the system representing the bulk  "
$ECHO " solution and semiconductor charges acting on the system. See: "
$ECHO
$ECHO " Q. Campbell and I. Dabo, Phys. Rev. B 95 (20), 205308 (2017). "
$ECHO " O. Andreussi and N. Marzari, Phys. Rev. B 90, 245101 (2014)"
$ECHO " O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, 064102 (2012)"
$ECHO
$ECHO

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
PSEUDO_LIST="Si.pbe-n-kjpaw_psl.1.0.0.UPF"

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

### ELECTROSTATIC EMBEDDING PARAMETERS ##############################
verbose=0             # if GE 1 prints debug informations;
                      # if GE 2 prints out gaussian cube files with 
                      # electrostatic potential, electrolyte charge, etc.
                      # WARNING: if GT 2 lot of I/O, much slower.
sc_permittivity='11.0' # The dielectric permitivitty of the surrounding 
                      # semiconductor
sc_carrier_density='1.D18' #The background carrier density of the 
                      # semiconductor in units of cm^-3
sc_electrode_chg='0.001' # The charge on the overall semiconductor 
                      # electrode. In units of (e). Positive charge indicates
                      # missing electrons, i.e. following QE charge convention
sc_chg_thr='1.D-5' # the threshold below which the dft charge  
                      # needs to change in the outer loop before 
                      # convergence
env_static_permittivity='78.3' # the dielectric permitivitty of the 
                      # solution in the system

### ELECTROLYTE PARAMETERS ##########################################
temperature=300             # temperature of the surrounding semiconductor
system_dim=2                # select dimensionality of the interface
                            # 0, 1 and 2: spherical, cylindrical or
                            #   planar interface.
system_axis=3               # select axis for the simplified interface.
                            # the axis is along the 1D direction or 
                            # normal to the 2D plane (1, 2 or 3 for 
                            # x, y or z axis).
sc_distance='3.0'          # distance from center of system where 
                            # cutoff between bulk and DFT semiconductor
                            # takes place
sc_spread='7.0'            # spread parameter of macroscopic averaging
                            # of planar potential for determining difference
                            # between charged and neutral interface potentials
### PERIODIC BOUNDARY CONDITIONS (PBC) ##############################
pbc_correction='ms-gcs'    # selecting the mott-schottky-guoy-chapman-stern
                           # pbc correction for extended semiconductor- 
                           # solution interfaces
pbc_dim=2                  # select the desired system dimensionality
                           # 0, 1 or 2: isolated, 1D or 2D system.
pbc_axis=3                 # set the axis along the 1D direction or 
                           # normal to the 2D plane (pbc_axis = 1, 2 
                           # or 3 for x, y or z axis).
#####################################################################

prefix='scf'
input=${prefix}'.in'

cat > $input << EOF 
&CONTROL
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR'
   outdir = '$TMP_DIR/'
   prefix = '$prefix'
/
&SYSTEM
  degauss = 0.0001,
  ecutwfc = 20,
  occupations = 'smearing',
  smearing = 'mv',
  ibrav = 0,
  nat = 4,
  ntyp = 1,
/
&ELECTRONS
  conv_thr = 1d-12, ! high convergence thresholds needed
                    ! for reliable voltage results
  diagonalization = 'david',
/
ATOMIC_SPECIES
  Si  28.0855 Si.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si 1.000000 0.500000 0.562500
  Si 0.500000 0.500000 0.437500
  Si 0.500000 1.000000 0.479167
  Si 1.000000 0.000000 0.520833
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  3.866975 0.000000 0.000000
  0.000000 3.866975 0.000000
  0.000000 0.000000 32.812368
EOF


cat > environ.in << EOF
&ENVIRON
   verbose = $verbose
   temperature = $temperature
   system_dim = $system_dim            
   system_axis = $system_axis   
   env_static_permittivity = $env_static_permittivity 
   env_electrolyte_ntyp = 2,
   electrolyte_linearized = .false., ! providing details for gcs
   zion(1) = 1,    ! portion of calculation 
   cion(1) = 0.01,
   zion(2) = -1,
   cion(2) = 0.01 ,
   sc_permittivity = $sc_permittivity
   sc_carrier_density = $sc_carrier_density
   sc_electrode_chg = $sc_electrode_chg
   sc_chg_thr = $sc_chg_thr
/
&BOUNDARY
   sc_distance = $sc_distance
   sc_spread = $sc_spread   
   electrolyte_spread = 0.01,
   electrolyte_distance = 15.2, ! in bohr
   electrolyte_mode = 'system',
/
&ELECTROSTATIC
   solver = 'fixed-point'
   auxiliary = 'full'
   pbc_correction = '$pbc_correction'
   pbc_dim = $pbc_dim
   pbc_axis = $pbc_axis
   tol = 1.D-11
/
EOF


$ECHO "  running an ms-gcs calculation to determine voltage associated with electrode"
$ECHO "  charge. Throughout, the calculation will be trying to determine the correct"
$ECHO "  charge distribution between the surface slab and the extended bulk semiconductor "
$ECHO "  with a carrier concentration of $sc_carrier_density cm^-3."
$ECHO "  In practice, runs similar to a ionic relaxation, with each scf step going down "
$ECHO "  a gradient to match the Fermi levels of the DFT and the bulk semiconductor. "
$ECHO

label="Si-ms-gcs-${sc_electrode_chg}"
$PW_COMMAND < $input > $label.out
check_failure $?
mv environ.in $label.environ.in
$ECHO "  done"
$ECHO

$ECHO
$ECHO "$EXAMPLE_DIR : done"
