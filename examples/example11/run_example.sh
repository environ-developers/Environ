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
$ECHO "system immersed in a continuum semiconductor  with open boundary "
$ECHO "conditions along the axis perpendicular to the slab plane, as described in "
$ECHO
$ECHO "     XXX " 
$ECHO
$ECHO "The electrochemical diffuse layer is modeled by solving the Poisson-Boltzmann "
$ECHO "equation under Mott-Schottky conditions. See "
$ECHO
$ECHO "  I. Borukhov, D. Andelman, and H. Orland, Phys. Rev. Lett. 79, 435 (1997);"
$ECHO "  I. Dabo, E. Cancès, Y.L. Li, and N. Marzari, arXiv 0901.0096 (2008); "
$ECHO "      S. Ringe, H. Oberhofer, C. Hille, S. Matera, and K. Reuter, "
$ECHO "                J. Comput. Theor. Chem. 12, 4052 (2016);"
$ECHO "    G. Fisicaro, L. Genovese, O. Andreussi, N. Marzari and S. Goedecker,"
$ECHO "                   J. Chem. Phys. 144, 014103 (2016). "
$ECHO
$ECHO "Details on the implementation are described in"
$ECHO
$ECHO "          XXXX  et al. "
$ECHO "             XXX (2022)"
$ECHO

# set the needed environment variables
. ../../../environment_variables

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
environ_thr='1.d2'    # electronic convergence threshold for the onset
                      # of solvation correction.
sc_permittivity='11.0' # The dielectric permitivitty of the surrounding 
                      # semiconductor

### ELECTROLYTE PARAMETERS ##########################################
temperature=300             # temperature of the surrounding semiconductor
system_dim=2                # select dimensionality of the interface
                            # 0, 1 and 2: spherical, cylindrical or
                            #   planar interface.
system_axis=3               # select axis for the simplified interface.
                            # the axis is along the 1D direction or 
                            # normal to the 2D plane (1, 2 or 3 for 
                            # x, y or z axis).
system_ntyp=1               # the geometrical center of the atoms of 
                            # all types up to system_ntyp defines the
                            # origin of the simplified interface.
sc_distance='15.0'          # distance parameter for the simplified
                            # erf-based interface of where semiconductor
                            # is embedded
sc_spread='0.01'            # spread parameter for the simplified
                            # erf-based interface.
### PERIODIC BOUNDARY CONDITIONS (PBC) ##############################
pbc_correction='ms'        # selecting the mott-schottky pbc correction
                           # for semiconductor embedded calculations
pbc_dim=2                  # select the desired system dimensionality
                           # 0, 1 or 2: isolated, 1D or 2D system.
pbc_axis=3                 # set the axis along the 1D direction or 
                           # normal to the 2D plane (pbc_axis = 1, 2 
                           # or 3 for x, y or z axis).
#####################################################################

prefix='scf'
input=${prefix}'.in'
rm -f results.txt

cat > $input << EOF 
&CONTROL
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR'
   outdir = '$TMP_DIR/'
   prefix = '$prefix'
   tprnfor = .true.
   verbosity = 'high'
/
&SYSTEM
  degauss = 0.0001,
  ecutwfc = 20,
  occupations = 'smearing',
  smearing = 'mv',
  ibrav = 0,
  nat = 4,
  ntyp = 1,
  tot_charge= 0.001, ! This seems like a very small amount of charge (and it is)
                     ! but recall that due to periodic boundary conditions, this
                     ! represents a fairly large concentration of dopants, i.e.
                     ! equivalent to 1.63e19 e/cm^-3 
/
&ELECTRONS
  conv_thr = 1d-08,
  diagonalization = 'david',
  mixing_mode = 'local-TF',
/
&IONS
/
&CELL
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
  -0.000000 3.866975 0.000000
  0.000000 0.000000 32.812368
EOF

for sc_carrier_density in '1.D16' '1.D17' '1.D18' ; do
  # testing three different carrier densities within the surrounding
  # semiconductor region, demonstrating the importance of concentrations
  # in sheilding semiconductor potentials

  $ECHO "  cleaning $TMP_DIR...\c"
  rm -rf $TMP_DIR/*
  $ECHO " done"

  cat > environ.in << EOF
&ENVIRON
   verbose = $verbose
   environ_thr = $environ_thr
   temperature = $temperature
   system_dim = $system_dim            
   system_axis = $system_axis   
   system_ntyp = $system_ntyp
   sc_carrier_density = $sc_carrier_density
/
&BOUNDARY
   sc_distance = $sc_distance
   sc_spread = $sc_spread   
/
&ELECTROSTATIC
   pbc_correction = '$pbc_correction'
   pbc_dim = $pbc_dim
   pbc_axis = $pbc_axis
   tol = 5.D-13
   inner_tol = 5.D-18
/
EOF

  sc_carrier_density=$(echo $sc_carrier_density | sed s/D//)

  $ECHO "  running the scf calculation with a diffuse ionic charge density;"
  $ECHO "  the PB equation is solved with the $pbc_correction boundary conditions"
  $ECHO "  using a carrier density of $sc_carrier_density cm^-3 for semiconductor region"

  label="Si-ms-carrier-density-${sc_carrier_density}"
  $PW_COMMAND < $input > $label.out
  check_failure $?
  mv environ.in $label.environ.in
  $ECHO "  done"
  $ECHO

  fermi=$(awk '/the Fermi energy is/ {en=$5}; END {print en}' $label.out )
  $ECHO "$sc_carrier_density cm^-3 potential = $fermi V (abs. scale)" >> results.txt
done
$ECHO >> results.txt
$ECHO
$ECHO "$EXAMPLE_DIR : done"
