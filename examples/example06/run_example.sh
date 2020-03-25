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
$ECHO "system immersed in a continuum electrolyte solution with open boundary "
$ECHO "conditions along the axis perpendicular to the slab plane, as described in "
$ECHO
$ECHO "     O. Andreussi and N. Marzari, Phys. Rev. B 90, 245101 (2014). " 
$ECHO
$ECHO "The electrochemical diffuse layer is modeled by solving different flavors "
$ECHO "of the Poisson-Boltzmann equation, see "
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
$ECHO "          F. Nattino, M. Truscott, N. Marzari, O. Andreussi,"
$ECHO "                   J. Chem. Phys. 150, 041722 (2019)."
$ECHO

# set the needed environment variables
. ../../../environment_variables

# compatibility with QE for versions prior to 6.4
if [ -z $NETWORK_PSEUDO ]; then
    NETWORK_PSEUDO=http://www.quantum-espresso.org/wp-content/uploads/upf_files/
fi

# required executables and pseudopotentials
BIN_LIST="pw.x"
PSEUDO_LIST="Ag.pbe-n-rrkjus_psl.1.0.0.UPF"

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
### ELECTROLYTE PARAMETERS ##########################################
env_electrolyte_ntyp=2      # number of electrolyte species in solution.
zion=1                      # ionic charge.
cion='1.0'                  # ionic concentration (in mol/l).
temperature=300             # temperature of the electrolyte solution.
electrolyte_mode='system'   # specify the method that is used to build
                            # the electrolyte cavity:
                            # electronic, ionic, full: same as for solvent
                            #   cavity, see previous examples. the 
                            #   corresponding parameters are 
                            #   electrolyte_rhomax, electrolyte_rhomin 
                            #   (for electronic and full), or
                            #   electrolyte_alpha and electrolyte_softness 
                            #   (for ionic);
                            # system: simplified erf-based interface. 
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
electrolyte_distance='30.0' # distance parameter for the simplified
                            # erf-based interface.
electrolyte_spread='0.01'   # spread parameter for the simplified
                            # erf-based interface.
### PERIODIC BOUNDARY CONDITIONS (PBC) ##############################
pbc_dim=2                  # select the desired system dimensionality
                           # 0, 1 or 2: isolated, 1D or 2D system.
pbc_axis=3                 # set the axis along the 1D direction or 
                           # normal to the 2D plane (pbc_axis = 1, 2 
                           # or 3 for x, y or z axis).
#####################################################################

prefix='scf'
input=${prefix}'.in'
output=${prefix}'.out'
rm -f results.txt

cat > $input << EOF 
&CONTROL
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '$PSEUDO_DIR/'
   outdir = '$TMP_DIR/'
   prefix = '$prefix'
   tprnfor = .true.
   verbosity = 'high'
/
&SYSTEM
   ntyp             = 1
   nat              = 8
   ibrav            = 0
   ecutwfc = 30
   ecutrho = 300
   occupations = 'smearing'
   degauss = 0.02
   smearing = 'mv'
   tot_charge = 0.5
/
&ELECTRONS
   conv_thr = 5.D-9
   diagonalization = 'davidson'
!   mixing_beta = 0.1
/
&IONS
/
&CELL
/

ATOMIC_SPECIES
Ag 107.8682 Ag.pbe-n-rrkjus_psl.1.0.0.UPF

K_POINTS automatic
2 2 1 0 0 0

CELL_PARAMETERS angstrom
2.93350761311527 0.00000000000000 0.00000000000000
0.00000000000000 2.93350761311527 0.00000000000000
0.00000000000000 0.00000000000000 54.44175455526538

ATOMIC_POSITIONS angstrom
Ag 1.4667538066 1.4667538066 20.0000000000
Ag 0.0000000000 0.0000000000 22.0336703459
Ag 1.4667538066 1.4667538066 24.1093139696
Ag 0.0000000000 0.0000000000 26.1848761000
Ag 1.4667538066 1.4667538066 28.2568784553
Ag 0.0000000000 0.0000000000 30.3324405857
Ag 1.4667538066 1.4667538066 32.4080836802
Ag 0.0000000000 0.0000000000 34.4417545553
EOF

for model in PB MPB ; do
  # two electrolyte models are tested:
  # PB  : standard Poisson-Boltzmann model, with point-like ions;
  # MPB : size-modified PB model, with finite-sized ions.
  #   the MPB model is activated through one of the following params:
  #   - cionmax : parameter setting the maximum electrolyte concentration
  #     (in mol/l); cionmax = 1/a**3, where a is the ionic size parameter
  #     originally employed by Borukhov et al.
  #   - rion : parameter corresponding to the ionic radius (in Bohr); it
  #     is related to cionmax through: cionmax = p / (4/3*pi*rion**3),
  #     where p=0.64 is the volume fraction for random close packing.

  if [ $model == "PB" ] ; then
      # standard PB model;
      # for high ionic strengths and low surface charge the linear-regime
      # approximation can be employed.
      # available implementations: numerical and (2D)-analytic;
      model="PB" ; cionmax=0
      linearized="true false"
      pbc_corrections="gcs parabolic"
  else
      # MPB model;
      # available implementations: numerical only;
      model="MPB" ; cionmax=10.0
      pbc_corrections="parabolic"
      linearized="false"
  fi

  for linear in $linearized ; do
    # .true. : consider the linear-regime approximation of the PB model;
    # .false. : solve the full non-linear PB equation.

    if [ $linear == "true" ]; then
       $ECHO "$model model (linearized) " >> results.txt
    else 
       $ECHO "$model model " >> results.txt
    fi

    for pbc_correction in $pbc_corrections ; do
      # "parabolic" : real-space parabolic PBC correction; the (M)PB
      #     equation is solved numerically using an iterative approach;
      # "gcs" : Gouy-Chapman-Stern (GCS) correction; the PB equation is
      #     solved analytically and imposed as boundary condition.
      #     NOTE: only available for 2D systems!

      # clean TMP_DIR
      $ECHO "  cleaning $TMP_DIR...\c"
      rm -rf $TMP_DIR/*
      $ECHO " done"

      $ECHO "  running the scf calculation with a diffuse ionic charge density;"
      $ECHO "  the $model equation is solved with the $pbc_correction PBC correction"
      $ECHO "  (linear-regime approximation: $linear);"

      label="Ag100_model-${model}_pbccorr-${pbc_correction}_linear-${linear}"
      cat > environ.in << EOF
&ENVIRON
   verbose = $verbose
   environ_thr = $environ_thr
   environ_type = 'input'
   env_static_permittivity = 80
   env_surface_tension = 0.D0
   env_pressure = 0.D0
   env_electrolyte_ntyp = $env_electrolyte_ntyp
   zion(1) = $zion
   zion(2) = -$zion
   cion(1) = $cion
   cion(2) = $cion
   cionmax = $cionmax
   temperature = $temperature
   system_dim = $system_dim            
   system_axis = $system_axis   
   system_ntyp = $system_ntyp
   electrolyte_linearized = .$linear.
/
&BOUNDARY
   solvent_mode = 'full'
   electrolyte_mode = '$electrolyte_mode'
   electrolyte_distance = $electrolyte_distance
   electrolyte_spread = $electrolyte_spread   
/
&ELECTROSTATIC
   pbc_correction = '$pbc_correction'
   pbc_dim = $pbc_dim
   pbc_axis = $pbc_axis
   tol = 5.D-13
   inner_tol = 5.D-18
/
EOF
      $PW_COMMAND < $input > $label.out
      check_failure $?
      mv environ.in $label.environ.in
      $ECHO " done"

      $ECHO
      $ECHO "$EXAMPLE_DIR : done"
    
      fermi=$(awk '/the Fermi energy is/ {en=$5}; END {print en}' $label.out )
      potshift=$(awk '/the potential shift due/ {en=$10}; END {print en}' $label.out )
      pot=$($ECHO "-($fermi+$potshift)" | bc -l)

      if [ $pbc_correction == "parabolic" ]; then
         method="numerical solution   :"
      else :
         method="2D-analytic solution :"
      fi
      $ECHO "$method potential = -($fermi + $potshift) V = $pot V (abs. scale)" >> results.txt
    done
    $ECHO >> results.txt
  done
done
