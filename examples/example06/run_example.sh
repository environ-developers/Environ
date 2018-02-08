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
$ECHO "system immersed in continuum solvent described by the SCCS model "
$ECHO
$ECHO "   O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, 064102 (2012) "
$ECHO
$ECHO "and with open boundary conditions along the axis perpendicular "
$ECHO "to the slab plane, as described in "
$ECHO
$ECHO "   O. Andreussi and N. Marzari, Phys. Rev. B 90, 245101 (2014) " 
$ECHO
$ECHO "The presence of electrolytes in the solution is accounted for through"
$ECHO "the linearized modified Poisson-Boltzmann model, see "
$ECHO 
$ECHO " I. Borukhov, D. Andelman, and H. Orland, Phys. Rev. Lett. 79, 435 (1997)"
$ECHO "   I. Dabo, E. Cancès, Y.L. Li, and N. Marzari, arXiv 0901.0096 (2008) "
$ECHO "       S. Ringe, H. Oberhofer, C. Hille, S. Matera, and K. Reuter, "
$ECHO "            J. Comput. Theor. Chem. 12, 4052 (2016)."
$ECHO "   G. Fisicaro, L. Genovese, O. Andreussi, N. Marzari and S. Goedecker,"
$ECHO "                  J. Chem. Phys. 144, 014103 (2016) "

# set the needed environment variables
. ../../../environment_variables

# required executables and pseudopotentials
BIN_LIST="pw.x"
PSEUDO_LIST="O.pbe-rrkjus.UPF C.pbe-rrkjus.UPF Pt.pbe-nd-rrkjus.UPF"

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
environ_thr='1.d0'    # electronic convergence threshold for the onset  
                      # of solvation correction
### ELECTROLYTE PARAMETERS ##########################################
env_electrolyte_ntyp=2  # number of electrolyte species in solution
zion=1                  # ionic charge
cion='5.0'              # ionic concentration (in mol/l)
cionmax='15.0'          # max ionic concentration (in mol/l) 
                        # NOTE: cionmax=1/a**3, where a is the ionic 
                        # size parameter originally employed by Borukhov. 
                        # one can alternatively set the ionic radius
                        # rion (in Bohr), related to cionmax according
                        # to: cionmax=p/(4/3*pi*rion**3), where p=0.64 
                        # is the random close packing volume fraction.
solvent_temperature=300 # temperature of the electrolyte solution 
stern_mode='system'     # specify the method that is used to build
                        # the electrolyte cavity
                        # electronic, ionic, full: same as for solvent
                        #   cavity, see previous examples. the 
                        #   corresponding parameters are stern_rhomax
                        #   and stern_rhomin (for electronic and full)
                        #   or stern_alpha and stern_softness (for 
                        #   ionic).
                        # system: simplified erf-based interface 
system_dim=2            # select dimensionality of the interface
                        # 0, 1 and 2: spherical, cylindrical or
                        #   planar interface
system_axis=3           # select axis for the simplified interface.
                        # the axis is along the 1D direction or 
                        # normal to the 2D plane (1, 2 or 3 for 
                        # x, y or z axis)
system_ntyp=1           # the geometrical center of the atoms of 
                        # all types up to system_ntyp defines the
                        # origin of the simplified interface
stern_distance='10.0'   # distance parameter for the simplified
                        # erf-based interface
stern_spread='0.5'      # spread parameter for the simplified
                        # erf-based interface 
### PERIODIC BOUNDARY CONDITIONS ####################################
pbc_correction='parabolic' # correction scheme to remove PBC 
                           # none: periodic calculation, no correction 
                           # parabolic: quadratic real-space correction
pbc_dim=2                  # select the desired system dimensionality
                           # 0, 1 or 2: isolated, 1D or 2D system
pbc_axis=3                 # set the axis along the 1D direction or 
                           # normal to the 2D plane (pbc_axis = 1, 2 
                           # or 3 for x, y or z axis)
#####################################################################

# clean TMP_DIR
$ECHO "  cleaning $TMP_DIR...\c"
rm -rf $TMP_DIR/*
$ECHO " done"

$ECHO "  running the scf calculation with a diffuse ionic charge"
$ECHO "  density and $pbc_correction periodic boundary correction"

prefix=PtCO_linearizedpb
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
   ibrav = 8
   celldm(1) = 10.6881
   celldm(2) = 0.866025
   celldm(3) = 3.95422
   nat = 10
   ntyp = 3
   ecutwfc = 35
   ecutrho = 280
   occupations = 'smearing'
   degauss = 0.03
   smearing = 'mv'
   nbnd = 80
   tot_charge = 1
   !
/
 &ELECTRONS
   !
   conv_thr = 5.D-9
   diagonalization = 'davidson'
   mixing_beta = 0.2
   electron_maxstep = 200
   !
 /
K_POINTS (automatic)
 1 1 1 0 0 0
ATOMIC_SPECIES
Pt 1 Pt.pbe-nd-rrkjus.UPF
C 1 C.pbe-rrkjus.UPF
O 1 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS (bohr)
Pt       8.061327071   0.098057998   8.992142901
Pt       2.608989366   0.098058283   8.992140585
Pt       0.000036609   4.720846294   8.968756935
Pt       5.335159557   4.721612729   9.380196435
Pt       0.000041121   7.802951963   4.604626508
Pt       5.335161233   7.697749113   4.753489408
Pt       2.697860636   3.152173889   4.688412329
Pt       7.972463687   3.152174491   4.688415209
C        5.335084148   4.646723426  12.901029877
O        5.335009643   4.619623254  15.079854269
EOF
cat > environ.in << EOF 
 &ENVIRON
   !
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
   solvent_temperature = $solvent_temperature
   system_dim = $system_dim            
   system_axis = $system_axis   
   system_ntyp = $system_ntyp
   !
 /
 &BOUNDARY
   !
   solvent_mode = 'full'
   stern_mode = '$stern_mode'
   stern_distance = $stern_distance
   stern_spread = $stern_spread   
   !
 /
 &ELECTROSTATIC
   !
   pbc_correction = '$pbc_correction'
   pbc_dim = $pbc_dim                 
   pbc_axis = $pbc_axis   
   tol = 5.D-13
   !
 /
EOF
   
$PW_COMMAND < $input > $output 
check_failure $?
$ECHO " done"

$ECHO
$ECHO "$EXAMPLE_DIR : done"
