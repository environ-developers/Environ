# Q note: I currently have thee as python scripts
# would not be hard to transfer to a shell script if I get my act together enough



import os
import shlex
import shutil
import subprocess
import numpy as np


def main():
	tot_charge = ['0.001','0.003','0.005']

	epsilon = 11.7
	prefix = 'SiO'

	for q in tot_charge:
		os.mkdir("q_%s"%(q))
		os.chdir("q_%s"%(q))
		make_pw_in(q)
		make_environ_in(q,epsilon)
		os.chdir("./../")


	return



def make_pw_in(q):
	f = open('scf.in','w')
	q_num = float(q)
	f.write("""&CONTROL
  calculation='scf',
  prefix='Si-O',
  !pseudo_dir='./',
  verbosity='high',
  restart_mode='from_scratch',
  nstep=200
  outdir = './tmp'
/
&SYSTEM
  ibrav=0,
  celldm(1)=7.2565749368d0,
  nat=15,
  ntyp=2,
  ecutwfc=50.0d0,
  ecutrho=200.0d0,
  occupations='smearing',
  smearing='mv',
  degauss=0.03000d0,
/
&ELECTRONS
  diagonalization='david',
  conv_thr=1d-10,
  mixing_mode='local-TF',
  mixing_beta=0.500d0,
  electron_maxstep = 250
/
&IONS
  !ion_dynamics = 'damp'
/
ATOMIC_SPECIES
  Si 28.085500d0 Si.upf
  O  15.999000d0 O.UPF
ATOMIC_POSITIONS {alat}
Si       0.000000000   0.000000000   1.802073905    0   0   0
Si       0.500000000   1.060724541   1.802078557    0   0   0
Si       0.000000000   0.353553391   2.302078557    0   0   0
Si       0.499909714   0.707106781   2.302060048    0   0   0
Si       0.000000000   0.000000000   2.802078557    0   0   0
Si       0.500000000   1.059951789   2.802078557    0   0   0
Si      -0.000852798   0.353553391   3.302078557    0   0   0
Si       0.500000000   0.707106781   3.302090970    0   0   0
Si       0.000000000   0.000000000   3.802078557    0   0   0
Si       0.500000000   1.060660172   3.802078557    0   0   0
Si       0.011264507   0.346111979   4.314003460
Si       0.492983512   0.724485161   4.333979006
Si       0.127614179  -0.162293223   4.668709662
Si       0.672109676   0.916491069   4.928620199
O        0.428409191   1.299131531   4.980203547
K_POINTS {automatic}
  2 2 1 1 1 1
CELL_PARAMETERS {alat}
  1.000000000000d0  0.000000000000d0  0.000000000000d0
  0.000000000000d0  1.414213562374d0  0.000000000000d0
  0.000000000000d0  0.000000000000d0  7.624988536903d0""")
	f.close()





	return

def make_environ_in(q,epsilon):
	f = open("environ.in",'w')
	q_num = float(q)
	f.write("""&Environ
    verbose = 3,
    environ_thr = 0.1,
    system_axis = 3,
    system_dim = 2,
    !env_static_permittivity = 2.0
    env_static_permittivity = 78.3,
    env_electrolyte_ntyp = 2,
    temperature = 300,
    electrolyte_linearized = .false.,
    zion(1) = 1,
    cion(1) = 0.01,
    zion(2) = -1,
    cion(2) = 0.01 ,
    sc_permittivity = %s,
    sc_carrier_density = 1.D18,  !! In units of cm^-3
    sc_electrode_chg = %s  ! In units of electrons (positive numbers correspond to positive charge)
    sc_chg_thr = 1.D-5 ! tells the outer semiconductor loop when to cut off its run
/
&BOUNDARY
    sc_spread = 0.01,
    sc_distance = 15.2
    electrolyte_spread = 0.01,
    electrolyte_distance = 15.2,
    electrolyte_mode = 'system',
/
&ELECTROSTATIC
    problem = 'generalized'
    solver = 'iterative'
    auxiliary = 'full'
    inner_solver = 'direct'
    pbc_correction = 'ms-gcs'
    pbc_dim = 2
    pbc_axis = 3
    tol = 1e-11,
/"""%(str(epsilon),str(q)))
	f.close()



	return









if __name__=='__main__':
    main()
