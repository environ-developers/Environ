import os
import shlex
import shutil
import subprocess
import numpy as np
import re


def main():
	prefix = 'Si-OH'	
	tot_charge = ['0.001','0.003','0.005','0.007','0.009','0.011','0.018']
	f = open('%s.q-v.dat'%prefix,'w')
	f.write(""" Potential (V-V_fb)  SS potential (V_cut-V_fb)  Electrode Charge (e)  Surface States Charge (e)    
 Electrode Charge per surface area (e/cm^2)     
 Surface State Charge per surface area (e/cm^2)\n""")

	for q in tot_charge:
		os.chdir("q_%s"%(q))
		
		if check_finished():
			g = open('q-v.dat','r')
			lines = g.readlines()
			f.write(lines[3])
			g.close()			
		else:
			print "q = %s did not converge"%(q)



		os.chdir("./../")

	f.close()

	return



def check_finished():

    output = open('scf.out')
    finished = False
    for line in output:
        line = line.rstrip()
        if re.search('q-v.dat', line) :
            finished = True
	    break



    return finished






if __name__=='__main__':
    main()



