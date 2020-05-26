import matplotlib.pyplot as plt
import numpy as np
from struct import pack
import glob
import re
import os
import math
import pickle
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Avenir Next Condensed'
mpl.rcParams['font.size'] = 16


#   A default graphing script to graph the planar average of a cube file
#   Assumes you're interested in the z direction of your slab, can be changed



def main():


#-- Define conversion factors
    Bohr2Ang = 0.529177249
    Ry2eV = 13.605698066


    fermi_0V = find_fermi()

    #Process Cube file input
    cube_0V='./velectrostatic.cube'		#Directory of potential you're interested in
    grid_0V, potential_0V, atom_0V, cell_0V = cube_reader(cube_0V)
    z = np.array(grid_0V[2, 1, 1, :]) * Bohr2Ang
    v_0V = np.array([ np.mean(potential_0V[:,:,i]) for i in range(potential_0V.shape[2]) ]) * Ry2eV
    macro_v_0V = macro_avg(z,v_0V,3.56)
    #pickle.dump( macro_v_0V, open( "macro_v_0V.p", "wb" ) )     	##Dumps planar average into a pickle file (*.p) that can be read in later without using read cube again
    #pickle.dump(z,open("z.p","wb"))
    #print shift_0V, fermi_0V

    plotname='q_0_potential.png'
    #-- Plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlim(0, z.max(axis=0))
    #plt.ylim(-22, 2)
    #plt.title(r'Potential')
    plt.xlabel(r'$z$, ($\mathrm{\AA}$)')
    plt.ylabel('electrostatic potential (V)')
    ax.plot(z,v_0V)						#Plots the potential
    plt.axhline(y=fermi_0V,linestyle='--',color='k')		#Plots a horizontal line for the fermi level
    #plt.legend(['0','0.2','0.4','0.6','0.8','1.0'],loc=2)


    #Uncomment and modify the next section to create "broken axis" graph which Dr. Dabo likes

    #ax.set_xticks([0,5,10,15])
    #ax.set_yticks([0,10,20])
    #ax.spines['left'].set_bounds(0,20)
    #ax.spines['bottom'].set_bounds(0,15)
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    #ax.yaxis.set_ticks_position('left')
    #ax.xaxis.set_ticks_position('bottom')
    #ax.xaxis.set_tick_params(width=1)
    #ax.yaxis.set_tick_params(width=1)


    plt.tight_layout()
    plt.savefig(plotname)
    plt.show()
    plt.close(fig)
    #print "Saved fig"


    return


### Finds the Fermi level within the scf.out file

def find_fermi():
    output = open('pw.out')
    for line in output:
        line = line.rstrip()
        if re.search('Fermi', line) :
            fermi_level = float(line.split()[4])
	    #print fermi_level
	    break

    try:
    	return  fermi_level
    except:
	return 0.0



### New Macroscopic averager

def macro_avg(z,v_avg,d):

    #z is the array of position values (in Angstrom)
    #v_avg is the planar average array
    #d is the macroscopic averaging distance (in Angstrom)
    #Ideally d should be the interplanar spacing of your slab



    ##Find the number of arrray units corresponding to the distance d
    for i in range(len(z)):
	if z[i] >d:
	    array_d = i
	    break

    ##Make sure array size is even
    if array_d%2 != 0:
	array_d = array_d + 1

    macro_avg = np.ones(len(z))

    #Creating macroscopic average
    for i in range(len(v_avg)):
	lower_bound = i - array_d/2
	upper_bound = i + array_d/2

	#Handling Averaging on the edge of the cell
	#This portion is most likely to change depending on how average is desired.
	if lower_bound < 0:
		lower_bound = 0
	if upper_bound > len(v_avg):
		upper_bound = len(v_avg)

	macro_avg[i] = np.mean(v_avg[lower_bound:upper_bound])

    return macro_avg

#########################
#  Do not modify below  #
# cube_reader stolen from Steve Weitzner
#########################

def cube_reader(cube_file):
   """
    Extract numerical data from Gaussian *.cube files.

    Parameters
    ----------
    cube_file : string
                Name of cube file to read

    Returns
    -------
    grid : np.ndarray, shape [3, Nx, Ny, Nz]
           4D tensor which contains the cartesian coordinates
           of the numerical grid in units of Bohr.

    Data3D : np.ndarray, shape [Nx, Ny, Nz]
             3D tensor which contains scalar field data on
             the numerical grid in atomic units.

             If charge data, in units of '[e]'.
             If potential data, in units of '[Ry/e]'.

    atom_data : np.ndarray, shape [N_at, 5]
                Each atom (row) contains the following (col):
                    [0] Atomic Number
                    [1] Atomic Mass (a.m.u.)
                    [2] x-component of position
                    [3] y-component of position
                    [4] z-component of position

                Example for an isolated water molecule (order O, H, H):
                    [[  8    15.9994   11.791224   12.043800   11.500000 ],
                     [  1    1.00794   13.443082   11.221856   11.500000 ],
                     [  1    1.00794   10.565694   10.664344   11.500000 ]]

    cell_data : np.ndarray, shape [3,3]
                Each column contains a basis vector of the supercell
                in units of Bohr.

    References
    ----------
    [1] http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm

    To Do
    -----
    -> Package as a module to enable system wide import statements.
    """

   #-- Read cube file into list
   with open(cube_file, 'r') as f:
       Data = f.readlines()

   #-- Parse the header of the cube file
   del Data[0:2] #-- remove first 2 lines
   Num_Atoms = int(Data[0].split()[0])
   Header = Data[0 : Num_Atoms + 4]
   Origin = map(float, Header[0].split()[1:4])
   N1 = int(Header[1].split()[0])
   N2 = int(Header[2].split()[0])
   N3 = int(Header[3].split()[0])
   R1 = map(float, Header[1].split()[1:4])
   R2 = map(float, Header[2].split()[1:4])
   R3 = map(float, Header[3].split()[1:4])

   #-- Get atom types and coordinates
   atom_data = [ ]
   for line in Header[4:]:
       atom_data += [ line.split() ]
   atom_data = np.array( atom_data, dtype='d' )
   #-- Get supercell dimensions
   cell_data = np.array( [ R1, R2, R3 ], dtype='d' ).T
   scalars = np.array( [N1, N2, N3], dtype='d' ).T
   cell_data = cell_data * scalars  #  broadcasting
   #-- Construct grid
   grid = np.empty([ N3, N2, N1, 3 ], dtype='d')
   for i in range(N1):
       for j in range(N2):
           for k in range(N3):
               X = Origin[0] + i*R1[0] + j*R2[0] + k*R3[0]
               Y = Origin[1] + i*R1[1] + j*R2[1] + k*R3[1]
               Z = Origin[2] + i*R1[2] + j*R2[2] + k*R3[2]
               grid[k, j, i] = [ X, Y, Z ]
   #-- Isolate potential data
   del Data[0 : Num_Atoms + 4]
   Data1D = [ ]
   for line in Data:
       Data1D += [ float(item) for item in line.split() ]
   #-- Re-map data into a 3D array
   Data3D = np.empty([N3, N2, N1], dtype='d')
   for i in range(N1):
       for j in range(N2):
           for k in range(N3):
               idx = k + N3*j + N3*N2*i
               Data3D[k,j,i] = Data1D[idx]

   return np.ascontiguousarray(grid.T), np.ascontiguousarray(Data3D.T), atom_data, cell_data


if __name__ == "__main__":
    main()
