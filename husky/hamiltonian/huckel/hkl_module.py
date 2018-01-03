import ctypes
from ctypes.util import find_library
from numpy.ctypeslib import ndpointer
import os

##############################
# find and load the library
##############################
path = os.path.dirname(os.path.realpath(__file__))
libname = path + '/hkl.so'
if not os.path.isfile(libname):
	raise FileNotFoundError('%s not found' %libname)

# load the library
lib_hkl = ctypes.cdll.LoadLibrary(libname)

# strin pointer type
ctypes_str = ctypes.c_char_p

############################################################
# set the argument type for number of orbtitals
############################################################
lib_hkl.compute_nb_orb.argtypes = [ctypes_str]
lib_hkl.compute_nb_orb.restype = ctypes.c_int

############################################################
# set the argument type for compute_hamiltonian
############################################################
lib_hkl.compute_huckel_hamiltonian_general.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ctypes.c_int, ctypes_str]
lib_hkl.compute_huckel_hamiltonian_general.restype = None

############################################################
# compute the Hamiltonian with huckel
############################################################
def hklHam(h,s,n,filename):
	return lib_hkl.compute_huckel_hamiltonian_general(h,s,n,filename)

############################################################
# compute the number of orbital in the system
############################################################
def nbOrb(filename):
	return lib_hkl.compute_nb_orb(filename)

############################################################
# write the file for individual fragments
############################################################
def writeInputFileFragment(xyz,param,filename,KHT):
	natom = len(xyz)
	f=open(filename,'w')

	# write the header
	f.write('nb_atom    \t %d\n' %natom)
	f.write('parameters \t %s\n' %param)
	f.write('Keht	     \t %f\n\n' %KHT)

	# write the positions
	for i in range(natom):
		at = xyz[i].symbol
		x,y,z = list(map(float,xyz.positions[i,:]))
		f.write('%s\t% f\t% f\t% f\n' %(at,x,y,z) )
	f.write('////////////////////////////////////////////////////\n\n')
	f.close()

############################################################
# write the file for total system
############################################################
def writeInputFileJunction(xyz_central,dict_elec,filename,KHT,hkl_param_file):

	# number of electrodes
	nelec = len(dict_elec)

	# compute the number total of atom
	natom = len(xyz_central)
	for ielec in range(nelec):
		natom += len(dict_elec[str(ielec)].xyz)
	
	# write the header
	f=open(filename,'w')
	f.write('nb_atom    \t %d\n' %natom)
	f.write('parameters \t %s\n' %hkl_param_file)
	f.write('Keht	     \t %f\n\n' %KHT)

	# write the electrode positions
	for ielec in range(nelec):
		xyz = dict_elec[str(ielec)].xyz
		nat = len(xyz) 
		for iat in range(nat):
			at = xyz[iat].symbol
			x,y,z = list(map(float,xyz.positions[iat,:]))
			f.write('%s\t% f\t% f\t% f\n' %(at,x,y,z) )

	#write the molecule positions
	nat = len(xyz_central)
	for iat in range(nat):
		at = xyz_central[iat].symbol
		x,y,z = list(map(float,xyz_central.positions[iat,:]))
		f.write('%s\t% f\t% f\t% f\n' %(at,x,y,z) )
	f.write('////////////////////////////////////////////////////\n\n')
	f.close()
