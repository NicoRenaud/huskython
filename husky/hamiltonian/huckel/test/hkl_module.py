import ctypes
from ctypes.util import find_library
from numpy.ctypeslib import ndpointer
import os


# find and load the library
path = os.getcwd()
lib_hkl = ctypes.cdll.LoadLibrary(find_library('/Users/nicolasrenaud/Documents/CODE/myCodes/huskython/huskython/hamiltonian/huckel/hkl'))

# set the argument type for number of orbtitals
lib_hkl.compute_nb_orb.argtypes = [ctypes.c_char_p]
lib_hkl.compute_nb_orb.restype = ctypes.c_int

# set the argument type for compute_hamiltonian
lib_hkl.compute_huckel_hamiltonian_general.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_char_p]
lib_hkl.compute_huckel_hamiltonian_general.restype = None

def hklHam(h,s,n,filename):
	return lib_hkl.compute_huckel_hamiltonian_general(h,s,n,filename)

def nbOrb(filename):
	return lib_hkl.compute_nb_orb(filename)