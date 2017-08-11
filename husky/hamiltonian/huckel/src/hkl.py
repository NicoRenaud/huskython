import numpy as np
import ase.io 

#######################################
# class for a single orbital
#######################################
class single_orb(object):

	def __init__(self,n=None,VSIP=None,exp=None,exp2=None,coef2=None):
		self.orb_of_e = n
		self.VSIP = VSIP
		self.exp = exp 
		self.exp2 = exp2
		self.coef2 = coef2

#######################################
# class for a single atoms
#######################################
class atom_parameter(object):

	def __init__(self):
		self.symbol = ''
		self.valence_electrons = 0
		self.norb_tot = 0
		self.orb = []
		for i in range(4):
			self.orb.append(single_orb())


class huckel(object):

	#############################################
	# init the class
	#############################################
	def __init__(self,xyz=None,param=None,param_format='new'):


		self.xyz = ase.io.read(xyz)
		self.param_format = param_format
		self.period = self.read_param(param)
		
		self.natom = len(self.xyz)		
		self.ntype = len(self.param)

		self.norb_tot = None
		self.atom_types = self.find_index_type()

		self.Keht = 1.75
		self.H = []
		self.S = []

	#############################################
	# read the parameter file
	#############################################
	def read_param(self,param):

		period = []
		nparam_per_orb = 2
		norb = [1,3,5,7]

		# read the data
		f = open(param,'r')
		alldata = f.readlines()
		
		# extract each line
		for i in range(len(alldata)):

			# split it
			data = alldata[i].split()

			# create the atom and read the first parameters
			atom = atom_parameter()
			atom.symbol = data[0]
			atom.valence_electrons = data[1]

			# extract the principal quantum number
			# and the slater coeffs
			pqn = data[2:6]
			slat_coef = data[6:]

			# read the principal quantum number of each orbital
			# and the parameter of each orbitals
			for i in range(4):
				atom.orb[i].orb_of_e = pqn[i]
				if pqb[i]>0:
					atom.orb[i].VSIP = slat_coef[i*nparam_per_orb]
					atom.orb[i].exp = slat_coef[i*nparam_per_orb+1]

			# total number of orbitals
			atom.norb_tot = np.sum((np.array(pqn)!=0)*np.array(norb))
			period.append(atom)

		return period

	#############################################
	# find the index of atom type for each atom
	# in the xyz data
	#############################################
	def find_index_type(self):
		atom_types = []
		for iatom in range(self.natom):
			for itype in range(self.ntype):
				if xyz[iatom].symbol.lower() == self.param[itype][0].lower():
					atom_types.append(itype)
					break
		return atom_types


	#############################################
	# Compute the overlap matrix of the system
	#############################################	
	def compute_matrices(self):

		# determine the size of the matrices
		self.norb_tot = 0
		for iat in range(self.natom):
			self.norb_tot += self.period[find_index_type[iat]].norb_tot

		# declare the matrices
		self.H = np.zeros((self.nb_orb,self.nb_orb))
		self.S = np.zeros((self.nb_orb,self.nb_orb))

		# loop over the atoms
		for iat in range(self.natom):

			for jat in range(self.natom):









