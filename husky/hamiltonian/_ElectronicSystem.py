
import numpy as np
import ase.io 
from huckel import hkl_module as hkl
from huckel.transform import lowdin
from huckel.orbitals import print_mo_mopac


###################################################
# Main class that defines the junction
###################################################

class ElectronicSystem(object):

	def __init__(self,mol = None, elec = None):

		# molecules
		self.xyz_mol = []
		self.norb_mol = []
		self.index_orb_mol = []
		self.index_orb_central = []
		self.nmol = 0
		
		# electrodes
		self.xyz_elec = []
		self.norb_elec = []
		self.index_orb_elec = []
		self.nelec = 0

		# hamiltonian of the molecules
		self.Hmol = {}
		self.Smol = {}

		# interaction between different molecules
		self.VinterMol = {}
		self.SinterMol = {}

		# interactions elctrodes molecules dictionaries
		self.VelecMol = {}
		self.SelecMol = {}

		# hamiltonian of the electrode
		self.Helec = {}
		self.Selec = {}

		# interaction between unit cell of the electrodes
		self.VintraElec = {}
		self.SintraElec = {}

		# size and parameters
		self.norb_total = 0
		self.hkl_param_file = None

	##########################################
	# ADD XYZ TO DEFINE THE SYSTEM
	##########################################

	# add a molecule to the junction
	def add_molecule(self,xyz_name):
		self.xyz_mol.append(ase.io.read(xyz_name))
		self.nmol += 1

	# add an electrode to the junction
	def add_electrode(self,xyz_name):
		self.xyz_elec.append(ase.io.read(xyz_name))
		self.nelec += 1

	# add an huckel parameter file
	def add_huckel_parameters(self,file_name):
		self.hkl_param_file = file_name

	##########################################
	# DIRECTLY ADD HAMILTONIAN MATRIX
	########################################## 

	# ELECTRODE HAMILOTNIAN
	def add_electrode_hamiltonian(self,h,s=None,ielec=None):

		if len(h.shape) == 1:
			h = np.matrix(h)

		if s == None:
			s = np.eye(len(h))

		if ielec ==None :
			ielec = self.nelec

		self.Helec[str(ielec)] = h
		self.Selec[str(ielec)] = s

		self.nelec += 1

	#  MOLECULE HAMILTONIAN
	def add_molecular_hamiltonian(self,h,s=None):
		if s == None:
			s = np.eye(len(h))
		self.Hmol['allMol'] = h
		self.Smol['allMol'] = s

	# INTERACTION ELECTRODE MOLECULE
	def add_electrode_molecule_interaction(self,v,s=None,ielec=None):

		if len(v.shape)==1:
			v = np.array(np.matrix(v))

		if s == None:
			s = np.zeros(v.shape)

		if ielec ==None :
			ielec = self.nelec

		self.VelecMol['%d->allMol' %ielec] = v
		self.SelecMol['%d->allMol' %ielec] = s

		self.nelec += 1
		self.norb_elec.append(len(v))


	# INTERACTION ELECTRODE MOLECULE
	def add_intra_electrode_interaction(self,v,s=None,ielec=None):

		if len(v.shape)==1:
			v = np.array(np.matrix(v))

		if s == None:
			s = np.zeros(v.shape)

		if ielec ==None :
			ielec = self.nelec

		self.VintraElec['%d' %self.nelec] = v
		self.SintraElec['%d' %self.nelec] = s

	
	########################################
	# compute the matrices with husky
	########################################
	def compute_matrices_huckel(self,KHT=1.75,lowdin_ortho=False,print_mo=True):

		# compute the orbitals of the molecules and orbitals
		hkl_temp = '__hkl_temp.in'
		nborb_tot = 0
		
		for i in range(self.nelec):
			hkl.writeInputFileFragment(self.xyz_elec[i],self.hkl_param_file,hkl_temp,KHT)
			norb_frag = hkl.nbOrb(hkl_temp)
			self.norb_elec.append(norb_frag)
			self.index_orb_elec.append(range(nborb_tot,nborb_tot+norb_frag))
			nborb_tot += norb_frag

		for i in range(self.nmol):
			hkl.writeInputFileFragment(self.xyz_mol[i],self.hkl_param_file,hkl_temp,KHT)
			norb_frag = hkl.nbOrb(hkl_temp)
			self.norb_mol.append(norb_frag)
			self.index_orb_mol.append(range(nborb_tot,nborb_tot+norb_frag))
			nborb_tot += norb_frag


		# write the input file for the junction
		hkl.writeInputFileJunction(self,hkl_temp,KHT)
		self.norb_total = hkl.nbOrb(hkl_temp)

		# compute the matrices
		h = np.zeros(self.norb_total*self.norb_total)
		s = np.zeros(self.norb_total*self.norb_total)
		hkl.hklHam(h,s,self.norb_total,hkl_temp)
		h = h.reshape(self.norb_total,self.norb_total)
		s = s.reshape(self.norb_total,self.norb_total)

		# lowdin if necessary
		if lowdin_ortho:
			h = lowdin(h,s)
			s = np.eye(self.norb_total)

		# extract the molecular submatrices
		for imol in range(self.nmol):
			self.Hmol['%d' %imol] = h[np.ix_(self.index_orb_mol[imol],self.index_orb_mol[imol])]
			self.Smol['%d' %imol] = s[np.ix_(self.index_orb_mol[imol],self.index_orb_mol[imol])]

		# extract the total molecular matrices
		self.index_orb_central = [index for molindex in self.index_orb_mol for index in molindex]
		self.Hmol['allMol'] = h[np.ix_(self.index_orb_central,self.index_orb_central)]
		self.Smol['allMol'] = s[np.ix_(self.index_orb_central,self.index_orb_central)]

		# extract the coupling between molecules
		for imol in range(self.nmol-1):
			v = h[np.ix_(self.index_orb_mol[imol],self.index_orb_mol[imol+1])]
			sv = s[np.ix_(self.index_orb_mol[imol],self.index_orb_mol[imol+1])]
			self.VinterMol['%d->%d' %(imol,imol+1)] = v
			self.SinterMol['%d->%d' %(imol,imol+1)] = sv

		# extract the electrode-molecule interactions
		for ielec in range(self.nelec):
			
			for imol in range(self.nmol):

				v = h[np.ix_(self.index_orb_elec[ielec],self.index_orb_mol[imol])]
				sv = s[np.ix_(self.index_orb_elec[ielec],self.index_orb_mol[imol])]

				self.VelecMol['%d->%d' %(ielec,imol)] = v
				self.SelecMol['%d->%d' %(ielec,imol)] = sv

			self.VelecMol['%d->allMol' %ielec] = h[np.ix_(self.index_orb_elec[ielec],self.index_orb_central)]
			self.SelecMol['%d->allMol' %ielec] = s[np.ix_(self.index_orb_elec[ielec],self.index_orb_central)]

		# print the molecular orbitals
		if print_mo:

			if self.nmol>1:
				for imol in range(self.nmol):
					w,u = np.linalg.eigh(self.Hmol['%d' %imol])
					fname = 'mo_mol%d.dat' %imol
					print_mo_mopac(w,u,self.Smol['%d' %imol],self.xyz_mol[imol],fname,self.hkl_param_file)

			w,u = np.linalg.eigh(self.Hmol['allMol'])
			print_mo_mopac(w,u,self.Smol['allMol'],self.xyz_mol[imol],'mo_central.dat',self.hkl_param_file)





