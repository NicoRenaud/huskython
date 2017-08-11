import os
import numpy as np
import ase.io 
from huckel import hkl_module as hkl
from huckel.transform import lowdin
from huckel.orbitals import print_mo_mopac
import scipy.linalg as spla
import time 

try:
	from tqdm import tqdm
except:
	def tqdm(a):
		return a

##################################################
# electrode class
##################################################
class Electrode(object):

	def __init__(self,xyz=None,h0=None,s0=None,hi=None,si=None,vmol=None,smol=None):


		# XYZ : the atomic coordinate of the elctrode
		# the xyz coordinates myst contains two unit cells. 
		# the first unit cell must be the one closest to the molecule

		# VMOL/SMOL : interaton/overlap netween the electrode and the molecule
		# when input by hand, the matrix must be norb_elec X norb_mol
		
		self.xyz = xyz

		# H/S og the unit cell
		self.h0 = h0
		self.s0 = s0

		# H/S of the interaction between cells
		self.hi = hi
		self.si = si

		# H/S interaction with molecules
		self.vmol = vmol
		self.smol = smol

		# dimension
		if h0 is not None:
			self.norb = len(h0)
		else:
			self.norb = 0

		self.index_orb = [[],[]]

		# band structure
		self.K = None
		self.bands = []

		# surface green function
		self.sgf = []
		self.sgf_energies = []

		#  channels for esqc
		self.channels = None
		self.index_in = None
		self.index_out = None

		# used in the ESQCsolver 
		# index of the non propagating channels 
		# they are obsolete as this esqc solver is not used
		self.nprop_channel = []
		self.prop_channel = None
		self.Umat = None
		self.Xmat = None



	#################################################################################
	#
	#		GREEN FUNCTION REALTED ROUTNIES
	#
	#################################################################################	

	# compute the green function for all energies and plot the result
	def surface_green_function(self,energies,eps=1E-2,tol=1E-8,itermax=5000,filename='surface_gf.png'):

		self.sgf = []
		t0 = time.time()
		print ' - Compute surface green fuction of the electrode'
		for iE,e in enumerate(tqdm(energies,ncols=50,desc='   ',leave=True)):
			self.sgf.append(self.compute_green_function(e,eps=eps,tol=tol,itermax=itermax))
		print ' - Surface Green function computed in %1.3f sec. ' %(time.time()-t0)
		self.plot_surface_green_function(energies)

	# compute the green function for one given energy
	def compute_green_function(self,E,eps=1E-2,tol=1E-8,itermax=5000):

		if self.hi is None:
			print "Interaction between layers not defined. Can't compute the green function of the electrode"
			return

		# init the GF
		green = np.linalg.inv((E+1j*eps)*self.s0 - self.h0)

		# reccursively compute the GF
		niter, err = 1, 1.
		while niter<itermax and err>tol :

			gf_old = green
			green = np.linalg.inv( (E+1j*eps)*self.s0 - self.h0 - np.dot((self.hi-E*self.si),np.dot(gf_old,(self.hi-E*self.si))))
			err = np.sum(np.abs(np.abs(green)-np.abs(gf_old))).real
			niter += 1

		if niter==itermax:
			print 'Surface green function has not converged after %d iterations.\nEnergy % 1.3f err/tol = % 1.3E/% 1.3E' %(itermax,E,err,tol)
		
		
		return green


	# interpolate the surface green energy
	def interpolate_surface_green_function(self,E):

		# see if we have a perfect match
		index_match = np.where(self.sgf_energies==E)[0]
		if len(index_match)>0:
			return self.sgf[index_match][0]

		# otherwise we linearly interpolate the sgf
		ind = np.searchsorted(self.sgf_energies,E)
		if ind == 0:
			print 'warning: Interpolation of the SGF outside precomputed range'
			return self.sgf[0]
		elif ind == len(self.sgf_energies):
			print 'warning: Interpolation of the SGF outside precomputed range'
			return self.sgf[-1]
		else:
			alpha = (E-self.sgf_energies[ind])/(self.sgf_energies[ind+1]-self.sgf_energies[ind])
			return self.sgf[ind] + alpha*(self.sgf[ind+1]-self.sgf[ind])

	# plot the surface green function
	def plot_surface_green_function(self,energies,filename='sgf.png'):
		
		import matplotlib.pyplot as plt
		gf_plot = np.array([ np.trace(self.sgf[i]) for i in range(len(self.sgf))])
		plt.plot(energies,gf_plot.real,color='blue')
		plt.plot(energies,gf_plot.imag,color='red')
		plt.xlabel('Energies')
		plt.ylabel('Surface GF')
		plt.savefig(filename)
		plt.close()

	#################################################################################
	#
	#		BAND STRUCTURE AND DOS RELATED ROUTINES
	#
	#################################################################################

	# compute the 1D band structure of the electrodes
	def compute_bands(self,nK = 100):
		
		if self.hi is None:
			print "Interaction between layers not defined. Can't cpmpute the bands"
			return

		self.K = np.linspace(0,np.pi,nK)
		for ik,k in enumerate(self.K):
			H = self.h0 + np.exp(1j*k) * self.hi + np.exp(-1j*k) * self.hi.T
			S = self.s0 + np.exp(1j*k) * self.si + np.exp(-1j*k) * self.si.T
			l,u = spla.eigh(H,b=S)
			self.bands.append(l)

		
		self.bands=np.array(self.bands).T
		self.nbands = len(self.bands)
		self.nK = nK

	# compute the densitu of states 
	def compute_dos(self,nE=100):

		if self.bands == []:
			print "The band structure is missing. Can't compute the DOS"
			return			

		emin,emax = np.amin(self.bands),np.amax(self.bands)
		self.Edos = np.linspace(emin,emax,nE)
		self.dos = np.zeros(nE)
		dK = self.K[1]-self.K[0]
		eps = 1e-3

		if self.nbands==1:
			invdE = np.asmatrix(1./(np.gradient(self.bands[0],dK)))
		else:
			invdE = 1./((np.gradient(self.bands,dK,edge_order=2))[0])

		for iE in range(self.nbands):
			for iK in range(self.nK):
				index = np.argmin(np.abs(self.Edos-self.bands[iE,iK]))
				self.dos[index] += (invdE[iE,iK]*np.conj(invdE[iE,iK])).real

	# plot the band structure alone
	def plot_bands(self,filename='bands.png'):
		import matplotlib.pyplot as plt 
		plt.plot(self.K,self.bands.T,color='b')
		plt.xlabel('K-points')
		plt.ylabel('Energies')
		plt.savefig(filename)
		plt.close()

	# plot the dos alone
	def plot_dos(self,filename='dos.png'):

		import matplotlib.pyplot as plt 
		
		plt.plot(self.Edos,self.dos,color='b')
		plt.ylabel('dos')
		plt.xlabel('Energies')
		plt.savefig(filename)

	# comute/represent the bands and the dos
	def band_structure(self,nK=100,nE=100,filename='bandstructure.png'):

		import matplotlib.pyplot as plt 
		from matplotlib import gridspec

		self.compute_bands(nK = nK)
		self.compute_dos(nE = nE)

		fig = plt.figure(figsize=(8,6))
		gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

		ax0 = plt.subplot(gs[0])
		ax0.plot(self.K, self.bands.T,color='black')
		plt.xlabel('k_points')
		plt.ylabel('Energies')

		ax1 = plt.subplot(gs[1])
		ax1.plot(self.dos,self.Edos)
		plt.xlabel('dos')
		plt.setp(ax1.get_xticklabels(), visible=False)

		plt.tight_layout()
		plt.savefig(filename)
		plt.close()

##################################################
# Molecule class
##################################################
class Molecule(object):

	def __init__(self,xyz=None,h0=None,s0=None):

		self.xyz = xyz
		self.h0 = h0
		self.s0 = s0

		if h0 is not None:
			self.norb = len(h0)
		else:
			self.norb = 0
		self.index_orb = []


##################################################
# Junction class
##################################################
class Junction(object):

	def __init__(self,electrodes = {}, molecules = {}):

		# the electrodes
		self.nelec = 0
		self.electrodes = electrodes		

		# the molecular fragments
		self.nmol = 0
		self.molecules = molecules
		
		# the total central region
		self.central_region = Molecule()

		# hkl parameter file
		self.hkl_param_file = None
		self.norb_tot = 0

		# interaction between molecules
		self.molecule_interaction_hamiltonian = {}
		self.molecule_interaction_overlap = {}

	##########################################
	# ADD XYZ TO DEFINE THE SYSTEM
	##########################################

	# add a molecule to the junction
	def add_molecule_xyz(self,xyz_name,imol=None):

		if imol is None:
			imol = self.nmol

		self.molecules[str(imol)] = Molecule(ase.io.read(xyz_name))
		self.nmol += 1

	# add an electrode to the junction
	def add_electrode_xyz(self,xyz_name,ielec=None,trans=[]):

		if ielec is None:
			ielec = self.nelec

		self.electrodes[str(ielec)] = Electrode()
		self.electrodes[str(ielec)].xyz = ase.io.read(xyz_name)

		# if we provide a translation vector it means
		# that we want to automatically construct
		# the second unitcell of the electrode
		if trans != []:
			tmp = self.electrodes[str(ielec)].xyz.copy()
			tmp.translate(trans)
			self.electrodes[str(ielec)].xyz.extend(tmp)

		self.nelec += 1

	# write a xyz file of the entire junction
	def write_xyz(self,filename='junction.xyz'):

		# number of electrodes
		nelec = len(self.electrodes)
		nmol = len(self.molecules)

		# compute the number total of atom
		natom = 0
		for imol in range(nmol):
			natom += len(self.molecules[str(imol)].xyz) 
		for ielec in range(nelec):
			natom += len(self.electrodes[str(ielec)].xyz)

		f = open(filename,'w')
		f.write('%d\n\n' %(natom))

		# write the electrode positions
		for ielec in range(nelec):
			xyz = self.electrodes[str(ielec)].xyz
			nat = len(xyz) 
			for iat in range(nat):
				at = xyz[iat].symbol
				x,y,z = map(float,xyz.positions[iat,:])
				f.write('%s\t% f\t% f\t% f\n' %(at,x,y,z) )

		#write the molecule positions
		for imol in range(nmol):
			xyz = self.molecules[str(imol)].xyz
			nat = len(xyz)
			for iat in range(nat):
				at = xyz[iat].symbol
				x,y,z = map(float,xyz.positions[iat,:])
				f.write('%s\t% f\t% f\t% f\n' %(at,x,y,z) )
		f.close()


	##########################################
	# ADD HKL PARAM
	##########################################

	# add an huckel parameter file
	def add_huckel_parameters(self,file_name):
		self.hkl_param_file = file_name


	############################################
	# ADD DIRECTLY HAMILTONAIN MATRICES
	############################################

	# electrode matrices
	def add_electrode_matrix(self,h0=None,s0=None,hi=None,si=None,vmol=None,smol=None,ielec=None):

		# deal with the index number
		if ielec is None:
			ielec = np.array([self.nelec])
				
		# deal with the overal of the cell
		if h0 is not None:
			if h0.ndim == 1:
				h0 = np.array(np.matrix(h0))
			if s0 is None:
				s0 = np.eye(len(h0))
			else :
				s0 = np.array(np.matrix(s0))			
					
		# deal with the overal of the interaction in the elec
		if hi is not None :
			if hi.ndim==1:
				hi = np.array(np.matrix(hi)) 
			if si is None:
				si  = np.zeros(hi.shape)
			else:
				si = np.array(np.matrix(si))

		# deal with the overal of electro/mol interaction
		if vmol is not None :
			if vmol.ndim == 1:
				vmol = np.array(np.matrix(vmol))
			if smol is None:
				smol  = np.zeros(vmol.shape)
			else:
				smol = np.array(np.matrix(smol))

		# store the matrices
		for i in ielec:
			key = str(i)
			self.electrodes[key] = Electrode(h0=h0,s0=s0,hi=hi,si=si,vmol=vmol,smol=smol)
			self.nelec += 1

	# molecular matrices
	def add_molecule_matrix(self,h0=None,s0=None):

		if s0 is None and h0 is not None:
			s0 = np.eye(len(h0))

		self.molecules[str(self.nmol)] = Molecule(h0=h0,s0=s0)
		self.nmol = 1


	# cenral region matrices
	def add_central_matrix(self,h0=None,s0=None):

		if s0 is None and h0 is not None:
			s0 = np.eye(len(h0))

		self.central_region = Molecule(h0=h0,s0=s0)


	# central region
	def create_central_region(self):

		# atomic positions and orbital indexes
		self.central_region.index_orb = self.molecules['0'].index_orb
		self.central_region.xyz = self.molecules['0'].xyz
		
		for imol in range(1,self.nmol):
			self.central_region.xyz += self.molecules[str(imol)].xyz
			self.central_region.index_orb += self.electrodes[str(imol)].index_orb



	########################################
	#
	#	compute the surface green function 
	#   of the electrodes
	#
	########################################
	def precompute_electrodes_surface_green_function(self,energies,eps=1E-2,tol=1E-8,itermax=5000,filename='surface_gf.png',identical_electrodes=True):
		first_electrode = True
		for index, elec in self.electrodes.iteritems():
			if first_electrode:
				elec.surface_green_function(energies,eps=eps,tol=tol,itermax=itermax,filename=filename)
				elec.sgf_energies = energies
				first_electrode = False
				key_first_electrode = index
			else:
				if identical_electrodes:
					elec.sgf = self.electrodes[key_first_electrode].sgf
				else:
					elec.surface_green_function(energies,eps=eps,tol=tol,itermax=itermax,filename=filename)
				elec.sgf_energies = energies

	# we can also save a sgf and load them in new calculations as it is an expensiv
	# step that only depends on the electrode stucture. usefull for example in MD calculations
	# would be better to use pickle instead of two differnet files
	def load_electrodes_surface_green_function(self,sgf_filename, sgf_energies_filename,ielec=None):
		sgf = np.loadtxt(sgf_filename)
		energies = np.loadtxt(sgf_energies_filename)
		if ielec==None:
			for index, elec in self.electrodes.iteritems():
				elec.sgf_energies = energies
				elec.sgf = sgf
		else:
			self.electrodes[str(ielec)].sgf = sgf 
			self.electrodes[str(ielec)].sgf_energies = energies

	# save the surface green function to load it later 
	# in order to save time
	# would be better to use pickle instead of two differnet files
	def save_electrodes_surface_green_function(self,sgf_filename, sgf_energies_filename):
		np.savetxt(sgf_filename,self.electrodes[str(ielec)].sgf)
		np.savetxt(sgf_energies_filename,self.electrodes[str(ielec)].sgf_energies)		

	########################################
	#
	#    compute the matrices with husky
	#
	########################################
	def compute_matrices_huckel(self,KHT=1.75,lowdin_ortho=False,print_mo=True,clean=True):



		# compute the orbitals of the molecules and orbitals
		hkl_temp = '__hkl_temp.in'
		nborb_tot = 0
		
		for i in range(self.nelec):
			hkl.writeInputFileFragment(self.electrodes[str(i)].xyz,self.hkl_param_file,hkl_temp,KHT)
			norb_frag = hkl.nbOrb(hkl_temp)/2
			self.electrodes[str(i)].norb = norb_frag
			self.electrodes[str(i)].index_orb[0] = (range(nborb_tot,nborb_tot+norb_frag))
			self.electrodes[str(i)].index_orb[1] = (range(nborb_tot+norb_frag,nborb_tot+2*norb_frag))
			nborb_tot += 2*norb_frag

		for i in range(self.nmol):
			hkl.writeInputFileFragment(self.molecules[str(i)].xyz,self.hkl_param_file,hkl_temp,KHT)
			norb_frag = hkl.nbOrb(hkl_temp)
			self.molecules[str(i)].norb = norb_frag
			self.molecules[str(i)].index_orb = range(nborb_tot,nborb_tot+norb_frag)
			nborb_tot += norb_frag

		# extract the total molecular matrices
		self.create_central_region()

		# write the input file for the junction
		hkl.writeInputFileJunction(self.central_region.xyz,self.electrodes,hkl_temp,KHT,self.hkl_param_file)
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
			self.molecules[str(imol)].h0 = h[np.ix_(self.molecules[str(imol)].index_orb,self.molecules[str(imol)].index_orb)]
			self.molecules[str(imol)].s0 = s[np.ix_(self.molecules[str(imol)].index_orb,self.molecules[str(imol)].index_orb)]

		# extract the central Hamilotnian/overlap
		self.central_region.h0 = h[np.ix_(self.central_region.index_orb,self.central_region.index_orb)]
		self.central_region.s0 = s[np.ix_(self.central_region.index_orb,self.central_region.index_orb)]
		self.central_region.norb = len(self.central_region.index_orb)		

		# extract the coupling between molecules
		for imol in range(self.nmol-1):
			v = h[np.ix_(self.index_orb_mol[imol],self.index_orb_mol[imol+1])]
			sv = s[np.ix_(self.index_orb_mol[imol],self.index_orb_mol[imol+1])]
			self.molecule_interaction_hamiltonian['%d->%d' %(imol,imol+1)] = v
			self.molecule_interaction_overlap['%d->%d' %(imol,imol+1)] = sv

		# extract the electrode-molecule interactions
		for ielec in range(self.nelec):
				v  = h[np.ix_(self.electrodes[str(ielec)].index_orb[0],self.central_region.index_orb)]
				sv = s[np.ix_(self.electrodes[str(ielec)].index_orb[0],self.central_region.index_orb)]
				self.electrodes[str(ielec)].vmol = v
				self.electrodes[str(ielec)].smol = sv

		# extract the hamiltonians of the electrodes
		for ielec in range(self.nelec):

			# unit cell
			self.electrodes[str(ielec)].h0 = h[np.ix_(self.electrodes[str(ielec)].index_orb[1],self.electrodes[str(ielec)].index_orb[1])]
			self.electrodes[str(ielec)].s0 = s[np.ix_(self.electrodes[str(ielec)].index_orb[1],self.electrodes[str(ielec)].index_orb[1])]

			# interaction between neighbouting cell
			self.electrodes[str(ielec)].hi = h[np.ix_(self.electrodes[str(ielec)].index_orb[0],self.electrodes[str(ielec)].index_orb[1])]
			self.electrodes[str(ielec)].si = s[np.ix_(self.electrodes[str(ielec)].index_orb[0],self.electrodes[str(ielec)].index_orb[1])]

		# print the molecular orbitals
		if print_mo:

			if self.nmol>1:
				for imol in range(self.nmol):
					w,u = np.linalg.eigh(self.molecules[str(imol)].h0)
					fname = 'mo_mol%d.dat' %imol
					print_mo_mopac(w,u,self.molecules[str(imol)].s0,self.molecules[str(imol)].xyz,fname,self.hkl_param_file)

			w,u = np.linalg.eigh(self.central_region.h0)
			print_mo_mopac(w,u,self.central_region.s0,self.central_region.xyz,'mo_central.dat',self.hkl_param_file)

		# clean the file system
		if clean:
			os.remove(hkl_temp)





