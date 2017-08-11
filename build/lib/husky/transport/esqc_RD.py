#
# This code is based on the following references
#
# [1] Electronic transport calculations for self-assembled monolayers of 1,4-phenylene diisocyanide on Au111 contacts
#      Robert Dahlke and Ulrich Schollwock Phys Rev B 69 085324 2004
#
# [2] Numerical calculations for electronic transport through molecular systems Robert Dahlke PHD THESIS 2004
#
#


import numpy as np
import scipy.linalg as spla
import sys
import warnings


np.set_printoptions(precision=3)
np.set_printoptions(suppress=False)
np.set_printoptions(linewidth=160)

class ESQCsolverRD(object):

	def __init__(self,junction,energies=[]):

		self.junction = junction
		self.energies = energies
		self.force = False

	def set_energy_range(self,emin=-10,emax=10,nE=101):
		self.energies = np.linspace(emin,emax,nE)


	# diagonalize the hamiltonian of each electrode unit cell
	def diagonalize_electrode(self,elec):

		# diagonalize the hamiltonian
		l,u = spla.eigh(elec.h0,b=elec.s0)

		# store  the unit cell H/S
		elec.h0 = np.diag(l)
		elec.s0 = np.eye(elec.norb)

		# store the inter-unitcell H/S
		elec.hi = np.dot(np.linalg.inv(u),np.dot(elec.hi,u))
		elec.si = np.dot(np.linalg.inv(u),np.dot(elec.si,u))
		
		#interaction with the molecule
		elec.vmol = np.dot(np.linalg.inv(u),elec.vmol)
		elec.smol = np.dot(np.linalg.inv(u),elec.smol)
		
		print elec.vmol

		elec.udiag = u

	def diagonalize_all_electrodes(self):
		for index, elec in self.junction.electrodes.iteritems():
			if elec.norb>1:
				self.diagonalize_electrode(elec)

	# diagonalize the hamiltonian of the molecule
	def diagonalize_central(self):

		l,u = spla.eigh(self.junction.central_region.h0,b=self.junction.central_region.s0)

		self.junction.central_region.h0 = l
		self.junction.central_region.s0 = np.eye(self.junction.central_region.norb)
		self.junction.central_region.udiag = u

 		for index, elec in self.junction.electrodes.iteritems():
			elec.vmol = np.dot(elec.vmol,u)	
			elec.smol = np.dot(elec.smol,u)	



	# get the propagator pf the electrode
	def get_electrode_propagator(self,ielec=0,E=0):

		# shortcut
		elec = self.junction.electrodes[str(ielec)]

		# build up the propagator 
		size = elec.norb
		I,Z = np.eye(size), np.zeros((size,size))
		mI = -np.dot(np.linalg.inv(elec.hi.T-E*elec.si.T),elec.hi-E*elec.si)
		p11 = -np.dot(np.linalg.inv(elec.hi.T-E*elec.si.T),elec.h0-E*elec.s0)
		elec.Pd = np.bmat('Z,I;mI,p11')

		# diagonalize the propagator
		L,U = np.linalg.eig(elec.Pd)

		# sort the +/- eigenvaues
		invL = 1./L[:,np.newaxis]
		m = np.absolute(L-invL)**2

		# get the indexes of the pairs
		ips,ims = np.where(m<1E-16)
		unique_index = np.where(ips<ims)
		ips = ips[unique_index]
		ims = ims[unique_index]
		index = np.concatenate([ips,ims])

		if len(index) != 2*size:
			print ''
			print elec.Pd
			print L
			print 1./(L)
			print L*np.conj(L)
			print index
			print m
			print ''
			sys.exit()		

		# store the data
		elec.L = L[index]
		elec.U = U[:,index]
		
		# store the propagating channels
		elec.prop_channel = np.abs(np.absolute(L[index]) - 1.0) < 1E-6 



	# get the total scattering matrix
	def get_total_matrix(self,ielec1=0,ielec2=1,E=0):

		# extract matrices for electrode 1
		h1 = self.junction.electrodes[str(ielec1)].hi.T-E*self.junction.electrodes[str(ielec1)].si.T
		M1 = self.junction.electrodes[str(ielec1)].h0-E*self.junction.electrodes[str(ielec1)].s0
		nelec1 = self.junction.electrodes[str(ielec1)].norb

		# extract matrices for electrode 2
		h2 = self.junction.electrodes[str(ielec2)].hi.T-E*self.junction.electrodes[str(ielec2)].si.T
		M2 = self.junction.electrodes[str(ielec2)].h0-E*self.junction.electrodes[str(ielec2)].s0
		nelec2 = self.junction.electrodes[str(ielec2)].norb

		# extract matrices for electrode/molecule coupling
		vmol1 = self.junction.electrodes[str(ielec1)].vmol-E*self.junction.electrodes[str(ielec1)].smol
		vmol2 = self.junction.electrodes[str(ielec2)].vmol-E*self.junction.electrodes[str(ielec2)].smol
		
		# extract matrix for molecule
		Mmol = self.junction.central_region.h0-E*self.junction.central_region.s0
		nmol = self.junction.central_region.norb

		# extended central system
		M0 = spla.block_diag(M1,M2,Mmol)
		M0[:nelec1,nelec1+nelec2:] = vmol1
		M0[nelec1+nelec2:,:nelec1] = vmol1.T
		M0[nelec1:nelec1+nelec2,nelec1+nelec2:] = vmol2
		M0[nelec1+nelec2:,nelec1:nelec1+nelec2] = vmol2.T
		nxt = len(M0)

		# coupling between the electrodes and the extend central 
		tau1 = np.zeros((nelec1,nxt))
		tau1[:,:nelec1] = h1

		# coupling between the electrodes and the extend central 
		tau2 = np.zeros((nelec2,nxt))
		tau2[:,nelec1:nelec1+nelec2] = h2

		#Zero matrices we need for assembly
		z1 = np.zeros((nelec1,nelec2))
		z2 = np.zeros((nelec2,nelec1))
		z3 = np.zeros((nxt,nelec1))
		z4 = np.zeros((nxt,nelec2))

		self.Meff = np.vstack( (
		np.hstack( ( h1, M1,     z1, z1,     tau1 )  ) ,
		np.hstack( ( z2, z2,     h2, M2,     tau2 )  ) ,
		np.hstack( ( z3, tau1.T, z4, tau2.T, M0   )  ) ) )

	# transform the total matrix
	def get_scattering_matrix(self,ielec1=0,ielec2=1):

		nelec1 = self.junction.electrodes['0'].norb
		nelec2 = self.junction.electrodes['1'].norb
		nmol = self.junction.central_region.norb+nelec1+nelec2
		ntot = 2*nelec1+2*nelec2+nmol

		Imol = np.eye(nmol)
		U = spla.block_diag(self.junction.electrodes[str(ielec1)].U,self.junction.electrodes[str(ielec2)].U,Imol)
		
		# indexes
		index_out = range(nelec1,2*nelec1)+range(2*nelec1+nelec2,2*nelec1+2*nelec2+nmol)
		index_in = range(nelec1)+range(2*nelec1,2*nelec1+nelec2)

		# transform the matrix
		S = np.dot(self.Meff,U)

		# in/ou matrices
		Mout = S[:,index_out]
		Min  = S[:,index_in]
		
		# scattering matrix
		self.Scat = -np.dot(np.linalg.inv(Mout),Min)
		#print Mout.shape, Min.shape, self.Scat.shape

	# transmission
	def esqc_transmission(self,ielec1=0,ielec2=1):

		# declare the te
		te = np.zeros(len(self.energies))		



		# size of the electrodes
		nelec1 = self.junction.electrodes[str(ielec1)].norb
		nelec2 = self.junction.electrodes[str(ielec2)].norb

		# loop over the energies
		for iE,e in enumerate(self.energies):

			# get the electrode propagators
			self.get_electrode_propagator(ielec=ielec1,E=e)
			self.get_electrode_propagator(ielec=ielec2,E=e)

			# form the total matrix
			self.get_total_matrix(ielec1=ielec1,ielec2=ielec2,E=e)

			# get the scattering matrix
			self.get_scattering_matrix(ielec1=ielec1,ielec2=ielec2)
			
			# form the mask matrix of propagating channels 
			prop_matrix = self.junction.electrodes[str(ielec1)].prop_channel[:nelec1] * self.junction.electrodes[str(ielec2)].prop_channel[nelec2:,np.newaxis]

			# compute S12 and the transmission coefficient
			te[iE] = np.linalg.norm( self.Scat[:nelec1,nelec1:nelec1+nelec2] * prop_matrix )
			

		return te




		
