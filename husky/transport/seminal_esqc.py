#
# This code is based on the following references
#
# [1] Electronic Transmission Through a Single Impurity in a Multi-configuration Scattering Matrix Approach
#     M Portaix C Joachim  L. Grill and C. Joachim (eds.), Imaging and Manipulating Molecular Orbitals,
#     Advances in Atom and Single Molecule Machines, DOI: 10.1007/978-3-642-38809-5_11, Springer-Verlag Berlin Heidelberg 2013
#
# [2] Intramolecular circuits connected to N electrodes using a scattering matrix approach
#     S Ami C joachim PHYSICAL REVIEW B, VOLUME 65, 155419
#
# [3] Electronic tranmission for single impurity problem P Sautet C Joachim PRB 38,12238, 1988
#
# [4] Conductance of molecular wires connected or bonded in parallel M. Magoga and C. Joachim PRB 59,16011, 1999 
#
# 		=============== WARNING =============== ================================
#			THIS CLASS ONLY WORK FOR  ELECTRODES WITH UNIT CELLS OF SIZE 1
#			DONT KNOW WHY THOUGH AND IT SHOULD BE FIXABLE
# 		=============== ------- =============== ================================
#
#


import numpy as np
import scipy.linalg as spla
import sys
import warnings


np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

class ESQCsolver(object):

	def __init__(self,junction,energies=[]):

		self.junction = junction
		self.energies = energies
		self.force = False

		for index, elec in self.junction.electrodes.iteritems():
			if elec.norb>1:
				print '================================'
				print '==    Warning     '
				print '==    this class only works for '
				print '==    electrodes with unit cell '
				print '==    of size 1' 
				print '================================'

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
	def get_electrode_propagator(self,ielec=0,E=0,n=0,return_inverse=False):


		# form the electrode propagator
		#
		#		| p11  -II |
		# P = 	| 		   |
		# 		|  II    0 |
		#
		# with p11 = -hi^{-1} . (h0-EI)
		#
 		# compute the U matrix that transform the electrode propagator in a Identity matrix
		# 
		#		 |  e^{iKn}           e^{-iKn}    |
		# Umat = |								  |
		#		 | e^{iK(n-1)}      e^{-iK(n-1)}  |
		#
		# where K in a diagonal matrix containing the arccos of the eigenvalues of p11
		#
		# also returns in the electrode instance a generating function 
		# that can compute Umat for different values of n (evaluate_umat)
		#

		print_error = True
		elec = self.junction.electrodes[str(ielec)]


		# build up the propagator 
		size = elec.norb
		I,mI,Z = np.eye(size), -np.eye(size), np.zeros((size,size))
		p11 = -np.dot(np.linalg.inv(elec.hi.T-E*elec.si.T),elec.h0-E*elec.s0)
		

		# diagonalize p11, create the diag matrix 
		lp11,Xp11 = np.linalg.eig(p11)
		Xp11 = np.linalg.inv(Xp11)
		

		# check for problematic channels
		elec.prop_channel = [True]*elec.norb
		for i in range(elec.norb):
			if np.abs(lp11[i]) > 2:
				lp11[i] = 0.0
				elec.prop_channel[i] = False

		# compute the K matrix
		K = np.arccos(0.5*lp11)	

		# store the matrices
		lp11 = np.diag(lp11)
		elec.Pd = np.bmat('lp11,mI;I,Z')
		

		# clean the Xp matrix and stor it
		Xp11[np.ix_(np.invert(elec.prop_channel),np.invert(elec.prop_channel))] = 0
		Xp11[np.invert(elec.prop_channel),np.invert(elec.prop_channel)] = 1
		elec.Xmat = np.kron(np.eye(2),Xp11)
		

		# create the U(p,p-1) matrix
		def up(p=1):

			# create the matrices
			U00, U01 = np.diag(np.exp(1j*K*p)), np.diag(np.exp(-1j*K*p))
			U10, U11 = np.diag(np.exp(1j*K*(p-1))), np.diag(np.exp(-1j*K*(p-1)))

			# remove the non conducting channels
			U00[np.ix_(np.invert(elec.prop_channel),np.invert(elec.prop_channel))] = 0
			U01[np.ix_(np.invert(elec.prop_channel),np.invert(elec.prop_channel))] = 0
			U10[np.ix_(np.invert(elec.prop_channel),np.invert(elec.prop_channel))] = 0
			U11[np.ix_(np.invert(elec.prop_channel),np.invert(elec.prop_channel))] = 0

			# return the total matrix
			return np.bmat('U00,U01;U10,U11')

		elec.evaluate_umat = up
		elec.Umat = up(n)

		if return_inverse:
			dbl_prop_channel = elec.prop_channel+elec.prop_channel
			elec.Umat[np.ix_(dbl_prop_channel,dbl_prop_channel)] = np.linalg.inv(elec.Umat[np.ix_(dbl_prop_channel,dbl_prop_channel)])
		
			

			
	def get_central_propagator(self,ielec1=0,ielec2=1,E=0):

		elec1 = self.junction.electrodes[str(ielec1)]
		elec2 = self.junction.electrodes[str(ielec2)]
		mol   = self.junction.central_region

		# inverse of the Hmol-EI
		eps = 1E-12
		invhm = - np.diag(1./(mol.h0-E+1j*eps))


		###############################################################################
		# first propagator 
		# 
		#  |c-1|   |  -(hi.T)^-1 x K      -(hi.T)^-1 x M 		| |c 0|
		#  |   | = |											| |   |
		#  |c 0|   |        II                  0				| |c 1|
		#
		#   with K = (elec.h0-EI) + elec1.hi x (- ( mol.h0-EI)^-1 ) elec1.hi.T
		#        M = elec1.hi x (- ( mol.h0-EI)^-1 ) elec2.hi.T
		#
		###############################################################################

		I = np.eye(elec1.norb)
		Z = np.zeros((elec1.norb,elec2.norb))

		K = elec1.h0-E*elec1.s0 + np.dot(elec1.vmol-E*elec1.smol,np.dot(invhm,elec1.vmol.T-E*elec1.smol.T))
		M = np.dot(elec1.vmol-E*elec1.smol,np.dot(invhm,elec2.vmol.T-E*elec2.smol.T))

		A = - np.dot(np.linalg.inv(elec1.hi.T-E*elec1.si.T),K) 
		B = - np.dot(np.linalg.inv(elec1.hi.T-E*elec1.si.T),M)

		Pm1 = np.bmat('A,B;I,Z')	


		###############################################################################
		# second propagator 
		# 
		#  |c 0|   |  - M^-1 x K            - M^-1 x hi 		| |c 1|
		#  |   | = |											| |   |
		#  |c 1|   |        II                  0				| |c 2|
		#
		#   with K = (elec.h0-EI) + elec2.hi x (- ( mol.h0-EI)^-1 ) elec2.hi.T
		#        M = elec2.hi x (- ( mol.h0-EI)^-1 ) elec1.hi.T
		#
		###############################################################################

		I = np.eye(elec2.norb)
		Z = np.zeros((elec2.norb,elec2.norb))
		K = elec2.h0-E*elec2.s0 + np.dot(elec2.vmol-E*elec2.smol,np.dot(invhm,elec2.vmol.T-E*elec2.smol.T))
		M = np.dot(elec2.vmol-E*elec2.smol,np.dot(invhm,elec1.vmol.T-E*elec1.smol.T))
		Mm1 = np.linalg.inv(M+1E-12*np.eye(elec2.norb))
	

		A = - np.dot(Mm1,K) 
		B = - np.dot(Mm1,elec2.hi-E*elec2.si)

		P1 = np.bmat('A,B;I,Z')	


		###############################################################################
		# final propagator
		#
		# 							M_eff = Pm1 x P1
		#
		###############################################################################

		self.Meff = np.dot(Pm1,P1)

		if 0:
			print 'Meff\n', self.Meff
			print 'Pm1\n', Pm1
			print 'P1\n', P1
			print 'Mm1\n', Mm1	
			print 'M\n', M
			print 'K\n',K


	# compute the transmission
	def esqc_transmission(self,ielec1=0,ielec2=1,diag_elec=False):

		te = np.zeros(len(self.energies))

		#self.assemble_test()

		#diagonalize the electrode subspace
		if diag_elec:
			for index, elec in self.junction.electrodes.iteritems():
				if elec.norb>1:
					self.diagonalize_electrode(elec)

		# diagonalize the central region
		self.diagonalize_central()

		print 'he\n', self.junction.electrodes['0'].h0
		print 'hi\n', self.junction.electrodes['0'].hi
		print 'VL\n', self.junction.electrodes['0'].vmol
		print 'H\n', self.junction.central_region.h0
		print 'VR\n', self.junction.electrodes['1'].vmol
		print 'he\n', self.junction.electrodes['1'].h0
		print 'hi\n', self.junction.electrodes['1'].hi


		# loop over the energies
		for iE,e in enumerate(self.energies):



			###########################
			#
			# first elecrode
			#
			###########################

			if 1:
				# get the electrode propagators
				self.get_electrode_propagator(ielec=ielec1,E=e,n=-2,return_inverse=False)

				# assemble the  matrix : U(0)-1 . X^-1 = (X.U)^-1
				mtemp_1 = np.dot( self.junction.electrodes[str(ielec1)].Xmat, self.junction.electrodes[str(ielec1)].Umat)
				mtemp_1 = np.linalg.inv(mtemp_1+1E-12*np.eye(len(mtemp_1)))
				
			else:

				# get the electrode propagator
				self.get_electrode_propagator(ielec=ielec1,E=e,n=0,return_inverse=True)

				# assemble the  matrix : U(0)-1 . X^-1 for elec1 
				# !! Note that we return the inverse of Umat for elec1
				mtemp_1 = np.dot(self.junction.electrodes[str(ielec1)].Umat, np.linalg.inv(self.junction.electrodes[str(ielec1)].Xmat))


			###########################
			#
			# Second elecrode
			#
			###########################

			# get the electrode propagators
			self.get_electrode_propagator(ielec=ielec2,E=e,n=-1,return_inverse=False)


			# assemble the matrix X . U(1) for elec2
			mtemp_2 = np.dot(self.junction.electrodes[str(ielec2)].Xmat,self.junction.electrodes[str(ielec2)].Umat)


			###########################
			#
			# Central part
			#
			###########################

			# get the central propagator
			self.get_central_propagator(ielec1=ielec1,ielec2=ielec2,E=e)

			
			###########################
			#
			# Assemble
			#
			###########################
		
			# assemble the total propagator and extract the part we are interested in
			norb_elec1 = self.junction.electrodes[str(ielec1)].norb
			norb_elec2 = self.junction.electrodes[str(ielec2)].norb
			F = np.dot(mtemp_1,np.dot(self.Meff,mtemp_2))[:norb_elec1,:norb_elec1]
			
			# inverse of the F matrix and TE
			index = self.junction.electrodes[str(ielec1)].prop_channel
			Fm1 = np.linalg.inv(F[np.ix_(index,index)])
			te[iE] = np.linalg.norm(Fm1)
			te[iE] = np.sum(np.sqrt(Fm1*np.conj(Fm1)))

		print 'Umat\n ', self.junction.electrodes['0'].Umat
		print 'Xmat\n ', self.junction.electrodes['0'].Xmat
		print 'F-1\n', Fm1
		print 'Meff\n', self.Meff
		return te



	def assemble_test(self):

		elec1 = self.junction.electrodes['0']
		elec2 = self.junction.electrodes['1']
		mol   = self.junction.central_region

		ntot = 2*elec1.norb+2*elec2.norb+mol.norb


		i = 0
		elec1_unit1_index = range(i,i+elec1.norb)
		i += elec1.norb 

		elec1_unit2_index = range(i,i+elec1.norb)
		i += elec1.norb 

		mol_index = range(i,i+mol.norb)
		i += mol.norb

		elec2_unit1_index = range(i,i+elec2.norb)
		i += elec2.norb

		elec2_unit2_index = range(i,i+elec2.norb)
		i += elec2.norb

		H  = np.zeros((ntot,ntot))
		H[np.ix_(elec1_unit1_index,elec1_unit1_index)] = elec1.h0
		H[np.ix_(elec1_unit2_index,elec1_unit2_index)] = elec1.h0

		H[np.ix_(elec1_unit1_index,elec1_unit2_index)] = elec1.hi
		H[np.ix_(elec1_unit2_index,elec1_unit1_index)] = elec1.hi.T

		H[np.ix_(elec1_unit2_index,mol_index)] = elec1.vmol
		H[np.ix_(mol_index,elec1_unit2_index)] = elec1.vmol.T


		H[np.ix_(elec2_unit1_index,elec2_unit1_index)] = elec2.h0
		H[np.ix_(elec2_unit2_index,elec2_unit2_index)] = elec2.h0

		H[np.ix_(elec2_unit1_index,elec2_unit2_index)] = elec1.hi
		H[np.ix_(elec2_unit2_index,elec2_unit1_index)] = elec1.hi.T

		H[np.ix_(elec2_unit2_index,mol_index)] = elec2.vmol
		H[np.ix_(mol_index,elec1_unit2_index)] = elec2.vmol.T

		H[np.ix_(mol_index,mol_index)] = mol.h0

		l,elec1.udiag = np.linalg.eigh(elec1.h0)
		l,elec2.udiag = np.linalg.eigh(elec2.h0)
		l,mol.udiag = np.linalg.eigh(mol.h0)

		U = spla.block_diag(elec1.udiag,elec1.udiag,mol.udiag,elec2.udiag,elec2.udiag)
		Hd = np.dot(np.linalg.inv(U),np.dot(H,U))

		l,u = np.linalg.eig(H)
		ld,ud = np.linalg.eig(Hd)


		print Hd[np.ix_(elec1_unit2_index,mol_index)]
		print Hd[np.ix_(elec2_unit2_index,mol_index)]

		


