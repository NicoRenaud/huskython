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
import time
import itertools

try:
	from tqdm import tqdm
except:
	def tqdm(a):
		return a

np.set_printoptions(precision=3)
np.set_printoptions(suppress=False)
np.set_printoptions(linewidth=160)

class CondChannel(object):

	def __init__(self):

		self.type = None
		self.index = None
		self.eigen = None
		self.inveigenv = None
		self.abs = None
		self.phase = None

class ESQCsolver(object):

	def __init__(self,junction,energies=[]):

		self.junction = junction
		self.energies = energies
		self.eps_prop = 1E-6
		self.force = False
		self.debug = False

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
	def get_electrode_propagator(self,ielec=0,E=0):

		# shortcut
		elec = self.junction.electrodes[str(ielec)]

		# build up the propagator 
		size = elec.norb
		I,Z = np.eye(size), np.zeros((size,size))

		hi = elec.hi.T-E*elec.si.T  #+ 1E-16*np.random.rand(elec.norb,elec.norb)
		if np.linalg.cond(hi) < 1./np.finfo(hi.dtype).eps: 
			invHi = spla.inv(hi)
		else:
			print '\n -> Warning : Interaction matrix between electrode cell ill conditionned'
			U,s,Vh = spla.svd(hi)
			invHi = np.dot(U,np.dot(np.diag(1./s),Vh))
			if not self.force:
				sys.exit()

		mI = -np.dot(invHi,elec.hi-E*elec.si)
		p11 = -np.dot(invHi,elec.h0-E*elec.s0)
		elec.Pd = np.bmat('Z,I;mI,p11')


		# diagonalize the propagator
		L,U = np.linalg.eig(elec.Pd)

		# sort the in/out channels
		channels,sort_index,prop = self.get_channels(L,eps=self.eps_prop,_force_=self.force,_debug_=self.debug)

		# store the data
		elec.L = L[sort_index]
		elec.U = U[:,sort_index]
		elec.channels = channels
		elec.index_in = range(elec.norb)
		elec.index_out = range(elec.norb,elec.norb*2)

		# store the propagating channels
		elec.prop_channel = np.array(prop)

		# store the number of prop channel
		elec.nprop_channel.append(np.sum(elec.prop_channel))


	# sort the in/out channels of the electrodes
	@staticmethod
	def get_channels(L,eps=1E-6,_force_=False,_debug_=False):

		# make copies and initiate
		l = np.copy(L)
		n = len(l)
		channels = []
		miss = range(n)

		# find the pairs that respect lp=1./lm^*
		for i in range(n):

			if np.isnan(l[i]) :
				continue

			for j in range(n):

				if i == j or np.isnan(l[j]):
					continue

				# if we find a pair
				rdiff1 = np.abs(l[i] - 1./l[j])/np.abs(l[i])
				rdiff2 = np.abs(l[j] - 1./l[i])/np.abs(l[j])

				if rdiff1 < eps and rdiff2 < eps:

					# create the channels
					c = CondChannel()
					c.index = [i,j]
					c.eigen = [l[i],l[j]]
					c.inveigen = [1./c.eigen[0],1./c.eigen[1]]
					c.abs = [np.abs(l[i]),np.abs(l[j])]

					if np.abs(c.abs[0]-1) < eps and np.abs(c.abs[1]-1) < eps:
						c.type = 'propagating'
					else:
						c.type = 'evanescent'

					c.phase = [np.angle(l[i]),np.angle(l[j])]

					# append the channel to the list
					channels.append(c)

					# print for debug
					if _debug_:

						print '\n---- Channel %d ----\n' %(i+1)
						print ' -> %s' %c.type
						print ' -> index : %d    %d'    %(c.index[0],      c.index[1])
						print ' -> eigen : {:.3E} {:.3E}'.format(c.eigen[0],      c.eigen[1])
						print ' -> inveg : {:.3E} {:.3E}'.format(c.inveigen[0],   c.inveigen[1])
						print ' -> norm  : %1.3E %1.3E' %(c.abs[0],        c.abs[1])
						print ' -> phase : %1.3f %1.3f' %(c.phase[0],      c.phase[1])
						#print ' -> diff  : %1.3E %1.3E' %(np.abs(c.eigen[0]-1./c.eigen[1]),np.abs(c.eigen[1]-1./c.eigen[0]))
						print ' -> rdiff : %1.3E %1.3E' %( rdiff1,rdiff2)

					# delete the index from the list
					miss.remove(i)
					miss.remove(j)


					# burn the two values
					l[i] = l[j] = None

					break

		# fix the channels if required
		nmiss = len(miss)
		if nmiss > 0:

			print ''
			print '='*60
			print '==  Warning : %d/%d  conduction channels are corrputed' %(nmiss,len(L))
			print '='*60

			l = np.array(L[miss])
			mat = np.abs(l-1./l[:,np.newaxis])
			mat += mat.T
			#mat += float('inf')*np.eye(nmiss)

			for i in range(nmiss/2):

				row_min,col_min = np.unravel_index(mat.argmin(),mat.shape)

				# create the channels
				c = CondChannel()
				c.type = 'evanescent'
				c.index = [col_min,row_min]
				c.eigen = [l[c.index[0]],l[c.index[1]]]
				c.inveigen = [1./c.eigen[0],1./c.eigen[1]]
				c.abs = [np.abs(c.eigen[0]),np.abs(c.eigen[1])]
				c.phase = [np.angle(c.eigen[0]),np.angle(c.eigen[1])]

				# append the channel
				channels.append(c)

				# burns the values
				mat[row_min,:] = float('inf')
				mat[col_min,:] = float('inf')
				mat[:,row_min] = float('inf')
				mat[:,col_min] = float('inf')

				# print for debug
				if _debug_:

					print '\n---- Fixed Channel %d ----\n' %(i+1)
					print ' -> %s' %c.type
					print ' -> index : %d    %d'    %(c.index[0],      c.index[1])
					print ' -> eigen : {:.3E} {:.3E}'.format(c.eigen[0],      c.eigen[1])
					print ' -> inveg : {:.3E} {:.3E}'.format(c.inveigen[0],   c.inveigen[1])
					print ' -> norm  : %1.3E %1.3E' %(c.abs[0],        c.abs[1])
					print ' -> phase : %1.3f %1.3f' %(c.phase[0],      c.phase[1])
					print ' -> rdiff  : % 1.3E %1.3E' %( np.abs(c.eigen[0]-1./c.eigen[1])/c.abs[0],np.abs(c.eigen[1]-1./c.eigen[0])/c.abs[1])

			if not _force_:
				sys.exit()


		# get the index list
		index_in = []
		index_out = []
		prop = []
		for c in channels:
			if c.type == 'evanescent':
				index_in.append(c.index[0])
				index_out.append(c.index[1])
				prop.append(0)
			elif c.phase[0] == 0 and c.phase[1] == 0:
				index_in.append(c.index[0])
				index_out.append(c.index[1])
				prop.append(1)
			elif c.phase[0] == np.pi and c.phase[1] == np.pi:
				index_in.append(c.index[0])
				index_out.append(c.index[1])
				prop.append(1)
			else:
				a = np.argwhere(np.array(c.phase)>0)[0][0]
				b = np.argwhere(np.array(c.phase)<0)[0][0]
				index_in.append(c.index[a])
				index_out.append(c.index[b])
				prop.append(1)

		index = index_in + index_out
		prop = prop + prop

		return channels,index,prop


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

		nelec1 = self.junction.electrodes[str(ielec1)].norb
		nelec2 = self.junction.electrodes[str(ielec2)].norb

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
		R = 1E-16*np.eye(len(Mout))

		# scattering matrix
		self.Scat = -np.dot(np.linalg.inv(Mout),Min)
		

	# transmission
	def compute_transmission(self,ielec1=0,ielec2=1):

		# declare the te
		te = np.zeros(len(self.energies))		

		# size of the electrodes
		nelec1 = self.junction.electrodes[str(ielec1)].norb
		nelec2 = self.junction.electrodes[str(ielec2)].norb

		# print the intro
		self.print_intro()

		t0 = time.time()
		
		# loop over the energies
		for iE,e in enumerate(tqdm(self.energies,ncols=50,desc='   ',leave=True)):

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
			te[iE] = np.linalg.norm( self.Scat[:nelec1,nelec1:nelec1+nelec2] * prop_matrix )**2
			
			
			
		print ' - Transmission computed in %1.3f sec. ' %(time.time()-t0)
		return te

	# print the detail of the calculation
	def print_intro(self):
		print '\n'
		print '='*40
		print '= Ellastic Scattering Quantum Chemistry Solver'
		print '= N. Renaud 2017'
		print '='*40
		print ''
		print ' - Energy range % 1.3f % 1.3f %d points' %(self.energies[0],self.energies[-1],len(self.energies))
		print ' - Compute Transmission'



		
