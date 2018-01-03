#
# This code is based on the following references
#
#[1] C Herman JACS 132,3682,2010 Organic Radicalas spin filters
#[2] Understanding transport through molecular orbitals Johannes S. Seldenthuis and Joseph M. Thijssen
#[3] S Data 

import numpy as np
import time

try:
	from tqdm import tqdm
except:
	def tqdm(a):
		return a

class NEGFsolver(object):

	def __init__(self,junction,energies=[],wbl=True):

		self.junction = junction
		self.energies = energies
		self.ldos = [1.0,1.0]
		self.wide_band_limit = wbl

		self.Vbias  = 0.001
		self.kT     = 0.025
		self.Lambda = 0.00

	# set energy range for the calculation
	def set_energy_range(self,emin=-10,emax=10,nE=101):
		self.energies = np.linspace(emin,emax,nE)

	# set ldos for WBL calculation
	def set_local_dos_electrode(self,ldos,ielec=None):
		if ielec==None:
			self.ldos = [ldos,ldos]
		else:
			self.ldos[ielec] = ldos

	# set the logical for WBL calculations
	def set_wide_band_limit(self,wbl):
		self.wide_band_limit = wbl

	# compute the self energy of the electrodes in the WBL
	def self_energy_elec(self,E,ielec):
		
		# dos of the elecrode
		if self.wide_band_limit:
			g = 1.0j*self.ldos[ielec]*np.eye(self.junction.electrodes[str(ielec)].norb)
		else:
			g = self.junction.electrodes[str(ielec)].interpolate_surface_green_function(E)

		# get the coupling
		Sigma = (E*self.junction.electrodes[str(ielec)].smol - self.junction.electrodes[str(ielec)].vmol)

		# form the self energy
		G = np.dot(Sigma.T,np.dot(g,Sigma))
		return G

	# compute the retarded Green function
	def green_function(self,E,G1,G2):
		return np.linalg.inv(E*self.junction.central_region.s0-self.junction.central_region.h0 - (G1+G2))
		
	# compute the transmission from [1]
	def compute_transmission(self,ielec1=0,ielec2=1):

		self.print_intro()
		te = np.zeros(len(self.energies))
		t0 = time.time()
		for iE,e in enumerate(tqdm(self.energies,ncols=50,desc='   ',leave=True)):

			G1 = self.self_energy_elec(e,ielec1)
			G2 = self.self_energy_elec(e,ielec2)
			Gr = self.green_function(e,G1,G2)

			m1 = np.dot(2.0*G1.imag,Gr)
			m2 = np.dot(2.0*G2.imag,Gr.conj())

			temp = np.trace(np.dot(m1,m2))
			te[iE] = np.sqrt((temp*np.conj(temp)).real)

		print(' - Transmission computed in %1.3f sec. ' %(time.time()-t0))
		return te

	# print the detail of the calculation
	def print_intro(self):
		print('\n')
		print('='*40)
		print('= Non Equilibrium Green Function Solver')
		print('= N. Renaud 2017')
		print('='*40)
		print('')
		print(' - Energy range % 1.3f % 1.3f %d points' %(self.energies[0],self.energies[-1],len(self.energies)))
		print(' - Wide Band Limit Approximation ', self.wide_band_limit)
		print(' - Compute Transmission')

	# compute the transmission from [2]
	def compute_transmission_sos(self,ielec1=0,ielec2=1):

		# decalre arrays
		te = np.zeros(len(self.energies))

		# Construct the retarded self-energy
		gL = (self.junction.electrodes[str(ielec1)].vmol[0,:].T)
		gR = (self.junction.electrodes[str(ielec2)].vmol[0,:].T)
		Sigma = -0.5j * (np.outer(gL, gL) + np.outer(gR, gR))

		# Diagonalize the (inverse) Green's function
		w, C = np.linalg.eig(self.junction.central_region.h0 + Sigma)

		# Sort the eigenvalues and eigenvectors
		ndx = np.argsort(np.real(w))
		eps = np.real(w[ndx])
		Gamma = -np.imag(w[ndx])
		C = C[:, ndx]

		# Normalize the eigenvectors
		for i in range(len(w)):
		    C[:, i] /= np.sqrt(np.dot(C[:, i], C[:, i]))

		# Calculate the effective couplings
		nuL = np.dot(C.T, gL)
		nuR = np.dot(C.T, gR)
		nu = nuL * nuR

		# Calculate the orbital transmissions
		for iE,e in enumerate(self.energies):
		    t = nu / (e - eps + 1j * Gamma)
		    te[iE] = np.abs(np.sum(t))**2

		return te


	# compute the transmission with dephasing [3]
	def compute_transmission_dephasing(self,imol='allMol',ielec1=0,ielec2=1,Lambda=0,iterMax=100,tol=1E-12,zplus=1.0j*1E-14):


		# iter  : max maximum iteration for the self consistent calculation of the Green function
		# zplus : infinitesimal for the Green function
		# tol   : tolerance for the self consistenant Green function 


		# decalre arrays
		te = np.zeros(len(self.energies))
		self.Lambda = Lambda

		# Construct the retarded self-energy
		gL = (self.junction.electrodes[str(ielec1)].vmol[0,:].T)
		gR = (self.junction.electrodes[str(ielec2)].vmol[0,:].T)

		sig1 = 0.5j * np.outer(gL, gL) 
		sig2 = 0.5j * np.outer(gR, gR)

		# extract the Hamiltonain
		H = self.junction.central_region.h0
		Np = len(H)		
		
		# fermi function
		f1 = 1./(1.+np.exp(-(self.Vbias/2)/self.kT))
		f2 = 1./(1.+np.exp(+(self.Vbias/2)/self.kT))

		# for all the energies
		for k,e in enumerate(self.energies):

			# compute the initial green function		
			G = np.linalg.inv( (e+zplus)*np.eye(Np,Np) - H - sig1 - sig2 )
			change = 1

			# self consistent loop
			nIter = 1
			while (change > tol) and (nIter < iterMax): 

				#compute the interaction self energy
				sigp = np.diag(self.Lambda*np.diag(G))

				# new green function
				S = np.linalg.inv( (e+zplus)*np.eye(Np,Np) - H - sig1 - sig2 -sigp )

				# change
				change = np.sum( np.sum(np.abs(G-S)))/np.sum(np.sum(np.abs(G)+np.abs(S)))


				# update the Green function
				G = 0.5*(G+S)

				# next iteration
				nIter += 1

				if (nIter >= iterMax):
					print(" Green function not converged : error = %e\n" %change) 

			
			# compute the spectral function
			G = S
			A = 1.0j*(G-np.conj(G.T))
			M = self.Lambda*np.multiply(G,np.conj(G))


			# electron density
			gamp = 1.0j*(sigp-np.conj(sigp.T))
			gam1 = 1.0j*(sig1-np.conj(sig1.T))
			gam2 = 1.0j*(sig2-np.conj(sig2.T))			
			gama = gam1+gam2+gamp

			sigin1 = f1*gam1
			sigin2 = f2*gam2

			Sin = np.dot(np.dot(G,(sigin1+sigin2)),np.conj(G.T))
			n = np.array(np.dot(np.linalg.inv(np.eye(Np,Np)-M),np.diag(Sin)))
			siginp = Lambda*np.diag(n)

			# compute the correlation function Gn
			Gn = np.dot(G,np.dot((sigin1+sigin2+siginp),np.conj(G.T)))

			# transmission
			te[k] = np.real(1/(f2-f1)*(np.trace(np.dot(gam1,Gn)) - np.trace(np.dot(sigin1,A))))
			
		return te

