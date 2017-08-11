#
# This code is based on the following references
#
#[1] C Herman JACS 132,3682,2010 Organic Radicalas spin filters
#[2] Understanding transport through molecular orbitals Johannes S. Seldenthuis and Joseph M. Thijssen
#[3] S Data 

import numpy as np

class negf(object):

	def __init__(self,sys,energies=[]):

		self.ElectronicSystem = sys
		self.energies = energies
		self.ldos = 1.0

		self.Vbias  = 0.001
		self.kT     = 0.025
		self.Lambda = 0.00

	def set_energy_range(self,emin=-10,emax=10,nE=101):
		self.energies = np.linspace(emin,emax,nE)

	# compute the self energy of the electrodes in the WBL
	def self_energy_elec(self,E,ielec,imol):
		
		# dos of the elecrode
		g = 1.0j*self.ldos*np.eye(self.ElectronicSystem.norb_elec[ielec])

		# get the coupling
		key = '%d->%s' %(ielec,str(imol))
		Sigma = E*self.ElectronicSystem.SelecMol[key] - self.ElectronicSystem.VelecMol[key]

		# form the self energy
		G = np.dot(Sigma.T,np.dot(g,Sigma))
		return -2*G.imag

	# compute the retarded Green function
	def green_function(self,E,G1,G2,imol):
		imol = str(imol)
		return np.linalg.inv(E*self.ElectronicSystem.Smol[imol]-self.ElectronicSystem.Hmol[imol]+ 0.5j*(G1+G2))
		
	# compute the transmission from [1]
	def negf_transmission(self,imol='allMol',ielec1=0,ielec2=1):

		self.te = np.zeros(len(self.energies))
		
		for iE,e in enumerate(self.energies):

			G1 = self.self_energy_elec(e,ielec1,imol)
			G2 = self.self_energy_elec(e,ielec2,imol)
			Gr = self.green_function(e,G1,G2,imol)

			m1 = np.dot(G1,Gr)
			m2 = np.dot(G2,Gr.conj())
			self.te[iE] = np.trace(np.dot(m1,m2).real)

	# compute the transmission from [2]
	def negf_transmission_sos(self,imol='allMol',ielec1=0,ielec2=1):

		# decalre arrays
		self.te = np.zeros(len(self.energies))

		# Construct the retarded self-energy
		key1 = '%d->%s' %(ielec1,str(imol))
		key2 = '%d->%s' %(ielec2,str(imol))
		gL = (self.ElectronicSystem.VelecMol[key1][0,:].T)
		gR = (self.ElectronicSystem.VelecMol[key2][0,:].T)
		Sigma = -0.5j * (np.outer(gL, gL) + np.outer(gR, gR))

		# Diagonalize the (inverse) Green's function
		w, C = np.linalg.eig(self.ElectronicSystem.Hmol[str(imol)] + Sigma)

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
		    self.te[iE] = np.abs(np.sum(t))**2


	# compute the transmission with dephasing [3]
	def negf_transmission_dephasing(self,imol='allMol',ielec1=0,ielec2=1,Lambda=0,iterMax=100,tol=1E-12,zplus=1.0j*1E-14):


		# iter  : max maximum iteration for the self consistent calculation of the Green function
		# zplus : infinitesimal for the Green function
		# tol   : tolerance for the self consistenant Green function 


		# decalre arrays
		self.te = np.zeros(len(self.energies))
		self.Lambda = Lambda

		# Construct the retarded self-energy
		key1 = '%d->%s' %(ielec1,str(imol))
		key2 = '%d->%s' %(ielec2,str(imol))
		gL = (self.ElectronicSystem.VelecMol[key1][0,:].T)
		gR = (self.ElectronicSystem.VelecMol[key2][0,:].T)

		sig1 = 0.5j * np.outer(gL, gL) 
		sig2 = 0.5j * np.outer(gR, gR)

		# extract the Hamiltonain
		H = self.ElectronicSystem.Hmol[str(imol)]
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
					print " Green function not converged : error = %e\n" %change 

			
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
			self.te[k] = np.real(1/(f2-f1)*(np.trace(np.dot(gam1,Gn)) - np.trace(np.dot(sigin1,A))))
			

