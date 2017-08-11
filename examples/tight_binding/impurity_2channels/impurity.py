from husky.hamiltonian.ElectronicSystem import Junction

from husky.transport.esqc import ESQCsolver
from husky.transport.negf import NEGFsolver

from husky.hamiltonian.model import tight_binding as tb

import matplotlib.pyplot as plt 
import numpy as np


############################################################################
# compute the transmission coefficient of an impurity
# connected to two electrodes that each contains 2 state per unit cell
#
#
#			 - o - o - o - o --- A --- o - o - o - o - 
#			   | x | x | x |  x  |  x  | x | x | x |  
#			 - o - o - o - o --- A --- o - o - o - o -
#
############################################################################

# define the junction 
system = Junction()

# impurity hamiltonian
e,alpha = 2,-2
h0 = np.array([[e,alpha],[alpha,-e]])

# add the molecule hamiltonian
system.add_central_matrix(h0=h0)

# inetractions between electrode layers
beta = -2.0
hi = beta*np.eye(2)


# add the electrodes
system.add_electrode_matrix(h0 = h0, hi = hi, vmol = hi)
system.add_electrode_matrix(h0 = h0, hi = hi, vmol = hi)

# compute the band structure of the right electrode
system.electrodes['0'].band_structure(nK = 250, nE=250, filename='bandstructure.png')

##########################################################
#				ESQC
##########################################################

# declare the ESQC solver and the energy range
trans=ESQCsolver(system)
trans.set_energy_range(emin=-8,emax=8,nE=500)

# compute the TE
te_esqc = trans.compute_transmission()

##########################################################
#				NEGF - WBL
##########################################################
trans=NEGFsolver(system)
trans.set_energy_range(emin=-8,emax=8,nE=500)


# set the wbl approx on and define the ldos of the elctrodes
trans.set_wide_band_limit(True)
trans.set_local_dos_electrode(0.25)

# compute the transport properties
te_negf_wbl = trans.compute_transmission()


##########################################################
#				NEGF - NO WBL
##########################################################
trans=NEGFsolver(system)
trans.set_energy_range(emin=-8,emax=8,nE=500)

# set the wbl approx OFF and define the 
# surface green function of the electrodes
trans.set_wide_band_limit(False)
trans.junction.precompute_electrodes_surface_green_function(np.linspace(-8.5,8.5,250),tol=1E-8,itermax=1E4,identical_electrodes=True)

# compute the surface green functions of the electrodes
te_negf = trans.compute_transmission()



##########################################################
#				Plot the results
##########################################################
plt.semilogy(trans.energies,te_esqc+1e-16,linewidth=4,color='black',label='ESQC')
plt.semilogy(trans.energies,te_negf_wbl+1e-16,linewidth=2,color='#DF1400',label='NEGF-WBL')
plt.semilogy(trans.energies,te_negf+1e-16,linewidth=2,color='#64D2FF',label='NEGF')
plt.ylim([1E-6,1E1])
plt.legend()
plt.savefig('te.png')






