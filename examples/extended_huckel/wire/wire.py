from husky.hamiltonian.ElectronicSystem import Junction
from husky.transport.negf import NEGFsolver
from husky.transport.esqc import ESQCsolver


import matplotlib.pyplot as plt 
import numpy as np
import scipy.linalg as spla

import sys


############################################################################
# compute the transmission coefficient of a monoatomic gold wire
############################################################################

# create the junction 
system = Junction()

# add the xyz files of the different part
system.add_molecule_xyz('./mol.xyz')
system.add_electrode_xyz('./elec1.xyz')
system.add_electrode_xyz('./elec2.xyz')

# add the exended Huckel parameters
system.add_huckel_parameters('./Au.param')

# compute the hamiltonians/overlaps using extended huckel theory
# the lowdin_ortho option controls the orthogonalization of the matrices
# if true a lowdin orthogonalization is performed.
system.compute_matrices_huckel(lowdin_ortho=False)

# compute/print the band structure of the electrode
system.electrodes['0'].band_structure(nK = 500, nE=150)

##########################################################
#				ESQC
##########################################################

trans = ESQCsolver(system)
trans.set_energy_range(emin=-17,emax=0,nE=500)
# compute the TE
te_esqc = trans.compute_transmission()

##########################################################
#				NEGF - WBL
##########################################################
trans=NEGFsolver(system)
trans.set_energy_range(emin=-17,emax=0,nE=500)


# set the wbl approx on and define the ldos of the elctrodes
trans.set_wide_band_limit(True)
trans.set_local_dos_electrode(0.5)

# compute the transport properties
te_negf_wbl = trans.compute_transmission()


##########################################################
#				NEGF - NO WBL
##########################################################
trans=NEGFsolver(system)
trans.set_energy_range(emin=-17,emax=0,nE=500)

# set the wbl approx OFF and define the 
# surface green function of the electrodes
trans.set_wide_band_limit(False)
trans.junction.precompute_electrodes_surface_green_function(np.linspace(-18,1,250),tol=1E-6,identical_electrodes=True)

# compute the surface green functions of the electrodes
te_negf = trans.compute_transmission()



##########################################################
#				Plot the results
##########################################################
plt.semilogy(trans.energies,te_esqc+1e-16,color='grey',linewidth=4,label='ESQC')
plt.semilogy(trans.energies,te_negf_wbl+1e-16,color='blue',label='NEGF-WBL')
plt.semilogy(trans.energies,te_negf+1e-16,color='red',label='NEGF')
plt.ylim([1E-6,10])
plt.legend()
plt.savefig('te.png')








