# Huskyton

## Installation

Clone the repository and cd into it. Then use pip to install the module

```
pip install ./
```

The installation will compile the huckel dynamic library. For that you must have gcc installed. 


## Example: Tight-binding model of a benzene ring

Once the installation is done you can use the module as in the following example. In this example we compute the transmission of a 'benzene' molecule using simple electrodes containing only one states per unit cell. 

![junction](./pics/tb/junction.png)

Husky allows to compute the electronic transmission of this system (and any other systems defined similarly) using different theoretical framework. For example the ESQC and NEGF transmission calculated with Husky of the system above is represented below.

![te](./pics/tb/te.png)

As we can see below the ESQC results agree with NEGF calculation including the SGF of the elecrodes. However the WBL approximation leads to slightly different results as expeced. One can also clearly see the quantum interference around E=0 that are so typical of this system.

### Defining the junction geometry

Before doing any simulation on the transport we need to define the junction geomerty. We need to import the ```Junction``` class defined in ```husky.hamiltonian.ElectronicSystem```. Some basic definition of tight-binding hamiltonian are defined in ```husky.hamiltonian.model``` as for example cyclic armoatic molecules. We can add the hamiltonian of the central part with the ```add_central_matrix(h0)```

```python
import numpy as np
from husky.hamiltonian.ElectronicSystem import Junction
from husky.hamiltonian.model import tight_binding as tb

# define the junction
system = Junction()

# add the molecule hamiltonian
# here a cyclic aromatic with 6 centers, 
# onsite energy e =  0 eV 
# coupling      a = -1 eV
system.add_central_matrix(h0=tb.cyclic_aromatic(N=6,e=0,a=-1.))
```

We now define the electrodes. Here the electrodes contains only one state per unit cell. We can add as many electrodes as we want via the ```add_elecrode_matrix(h0,hi,vmol)```. Here ```h0``` is the hamiltonian of each unit cell, ```hi``` the interation between neighboring unit cells and ```vmol``` the interaction matrix between the last site of the electrode and the central part of the junction.

```python
# add single state electrodes
# onsite energy e0 = 0
# coupling between unit cell of the electrodes vi = -2.5
# coupling between the last site of the electorde and the molecule vmol 1
e0,vi,vmol = 0.,-2.5, 1.

# first electrode
system.add_electrode_matrix(h0 = np.array([e0]), hi = np.array([vi]),vmol = np.array([vmol,0,0,0,0,0]))

# second electrode
system.add_electrode_matrix(h0 = np.array([e0]), hi = np.array([vi]),vmol = np.array([0,0,vmol,0,0,0]))
```

The ```Elecrode``` class allow to easily compute the band structure of the electrodes via the ```band_structure(nK,nE,filename)``` method. Here ```nK``` an ```nE```  are the number of points to be used in K and energy space respectively. ```filename``` specify the .png file where the band structure is plotted.

```python
# compute the band structure of the electrode
system.electrodes['0'].band_structure(nK = 250, nE=250, filename='bandstructure_0.png')
```

![band structure](./pics/tb/bandstructure_0.png)

### Computing the transmission with ESQC

ESQC or Ellastic Scattering Quantum Chemistry is a powerfull method to compute the electronic transmission of molecular junction. We have used here the formalism defined in the references 

  * [1] Electronic transport calculations for self-assembled monolayers of 1,4-phenylene diisocyanide on Au111 contacts. Robert Dahlke and Ulrich Schollwock Phys Rev B 69 085324 2004

  * [2] Numerical calculations for electronic transport through molecular systems Robert Dahlke. PHD THESIS 2004

This method is implemented in the ```ESQCsolver``` defined in ```husky.transport.esqc```. To define an instance of the solver we must pass a Junction instance as defined above. The calculation of the transmission is then simply achieved with the ```compute_transmission()``` method.

```python
##########################################################
#				ESQC
##########################################################
from husky.transport.esqc import ESQCsolver

# declare the esqc solver and the enrgy range
trans=ESQCsolver(system)
trans.set_energy_range(emin=-5,emax=5,nE=500)

# compute the transmission and plot it
te_esqc = trans.compute_transmission()
```


### Computing the transmission with NEGF

Non-equlibrium Green functions are a particularly efficient method for quantum transport. We have here implemented this method in a ```NEGFsolver``` defined in ```husky.transport.negf```. As for the ```ESQCsolver``` a ```Junction``` instance must be passed during the definition of a ```NEGFsolver``` instance. 


#### Wide band limit approximation

The wide-band-limit approximation is a very popular simplification of the transport problem where the density of state of the electrodes is assumed constant. To use the WBL approximation we must set it true in the solver instance and give a value of the local DOS of the electrode. Once this is done the calculation of the transport is simply done by the ```compute_transmission()``` method.

```python
##########################################################
#				NEGF - WBL
##########################################################
from husky.transport.negf import NEGFsolver
trans=NEGFsolver(system)
trans.set_energy_range(emin=-5,emax=5,nE=500)


# set the wbl approx on and define the ldos of the elctrodes
trans.set_wide_band_limit(True)
trans.set_local_dos_electrode(0.4)

# compute the transport properties
te_negf_wbl = trans.compute_transmission()
```

#### Beyond the WBL approximation

However once can also compute the DOS of the electrode via their surface green function to obtain a more accurate values of the transmission. The calculation of the SGF is done via a method of the ```Junction``` instance. The SGF is precomputed before the calculation of the transport for efficiency reasons. Once the SGF calculated, the transmission is computed as previously via the ```compute_transmission()``` method.

```python
##########################################################
#				NEGF - NO WBL
##########################################################
trans=NEGFsolver(system)
trans.set_energy_range(emin=-5,emax=5,nE=500)

# set the wbl approx OFF and define the 
# surface green function of the electrodes
trans.set_wide_band_limit(False)
trans.junction.precompute_electrodes_surface_green_function(np.linspace(-8,8,250),tol=1E-6,identical_electrodes=True)

# compute the surface green functions of the electrodes
te_negf = trans.compute_transmission()
```

![sgf](./pics/tb/sgf.png)

### Plotting the results

Once everything is calculated we can plot the result to compare the different methods.

```python
##########################################################
#				Plot the results
##########################################################
import matplotlib.pyplot as plt 

plt.semilogy(trans.energies,te_esqc,linewidth=4,color='black',label='ESQC')
plt.semilogy(trans.energies,te_negf_wbl,linewidth=2,color='#DF1400',label='NEGF-WBL')
plt.semilogy(trans.energies,te_negf,linewidth=2,color='#64D2FF',label='NEGF')
plt.ylim([1E-6,2])
plt.legend()
plt.savefig('te.png')

```



