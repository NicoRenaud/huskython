import hkl_module as hkl
import numpy as np
import ctypes

#dir(hkl)
filename = './benzene.in'

nb_orb = hkl.nbOrb(filename)

h = np.zeros(nb_orb*nb_orb)
s = np.zeros(nb_orb*nb_orb)

hkl.hklHam(h,s,nb_orb,filename)
h = h.reshape(nb_orb,nb_orb)
s = s.reshape(nb_orb,nb_orb)
print h.shape, s.shape
