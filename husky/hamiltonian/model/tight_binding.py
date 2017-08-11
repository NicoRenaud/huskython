
import numpy as np 

def cyclic_aromatic(N=6,e=0,a=1):
	h = e*np.eye(N)
	for i in range(N-1):
		h[i,i+1] = a 
		h[i+1,i] = a 
	h[0,-1] = a
	h[-1,0] = a
	return h

def linear_aromatic(N=6,e=0,a=1):
	h = e*np.eye(N)
	for i in range(N-1):
		h[i,i+1] = a 
		h[i+1,i] = a 
	return h