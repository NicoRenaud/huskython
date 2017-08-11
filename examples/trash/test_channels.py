import numpy as np


def missing_elements(index,ind_max):
	miss = []
	for i in range(ind_max):
		if i not in index:
			miss.append(i)
	return miss

def sort_inout_channels(L):

	eps = 1E-12
	l = np.copy(L)
	n = len(l)
	index = []

	# find the pairs that respect lp=1./lm^*
	for i in range(n):

		if l[i] == None:
			continue

		for j in range(n):

			if i == j or l[j] == None :
				continue

			if np.abs(l[i] - 1./np.conj(l[j]))<eps:
				index.append([i,j])
				l[i] = None
				l[j] = None

	# fint the potential missing dudes
	miss = []
	for i in range(len(L)):
		if i not in np.array(index).flatten().tolist():
			miss.append(i)

	# create the total list
	if len(index) == 0:
		index =  miss
	else:
		index = np.array(index)[:,0].tolist()+ miss[:len(miss)/2] + np.array(index)[:,1].tolist() + miss[len(miss)/2:]

	print index
	print miss
	return index

r = np.random.rand(5)
r = np.append(r,r[-1])

#R = np.concatenate([np.exp(-1j*r),1./np.exp(1j*r)])
R = np.array([])
R = np.append(R,1.3)
R = np.append(R,-2)
R = np.append(R,8)
R = np.append(R,-12)
np.random.shuffle(R)

print R
pairs=sort_inout_channels(R)


