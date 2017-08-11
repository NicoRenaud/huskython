import numpy as np
import math

###########################################################
#
# block diagonalization
#
###########################################################
def compute_bd_transformation(H,index1,index2):

	l,S = np.linalg.eigh(H)
	S11 = S[np.ix_(index1,index1)]
	S12 = S[np.ix_(index1,index2)]
	S21 = S[np.ix_(index2,index1)]
	S22 = S[np.ix_(index2,index2)]

	if len(index2)>1:
		iS22 = np.linalg.inv(S22)
	else:
		iS22 = 1./S22

	X = np.dot(S12,iS22)
	U = np.eye(len(H))

	U[np.ix_(index1,index2)] = X
	U[np.ix_(index2,index1)] = -X.T

	F = np.linalg.inv(spla.sqrtm(np.dot(U.T,U),-0.5))
	T = np.dot(U,F)

	return T

###########################################################
#
# compute the transport
#
###########################################################
def lowdin(H,S):

	# Diagonalize the overlap matrix
	n = S.shape[0]
	s, U = np.linalg.eigh(S)

	#check for lienar dependence in the Fock matrix
	check = sum(s<1E-8)
	if check > 0:
		print(' \t \t ====== Warning ====== \n \t \t == Linear dependence found in the overlap matrix')
		ind = where(s<1E-8)
		s[ind] = 1E16
		U[:,ind] = 0
		U[ind,ind] = 1
		H[:,ind] = 0
		H[ind,:] = 0
		print(' \t \t == %d basis function discarded in the Fock matrix \n \t \t ====== ======= ====== ' % len(ind))


	V = np.dot(np.dot(U, np.diag(1 / np.sqrt(s))), U.T)

	# Transform the Hamiltonian
	H = np.dot(np.dot(V, H), V.T)

	return H

###########################################################
#
# define the rotation matrix
#
###########################################################
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])