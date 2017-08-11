import numpy as np

###########################################################
#
# find the index that corresponds to a given atom name
#
###########################################################
def find_index(name):
	#list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
	# replaced Fe by Pb otherwise Jmol can't print MO on Pb ....
	list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Pb", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Fe", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
	index= [i for i, x in enumerate(list) if x == name]
	return index[0]+1

###########################################################
#
# find the slater expoants of a given atom
#
###########################################################
def get_slater_coeff(name,param):

	f =open(param,'r')
	data = f.readlines()
	f.close

	for i in range(len(data)):
		data[i] = data[i].split()
		if data[i][0] == name:
			sc1 = data[i][7]
			sc2 = data[i][13]
			sc3 = data[i][19]
			break
	return float(sc1),float(sc2),float(sc3)


###########################################################
#
# print the MO in mopac format
#
###########################################################	
def print_mo_mopac(wH,uH,s,system,filename,PARAM):

	natom = len(system)
	natom_print = natom

	nb_orb = len(uH)
	count = 0

	# comute the inverse of S
	invS = np.linalg.inv(s)


	# open the output file
	f = open(filename,'w')

	# header
	f.write("        %d MOPAC-Graphical data Version 2012.13.084W\n" %natom_print)

	# print the atoms
	for i in range(natom):
		at = find_index(system[i].symbol)
		f.write('%4d %*f%*f%*f  0.0000\n' %(at,12,system.positions[i,0],12,system.positions[i,1],12,system.positions[i,2]))

	# print the slater exponents
	for i in range(natom):
		sc1,sc2,sc3 = get_slater_coeff(system[i].symbol,PARAM)
		f.write("  %1.7f  %1.7f  %1.7f\n" %(sc1,sc2,sc3))

	occ = np.zeros(nb_orb)

	# print the orbitals
	for iorb in range(nb_orb):

		
		f.write(" ORBITAL %3d  A  %2.5g\n" %(occ[iorb],wH[iorb]))

		for jcomp in range(nb_orb):
			f.write("% 1.8E" %(uH[jcomp,iorb]).real)
			count += 1
			if count == 5:
				f.write("\n")
				count  = 0
		if count>0:
				f.write("\n")
				count  = 0

	# print the inverse matrix
	count = 0
	f.write("INVERSE_MATRIX[%dx%d]=\n"%(nb_orb,nb_orb))
	for i in range(nb_orb):
		for j in range(i+1):
			f.write("% 1.8E" %(invS[j,i]))
			count+=1
			if count == 5:
				f.write("\n")
				count  = 0
		if count>0:
				f.write("\n")
				count  = 0
	f.write(" Keywords: SYMMETRY GRAPHF")
	f.close()

