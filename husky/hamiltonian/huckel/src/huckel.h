
#include "./header.h"

#ifndef _comp_HKL_H_h
#define _comp_HKL_H_h


int findIndex_reduced(char *name, char *PARAM);
void read_molecule_input_file(char *file_name, atom **MOL, int *no_atoms, double *KEHT, char *PARAM );
void extract_atom_mol(atom *mol_only,  int *ind_atom_mol, int no_atom_mol, atom *MOL);
void overlap(int atom_row,int atom_col,double delx,double dely,double delz,double S[16][16],double H[16][16],double KEHT);
void mov(double *sigma,double *pi,double *delta,double *phi,int atom_col,int atom_row,double rr,int n1,int n2,int l1,int l2);
void abfns(double *a,double *b,double sk1,double sk2,double rr,int maxcal);
double lovlap(double *a,double *b,double sk1,double sk2,double r,int l1,int l2,int m1,int n1,int n2);
void read_atomic_parameters(char * HUCKEL_PARAM);
int compute_nb_orb(char *file_name);
void compute_huckel_hamiltonian_general(double *hmat, double *smat, int nb_orb_precomp, char *file_name);


/*
OLD DEFINITION JUST FOR LEGACY AND CHECK
void compute_huckel_hamiltonian(double *Hout, double *Sout,
				int sizeH, int no_atoms, atom *molecule,char *orb_elec,char *type_struct_elec, char *use_overlap);
void compute_huckel_hamiltonian_general(int *sizeH, int *size_sys, atom **MOL,char *file_name,char *out_path);
void read_atomic_parameters(char *PARAM);
void overlap(int atom_row,int atom_col,double delx,double dely,double delz,double S[16][16],double H[16][16],double KEHT);
void mov(double *sigma,double *pi,double *delta,double *phi,int atom_col,int atom_row,double rr,int n1,int n2,int l1,int l2);
void abfns(double *a,double *b,double sk1,double sk2,double rr,int maxcal);
double lovlap(double *a,double *b,double sk1,double sk2,double r,int l1,int l2,int m1,int n1,int n2);
// not sure if these ones are still used anywhere
void defInteraction(double *V1, double *S1, double *V2, double *S2, 
		    int *index_orb_mol, int *index_orb_elec_1, int *index_orb_elec_2, int nb_orb_elec_1, int nb_orb_elec_2, int nb_orb_mol,
		    double *H, double *S, int nb_orb_tot				);
int compute_nb_orb(char *file_name);
void compute_orb_index(int *index_orb_mol, int *index_orb_elec_1, int *index_orb_elec_2, int nb_orb_mol, int nb_orb_elec_1, int nb_orb_elec_2);
void orbitale_information(ORB *orbitale, atom *molecule, int nb_atm, int nb_orb );
*/
#endif