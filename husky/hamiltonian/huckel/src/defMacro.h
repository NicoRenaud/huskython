#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>


#ifndef _defMacro_H

#define _defMacro_H
#define KB 8.67 * pow(10,-5) 
#define CAHRGE 1.602176468 * pow(10,-19)
#define HBAR 1.05 * pow(10,-34)
#define FS2UA 1.5192678
#define PI 3.14116457148 
#define DBG 1
#define AFF_SCREEN 1
#define NEGF_VERBOSE 0
#define DIST_FIRST_NEIGHBOR 1.5
#define MAX_CLUSTER 10
#define VERBOSE 1


// distance between S and Au
#define distSAu 2.25

// coupling between between S and Au
#define CPLG_SAU   -5.5	// s orbitals
#define CPLG_SAU_s -14.3
#define CPLG_SAU_p -6.3	// p orbitals

// coupling between Au and Au
#define CPLG_AUAU -6  // s orbitals
#define CPLG_AUAU_s -10.7  // s orbitals
#define CPLG_AUAU_p -4.25  // p orbitals
#define CPLG_AUAU_d -2.5  // d orbitals

///////////////////////////////////////
// 	 STRUCTURES : SYSTEM
///////////////////////////////////////


typedef struct{
  

  // position of the atoms
  char pos[100];
  
  // metho to compute hamiltonian
  char elec_struct[100];
  
  // use or not the overlap
  char use_overlap[10];
  
  // number MO to output
  int nbr_mo_1;
  int nbr_mo_2;
  int nbr_mo_3;
  int *index_mo_1;
  int *index_mo_2;
  int *index_mo_3;
  
  // nb atome
  int index_contact[2]; 
  double cplg_contact[2]; // obsolete ?
  char orb_contact[100];
  
  // size of the final cluster
  char cluster[100];
  char orb_elec[100];
  
  
  // output options
  char export_MO[5];
   
} SYS_XYZ;



  



///////////////////////////////////////
//
//    Proto des fonctions Lapack
//
///////////////////////////////////////

extern void ssyev_( char *jobz, char *uplo, int *n, float *a, int *lda,
        float *w, float *work, int *lwork, int *info );

extern int dsyev_(char *jobz, char *uplo, int *n, double *a,  int *lda, 
	double *w, double *work, int *lwork, int *info);

extern void cgeev_( char* jobvl, char* jobvr, int* n, complex float* a,
                int* lda, complex float* w, complex float* vl, int* ldvl, complex float* vr, int* ldvr,
                complex float* work, int* lwork, float* rwork, int* info );

extern int dggev_(char *jobvl, char *jobvr, int *n, double *
                  a, int *lda, double *b, int *ldb, double *alphar, 
                  double *alphai, double *beta, double *vl, int *ldvl,
                  double *vr, int *ldvr, double *work, int *lwork,
                  int *info);

extern void zgeev_( char* jobvl, char* jobvr, int* n, complex double* a,
                   int* lda, complex double* w, complex double* vl, int* ldvl, complex double* vr, int* ldvr,
                   complex double* work, int* lwork, double* rwork, int* info );

extern int cgetrf_(int *m, int *n, complex *a, int *lda,
                   int *ipiv, int *info);

extern int cgetri_(int *n, complex *a, int *lda, int *
                   ipiv, complex *work, int *lwork, int *info);

extern int zgetrf_(int *m, int *n, complex double *a,
                   int *lda, int *ipiv, int *info);

extern int zgetri_(int *n,complex double *a, int *lda,
                   int *ipiv, complex double *work, int *lwork, int *info);

int dsyev_(char *jobz, char *uplo, int *n, double *a,
           int *lda, double *w, double *work, int *lwork,
           int *info);


int dgeev_(char *jobvl, char *jobvr, int *n, double *
           a, int *lda, double *wr, double *wi, double *vl,
           int *ldvl, double *vr, int *ldvr, double *work,
           int *lwork, int *info);

#endif