#ifndef HEADER_H
#define HEADER_H

  #define SQRT3        1.73205080756888
  #define SQRT6        2.44948974278318
  #define SQRT10       3.16227766016838
  #define SQRT15       3.87298334620742
  #define AUI          1.889644746

  #define TABLESIZE    103
  #define MAXCOL       80
  #define MAXATOM      100
  #define MAXSIZE      10000

  typedef struct
		  {
		  int orb_of_e;
		  double VSIP;
		  double exp;
		  double exp2[2];
		  double coef2[2];
		  } single_orb;
			  


  typedef struct
		  {
		  char symbol[3];
		  int  valence_electron; 
		  single_orb orb[4];
		  } atom_parameter;

		  
  typedef struct
		  {
		  char atomTypeChar[5];
		  int atomtype;
		  double x;
		  double y;
		  double z;
		  } atom;

	
  typedef struct{
    int atm_index_table;
    char atm[5];
    int atmNbr;
    char type[3];
  } ORB;

#endif
