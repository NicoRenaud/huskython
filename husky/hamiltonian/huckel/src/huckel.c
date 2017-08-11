// This is the new code
// This was part of Yaehmop and the doc has been lost in history
// would be nice to have a pure python implementation not to have 
// to create a dynamic library and call it in python. 

#include "./defMacro.h"
#include "./header.h"


// global variable I switched the static on
static  atom_parameter table[TABLESIZE];
static  atom_parameter *period=table-1;


//static  double hmat[MAXSIZE][MAXSIZE],smat[MAXSIZE][MAXSIZE];
//i dont want that sort of allocation
//why are they global they only appear in 1 function !!
//static  double *hmat, *smat;
static  int no_atoms;
static  int no_orbitals;




////////////////////////////////////////
//	find index from the reduced PARAM LIST
////////////////////////////////////////

int findIndex_reduced(char *name, char *PARAM)
{
  
  int index = 0;
	FILE *Fparameter;
	int i,j;
	char atom_name[3];

	// open the parameter file and read it
	Fparameter=fopen(PARAM,"r");
  
	if(Fparameter == NULL)
	{
	  //printf("Huckel parameters not found at:  /home/nico/Bureau/husky.1/SRC/huckel/parameters \n");
    printf("Huckel parameters not found at:  %s \n",PARAM);
	  exit(1);
	}


	for(i=1;i<=250;i++)
	{
		
		// read the atom name
		fscanf(Fparameter,"%s",atom_name);
		

		
		if (strcmp(name,atom_name) == 0)
		{
			index = i;
			break;
		}
		
		// read the rest
		fscanf(Fparameter,"%*i");
	
		for(j=0;j<4;j++)
		{
			fscanf(Fparameter,"%*i ");
		}
	
		for(j=0;j<4;j++)
			fscanf(Fparameter,"%*f %*f %*f %*f %*f %*f");
		fscanf(Fparameter,"\n");

  }
	fclose(Fparameter);

	// return the index 
  return(index);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
//			Read the input file
////////////////////////////////////////////////////////////////////////////////////////////////////

void read_molecule_input_file(char *file_name, atom **MOL, int *no_atoms, double *KEHT, char *PARAM )
{
  
  int i;
  int INT;
  int INT2;
  float X,Y,Z;
  
  int current_cell = 0;
  int read_cell[2];
  
  int no_atom_tot=-1;
  int nb_atom_tot_read=0;
  
		
  int *no_atom_cell;
  no_atom_cell = calloc(2,sizeof(int));
  no_atom_cell[0]  = -1;
  no_atom_cell[1]  = -1;
  
  int *no_atom_cell_read;
  no_atom_cell_read = calloc(2,sizeof(int));
  no_atom_cell_read[0]  = 0;
  no_atom_cell_read[1]  = 0;
  

  
  
  atom *MOL_TEMP=NULL;
  atom *MORE_MOL;
  
  char str_cmt[5];
  char str[100];
	

  int *index_cell_1=NULL;
  int *index_cell_2=NULL;
  
  int *more_index_cell_1;
  int *more_index_cell_2;
  
  int index_check;
  
  FILE *f;
  
  
  // check if the file exists
   f = fopen(file_name,"r");
   if (!f) {
      printf("couldn't read input file name %s\n", file_name);
      exit(1);
   }
   
   //////////////////////////////////////////
   // for all lines of the input file
   //////////////////////////////////////////
   while (!feof(f))
   //while(fgets(str,20,f) != NULL)
   {
     
			///////////////////////
			// read the string
			///////////////////////
			fscanf(f,"%s",str);
			

			
			////////////////////////
			// if it's a comment
			////////////////////////
			strncpy(str_cmt,str,2);
			str_cmt[2] = '\0';
			if (!strcmp(str_cmt,"//"))
				fgets(str,100,f);
			
			///////////////////////////
			// if it's not a comment
			///////////////////////////
			else
			{
				
		 ////////////////////////////////////////////////////////////////////////////////
		 //				GENERAL INFORMATIONS
		 ////////////////////////////////////////////////////////////////////////////////
		 
					////////////////////////
					// total number of atoms
					///////////////////////    
					if( !strcmp(str,"nb_atom_tot") || !strcmp(str,"NB_ATOM_TOT")   || 
					!strcmp(str,"nb_atom")  || !strcmp(str,"NB_ATOM")  ) 
					{
						fscanf(f,"%d",&no_atom_tot);		  
						MOL_TEMP = calloc(no_atom_tot,sizeof(atom));
					}  
				
					////////////////////////
					// mixing term for ehmo
					///////////////////////
					else if( !strcmp(str,"Keht") || !strcmp(str,"KEHT")   ||
					!strcmp(str,"keht")  || !strcmp(str,"kEHT")  ) 
					{
						fscanf(f,"%lf",KEHT);
					}
					
					////////////////////////
					// EHMO parameters
					///////////////////////    
					else if( !strcmp(str,"parameters") || !strcmp(str,"PARAMETERS")   ||
					!strcmp(str,"param")  || !strcmp(str,"PARAM")  ) 
					{
						fscanf(f,"%s",PARAM);
					}
				
		 ////////////////////////////////////////////////////////////////////////////////
		 //			READING POSITIONS
		 ///////////////////////////////////////////////////////////////////////////////
								
					else
					{
						 
						 index_check = findIndex_reduced(str,PARAM);
						 //printf("str : %s index : %d\n",str,index_check);
						 if (index_check > 0)
						 {
				
							// read the data
							fscanf(f,"%f%f%f",&X,&Y,&Z);
							
							// realloc the memory if not dine previously
							if(no_atom_tot == -1 || nb_atom_tot_read >= no_atom_tot)
							{
								MORE_MOL = realloc(MOL_TEMP,(nb_atom_tot_read+1)*sizeof(atom));
								MOL_TEMP = MORE_MOL;
							}
							
							// store the data
							strcpy(MOL_TEMP[nb_atom_tot_read].atomTypeChar,str);
							MOL_TEMP[nb_atom_tot_read].atomtype = index_check;
							MOL_TEMP[nb_atom_tot_read].x=X;
							MOL_TEMP[nb_atom_tot_read].y=Y;
							MOL_TEMP[nb_atom_tot_read].z=Z;
								
							// one more atom read!
							nb_atom_tot_read ++;
				
						}
					}
				}
			}

      // export the data
      *MOL = MOL_TEMP;
      *no_atoms = nb_atom_tot_read;
	
      if(nb_atom_tot_read != no_atom_tot)
					printf("\nWarning: input claims there is %d atom in the system but %d have been found \n Please check %s\n",no_atom_tot,nb_atom_tot_read,file_name);
	
       
      // close the file
      fclose(f);

}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void extract_atom_mol(atom *mol_only,  int *ind_atom_mol, int no_atom_mol, atom *MOL)
{
  
  int i;
  for(i=0;i<no_atom_mol;i++)
  {
   	strcpy(mol_only[i].atomTypeChar,MOL[ind_atom_mol[i]].atomTypeChar);
	mol_only[i].atomtype = MOL[ind_atom_mol[i]].atomtype;
	mol_only[i].x = MOL[ind_atom_mol[i]].x;
	mol_only[i].y = MOL[ind_atom_mol[i]].y;
	mol_only[i].z = MOL[ind_atom_mol[i]].z; 
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////
//
//  COMPUTE THE OVERLAP BETWEEN ORBITALS
//
//////////////////////////////////////////////////////////////////////////////////////////////
void overlap(int atom_row,int atom_col,double delx,double dely,double delz,double S[16][16],double H[16][16],double KEHT)
{
	double rt2,r,t,ca,cb,sa,sb,ca2,sa2,cb2,sb2,cbsb,casa,cb2sb2,s2b,sa3,c2b,s3b,c3b,s2a,c2a,VSIP_row,VSIP_col;
	double pt[9],dt[25],ft[49],*ptr,*dtr,*ftr;
	double sigma,pi,delta,phi;
	int j,k,nsrow,nprow,ndrow,nfrow,nscol,npcol,ndcol,nfcol;
	
	void mov(double *,double *,double *,double *,int,int,double,int,int,int,int);

	
	ptr=pt-1;
	dtr=dt-1;
	ftr=ft-1;


	rt2=delx*delx+dely*dely;
	r=sqrt(rt2+delz*delz);


	if(rt2 <= 1e-10)
	{
		cb=1.0;
		sb=0.0;
		sa=0.0;
	}
	else
	{
		t=sqrt(rt2);
		cb=delx/t;
		sb=dely/t;
		sa=t/r;
	}
	
	ca=delz/r;

	ptr[1]=sa*cb;
	ptr[2]=sa*sb;
	ptr[3]=ca;
	ptr[4]=ca*cb;
	ptr[5]=ca*sb;
	ptr[6]=-sa;
	ptr[7]=-sb;
	ptr[8]=cb;
	ptr[9]=0.0;

	if(period[atom_row].orb[2].orb_of_e+period[atom_col].orb[2].orb_of_e>0) 
	{
		ca2=ca*ca;
		sa2=sa*sa;
		cb2=cb*cb;
		sb2=sb*sb;
		cbsb=cb*sb;
		casa=ca*sa;
		cb2sb2=cb2-sb2;
		dtr[1]=SQRT3*0.5*sa2*cb2sb2;
		dtr[2]=1.0-1.5*sa2;
		dtr[3]=SQRT3*cbsb*sa2;
		dtr[4]=SQRT3*casa*cb;
		dtr[5]=SQRT3*casa*sb;
		dtr[6]=casa*cb2sb2;
		dtr[7]=-SQRT3*casa;
		dtr[8]=2.0*casa*cbsb;
		dtr[9]=cb*(ca2-sa2);
		dtr[10]=sb*(ca2-sa2);
		dtr[11]=-2.0*sa*cbsb;
		dtr[12]=0.0;
		dtr[13]=sa* cb2sb2;
		dtr[14]=-ptr[5];
		dtr[15]=ptr[4];

		if(period[atom_row].orb[2].orb_of_e*period[atom_col].orb[2].orb_of_e>0) 
		{
			dtr[16]=0.5*(1.0+ca2)*cb2sb2;
			dtr[17]=0.5*SQRT3*sa2;
			dtr[18]=cbsb*(1.0+ca2);
			dtr[19]=-casa*cb;
			dtr[20]=-casa*sb;
			dtr[21]=-2.0*ca*cbsb;
			dtr[22]=0.0;
			dtr[23]=ca*cb2sb2;
			dtr[24]=ptr[2];
			dtr[25]=-ptr[1];
		}
	}


	if(period[atom_row].orb[3].orb_of_e+period[atom_col].orb[3].orb_of_e>0) 
	{
		s2b=2.0*sb*cb;
		sa3=sa2*sa;
		c2b=(cb2-sb2);
		s3b=(c2b*sb + s2b*cb);
		c3b=(c2b*cb - s2b*sb);
		s2a=2.0*sa*ca;
		ftr[1]=0.5*ca*(5.0*ca2 -3.0);
		ftr[2]=SQRT6*0.25*cb*sa*(5.0*ca2 -1.0);
		ftr[3]=SQRT6*0.25*sb*sa*(5.0*ca2 -1.0);
		ftr[4]=SQRT15*0.5*s2b*ca*sa2;
		ftr[5]=SQRT15*0.5*c2b*ca*sa2;
		ftr[6]=SQRT10*0.25*c3b*sa3;
		ftr[7]=SQRT10*0.25*s3b*sa3;
		ftr[8]=-SQRT6*0.25*sa*(5.0*ca2 -1.0);
		ftr[9]=0.25*cb*ca*(15.0*ca2 -11.0);
		ftr[10]=0.25*sb*ca*(15.0*ca2 -11.0);
		ftr[11]=SQRT10*0.25*s2b*sa*(3.0*ca2 - 1.0);
		ftr[12]=SQRT10*0.25*c2b*sa*(3.0*ca2 - 1.0);
		ftr[13]=SQRT15*0.25*c3b*ca*sa2;
		ftr[14]=SQRT15*0.25*s3b*ca*sa2;
		ftr[15]=0.0;
		ftr[16]=-0.25*sb*(5.0*ca2 -1.0);
		ftr[17]=0.25*cb*(5.0*ca2 -1.0);
		ftr[18]=SQRT10*0.25*c2b*s2a;
		ftr[19]=-SQRT10*s2b*s2a*0.25;
		ftr[20]=-SQRT15*0.25*s3b*sa2;
		ftr[21]=SQRT15*0.25*c3b*sa2;

		if(period[atom_row].orb[2].orb_of_e*period[atom_col].orb[2].orb_of_e>0) 
		{
			c2a=ca2-sa2;
			ftr[22]=0.0;
			ftr[23]=SQRT10*0.5*sb*ca*sa;
			ftr[24]=-SQRT10*0.5*cb*ca*sa;
			ftr[25]=c2b*c2a;
			ftr[26]=-s2b*c2a;
			ftr[27]=-SQRT6*0.25*s3b*s2a;
			ftr[28]=SQRT6*0.25*c3b*s2a;
			ftr[29]=SQRT15*0.5*ca*sa2;
			ftr[30]=SQRT10*0.25*cb*sa*(1.0 -3.0*ca2);
			ftr[31]=SQRT10*0.25*sb*sa*(1.0 -3.0*ca2);
			ftr[32]=0.5*s2b*ca*(3.0*ca2 -1.0);
			ftr[33]=0.5*c2b*ca*(3.0*ca2 -1.0);
			ftr[34]=SQRT6*0.25*c3b*sa*(1.0 + ca2);
			ftr[35]=SQRT6*0.25*s3b*sa*(1.0 + ca2);
		}

		if(period[atom_row].orb[3].orb_of_e+period[atom_col].orb[3].orb_of_e>0) 
		{
			ftr[36]=-SQRT10*0.25*sa3;
			ftr[37]=SQRT15*0.25*cb*ca*sa2;
			ftr[38]=SQRT15*0.25*sb*ca*sa2;
			ftr[39]=-SQRT6*0.25*s2b*sa*(1.0 + ca2);
			ftr[40]=-SQRT6*0.25*c2b*sa*(1.0 + ca2);
			ftr[41]=0.25*c3b*ca*(3.0 + ca2);
			ftr[42]=0.25*s3b*ca*(3.0 + ca2);
			ftr[43]=0.0;
			ftr[44]=-SQRT15*0.25*sb*sa2;
			ftr[45]=SQRT15*0.25*cb*sa2;
			ftr[46]=-SQRT6*0.25*c2b*s2a;
			ftr[47]=SQRT6*0.25*s2b*s2a;
			ftr[48]=-0.25*s3b*(1.0 +3.0*ca2);
			ftr[49]=0.25*c3b*(1.0 + 3.0*ca2);
		}
	}


	r=AUI*r;


	nsrow=period[atom_row].orb[0].orb_of_e;
	nprow=period[atom_row].orb[1].orb_of_e;
	ndrow=period[atom_row].orb[2].orb_of_e;
	nfrow=period[atom_row].orb[3].orb_of_e;
	nscol=period[atom_col].orb[0].orb_of_e;
	npcol=period[atom_col].orb[1].orb_of_e;
	ndcol=period[atom_col].orb[2].orb_of_e;
	nfcol=period[atom_col].orb[3].orb_of_e;


/*	(s_row:s_col)	*/

	if(nsrow*nscol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nscol,nsrow,0,0);
		S[0][0]=sigma;
	}


/*	(s_row:p_col)	*/

	if(nsrow*npcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,npcol,nsrow,1,0);
		sigma=-sigma;

		for(k=1;k<=3;k++)
			S[0][k]=ptr[k]*sigma;
	}


/*	(s_row:d_col)	*/

	if(nsrow*ndcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,ndcol,nsrow,2,0);

		for(k=1;k<=5;k++)
			S[0][3+k]=dtr[k]*sigma;
	}


/*	(s_row:f_col)	*/

	if(nsrow*nfcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nfcol,nsrow,3,0);
		sigma=-sigma;

		for(k=1;k<=7;k++)
			S[0][8+k]=ftr[k]*sigma;
	}



/*	(p_row:s_col)	*/

	if(nprow*nscol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nscol,nprow,0,1);

		for(j=1;j<=3;j++)
			S[j][0]=ptr[j]*sigma;
	}



/*	(p_row:p_col)	*/

	if(nprow*npcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,npcol,nprow,1,1);
		sigma=-sigma;

		for(j=1;j<=3;j++)
		{
			for(k=1;k<=3;k++)
			{
				S[j][k]=ptr[j]*ptr[k]*sigma+(ptr[j+3]*ptr[k+3]+ptr[j+6]*ptr[k+6])*pi;
				S[k][j]=S[j][k];
			}
		}
	}



/*	(p_row:d_col)	*/

	if(nprow*ndcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,ndcol,nprow,2,1);
		pi=-pi;

		for(j=1;j<=3;j++)
		{
			for(k=1;k<=5;k++)
			{
				S[j][3+k]=ptr[j]*dtr[k]*sigma+(ptr[j+3]*dtr[k+5]+ptr[j+6]*dtr[k+10])*pi;
			}
		}
	}



/*	(p_row:f_col)	*/

	if(nprow*nfcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nfcol,nprow,3,1);
		sigma=-sigma;

		for(j=1;j<=3;j++)
		{
			for(k=1;k<=7;k++)
			{
				S[j][8+k]=ptr[j]*ftr[k]*sigma+(ptr[j+3]*ftr[k+7]+ptr[j+6]*ftr[k+14])*pi;
			}
		}
	}



/*	(d_row:s_col)	*/

	if(ndrow*nscol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nscol,ndrow,0,2);
		
		for(j=1;j<=5;j++)
			S[3+j][0]=dtr[j]*sigma;
	}



/*	(d_row:p_col)	*/

	if(ndrow*npcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,npcol,ndrow,1,2);
		sigma=-sigma;

		for(j=1;j<=5;j++)
		{
			for(k=1;k<=3;k++)
			{
				S[3+j][k]=dtr[j]*ptr[k]*sigma+(dtr[j+5]*ptr[k+3]+dtr[j+10]*ptr[k+6])*pi;
			}
		}
	}



/*	(d_row:d_col)	*/

	if(ndrow*ndcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,ndcol,ndrow,2,2);
		pi=-pi;

		for(j=1;j<=5;j++)
		{
			for(k=1;k<=5;k++)
			{
				S[3+j][3+k]=dtr[j]*dtr[k]*sigma+(dtr[j+5]*dtr[k+5]+dtr[j+10]*dtr[k+10])*pi+(dtr[j+15]*dtr[k+15]+dtr[j+20]*dtr[k+20])*delta;
				S[3+k][3+j]=S[3+j][3+k];
			}
		}
	}



/*	(d_row:f_col)	*/

	if(ndrow*nfcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nfcol,ndrow,3,2);
		sigma=-sigma;
		delta=-delta;

		for(j=1;j<=5;j++)
		{
			for(k=1;k<=7;k++)
			{
				S[3+j][8+k]=dtr[j]*ftr[k]*sigma+(dtr[j+5]*ftr[k+7]+dtr[j+10]*ftr[k+14])*pi+(dtr[j+15]*ftr[k+21]+dtr[j+20]*ftr[k+28])*delta;
			}
		}
	}



/*	(f_row:s_col)	*/

	if(nfrow*nscol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nscol,nfrow,0,3);

		for(j=1;j<=7;j++)
			S[8+j][0]=ftr[j]*sigma;
	}



/*	(f_row:p_col)	*/

	if(nfrow*npcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,npcol,nfrow,1,3);
		sigma=-sigma;

		for(j=1;j<=7;j++)
		{
			for(k=1;k<=3;k++)
			{
				S[8+j][k]=ftr[j]*ptr[k]*sigma+(ftr[j+7]*ptr[k+3]+ftr[j+14]*ptr[k+6])*pi;
			}
		}
	}



/*	(f_row:d_col)	*/

	if(nfrow*ndcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,ndcol,nfrow,2,3);
		pi=-pi;

		for(j=1;j<=7;j++)
		{
			for(k=1;k<=5;k++)
			{
				S[8+j][3+k]=ftr[j]*dtr[k]*sigma+(dtr[k+5]*ftr[j+7]+dtr[k+10]*ftr[j+14])*pi+(ftr[j+21]*dtr[k+15]+ftr[j+28]*dtr[k+20])*delta;
			}
		}
	}



/*	(f_row:f_col)	*/

	if(nfrow*nfcol>0)
	{
		mov(&sigma,&pi,&delta,&phi,atom_col,atom_row,r,nfcol,nfrow,3,3);
		sigma=-sigma;
		delta=-delta;

		for(j=1;j<=7;j++)
		{
			for(k=1;k<=7;k++)
			{
				S[8+j][8+k]=ftr[j]*ftr[k]*sigma+(ftr[j+7]*ftr[k+7]+ftr[j+14]*ftr[k+14])*pi+(ftr[j+21]*ftr[k+21]+ftr[j+28]*ftr[k+28])*delta+(ftr[j+35]*ftr[k+35]+ftr[j+42]*ftr[k+42])*phi;
				S[8+k][8+j]=S[8+j][8+k];
			}
		}
	}


	for(j=0;j<16;j++)
	{
		if(j==0) 		VSIP_row=period[atom_row].orb[0].VSIP;
		if((3>=j)&&(j>0)) 	VSIP_row=period[atom_row].orb[1].VSIP;
		if((8>=j)&&(j>3)) 	VSIP_row=period[atom_row].orb[2].VSIP;
		if((15>=j)&&(j>8)) 	VSIP_row=period[atom_row].orb[3].VSIP;

		for(k=0;k<16;k++)
		{

			if(k==0) 		VSIP_col=period[atom_col].orb[0].VSIP;
			if((3>=k)&&(k>0)) 	VSIP_col=period[atom_col].orb[1].VSIP;
			if((8>=k)&&(k>3)) 	VSIP_col=period[atom_col].orb[2].VSIP;
			if((15>=k)&&(k>8)) 	VSIP_col=period[atom_col].orb[3].VSIP;

			H[j][k]=(VSIP_row+VSIP_col)*S[j][k]*KEHT/2.;
		}
	}

}



//////////////////////////////////////////////////////////////////////////////////////////////
//
//  DOES SOMETHING IMPORTANT
//
//////////////////////////////////////////////////////////////////////////////////////////////
void mov(double *sigma,double *pi,double *delta,double *phi,int atom_col,int atom_row,double rr,int n1,int n2,int l1,int l2)
{
	int i,nn,ia,ib,ik,il,ij,maxcal;
	double sk1,sk2,rll[4],xx,yy;
	double aa[30],bb[30],*a=aa-1,*b=bb-1;
	
	void abfns(double *,double *,double,double,double,int);

	double lovlap(double *,double *,double,double,double,int,int,int,int,int);

	nn=(l1<l2)?l1+1:l2+1;

	*sigma=0.0;
	*pi=0.0;
	*delta=0.0;
	*phi=0.0;

	ia=1;
	ib=1;
	maxcal=n1+n2;

	
	for(i=1;i<30;i++)
	{
		a[i]=0.0;
		b[i]=0.0;

	}
	for(i=0;i<3;i++)
	{
		rll[i]=0.0;
	}



	if (period[atom_col].orb[l1].coef2[1]!=0) ia=2;
	if (period[atom_row].orb[l2].coef2[1]!=0) ib=2;

	for(ik=0;ik<ia;ik++)
	{
		for(il=0;il<ib;il++)
		{
			sk1=period[atom_col].orb[l1].exp2[ik];

			sk2=period[atom_row].orb[l2].exp2[il];

			abfns(a,b,sk1,sk2,rr,maxcal);
			
			for(ij=0;ij<nn;ij++)
			{
				rll[ij]=lovlap(a,b,sk1,sk2,rr,l1,l2,ij,n1,n2);
			}

			xx=period[atom_col].orb[l1].coef2[ik];
			yy=period[atom_row].orb[l2].coef2[il];

			*sigma=*sigma+xx*yy*rll[0];
			*pi=*pi+xx*yy*rll[1];
			*delta=*delta+xx*yy*rll[2];
			*phi=*phi+xx*yy*rll[3];

		}
	}
}





//////////////////////////////////////////////////////////////////////////////////////////////
//
// THAT ONE TOO
//
//////////////////////////////////////////////////////////////////////////////////////////////
void abfns(double *a,double *b,double sk1,double sk2,double rr,int maxcal)
{
	int i,j,il,ix,ir,is,k,in;
	double rho1,rho2,c,d,h,r,ra,rho22,t,tr,temp;

	j=maxcal+1;

	rho1=.5*(sk1+sk2)*rr;
	rho2=.5*(sk1-sk2)*rr;

	if((rho1>165)||(rho2>165))
	{
		for(i=1;i<=20;i++)
		{
			a[i]=0;
			b[i]=0;
		}
		return;
	}

	c=exp(-rho1);
	a[1]=c/rho1;

	for(i=2;i<=j;i++)
		a[i]=((double)(i-1)*a[i-1]+c)/rho1;
	
	ix=j;
	ir=(rho2>0)? (2*rho2) : -(2*rho2);
	is=(ir+1<19)?(ir+1):19;


	if(rho2==0)
	{
		for(i=1;i<=ix;i=i+2)
		{
			b[i]=2.0/i;
			b[i+1]=0;
		}
		return;
	}

	d=exp(rho2);
	h=1/d;
	
	r=d-h;
	temp=(r>0)?r:-r;
	if(temp<.1)
	{
		ra=rho2;
		rho22=rho2*rho2;
		t=rho2;

		for(i=2;i<=50;i=i+2)
		{
			t=t*rho22/(double)(i*i+i);
			ra=ra+t;
			if(t<1e-30) break;
		}
		r=ra+ra;
		
	}
	
	b[1]=r/rho2;

	for(i=2;i<=ix;i=i+is)
	{
		if(ir!=0)
		{
			il=is-1;
			if(1<=il)
			{
				for(j=1;j<=il;j++)
				{
					k=i+j-1;
					if((float)k/2 !=(int)(k/2))
					{
						b[k]=(r+(double)(k-1)*b[k-1])/rho2;
					}
					else
					{
						b[k]=-(d+h-(double)(k-1)*b[k-1])/rho2;
					}
				}
			}
		}


		in=i+is-1;

		if(in-ix >0) 
			return;
			
		
		if((float)in/2 ==(int)(in/2)) 
		{
			tr=rho2;
			b[in]=-2.0*tr/(double)(in+1);

			for(j=1;j<=500;j++)
			{
				tr=tr*rho2*rho2/(double)((2*j)*(2*j+1));
				temp=tr/b[in];
				temp=(temp>0)?temp:-temp;
				if(temp-1.0e-7 <=0)
					break;
				b[in]=b[in]-2.0*tr/(double)(in+1+2*j);
			}
		}
		else
		{
			tr=1.0;

			b[in]=2.0*tr/(double)(in);

			for(j=1;j<=500;j++)
			{
				tr=tr*rho2*rho2/(double)((2*j)*(2*j-1));
				temp=tr/b[in];
				temp=(temp>0)?temp:-temp;
				if(temp-1.0e-7 <= 0) 
					break;
				b[in]=b[in]+2.0*tr/(double)(in+2*j);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//
//  EVEN MORE IMPORTANT
//
//////////////////////////////////////////////////////////////////////////////////////////////
double lovlap(double *a,double *b,double sk1,double sk2,double r,int l1,int l2,int m1,int n1,int n2)
{
	double fact[25];

	int m2,jend,kend,ieb,i,j,k,i1,i2,i3,i4,i5,i6,iev,iab,ibb,icb,idb,ju,ku,ip,ir;

	double strad=0.0,value,rhoa,rhob,rhoap,rhoab,rhopo,terma,term,value1,value2,value3,value4,con1,con12;

	int bincoe[8][8]={1,1,1,1,1,1,1,1, 0,1,2,3,4,5,6,7, 0,0,1,3,6,10,15,21, 0,0,0,1,4,10,20,35, 0,0,0,0,1,5,15,35, 0,0,0,0,0,1,6,21, 0,0,0,0,0,0,1,7, 0,0,0,0,0,0,0,1};


	fact[0]=1;
	
	for(i=1;i<25;i++)
	{
		fact[i]=i*fact[i-1];
	}


	m2=m1;
	rhoa=r*sk1;
	rhob=r*sk2;
	rhoap=pow(rhoa,2*n1+1);
	rhoab=pow(rhob,2*n2+1);
	rhopo=rhoap*rhoab;

        terma=pow(0.5,l1+l2+1)*sqrt(((l1+l1+1)*(l2+l2+1))*fact[l1-m1]*fact[l2-m1]/(fact[n1+n1]*fact[n2+n2]*fact[l1+m1]*fact[l2+m1])*rhopo);  


	jend=1+(int)((l1-m1)/2);
	kend=1+(int)((l2-m2)/2);
	ieb=m1+1;

	for(j=1;j<=jend;j++)
	{
		ju=j-1;
		iab=n1-l1+ju+ju+1;
		icb=l1-m1-ju-ju+1;
		con1=fact[l1+l1-ju-ju]/(double)(fact[l1-m1-ju-ju]*fact[ju]*fact[l1-ju]);
		for(k=1;k<=kend;k++)
		{
			ku=k-1;
			con12=con1*fact[l2+l2-ku-ku]/(double)(fact[l2-m2-ku-ku]*fact[ku]*fact[l2-ku]);
			iev=ju+ku+l2;
			if((int)(iev/2)!=(float)iev/2) con12=-con12;
			ibb=n2-l2+ku+ku+1;
			idb=l2-m2-ku-ku+1;
			value=0.0;
			for(i6=1;i6<=ieb;i6++)
			{
				for(i5=1;i5<=ieb;i5++)
				{
					value1=bincoe[i6-1][ieb-1]*bincoe[i5-1][ieb-1];
					iev=i5+i6; 
					if((int)(iev/2)!=(float)iev/2) value1=-value1;

					for(i4=1;i4<=idb;i4++)
					{
						value1=-value1;
						value2=bincoe[i4-1][idb-1]*value1;

						for(i3=1;i3<=icb;i3++)
						{
							value3=bincoe[i3-1][icb-1]*value2;

							for(i2=1;i2<=ibb;i2++)
							{
								value3=-value3;
								value4=bincoe[i2-1][ibb-1]*value3;

								for(i1=1;i1<=iab;i1++)
								{
									term=value4*bincoe[i1-1][iab-1];
									ir=i1+i2+ieb+ieb-i6-i6-i3+idb-i4+icb-1;
									ip=iab-i1+ibb-i2+2*ieb-2*i5+icb-i3+idb-i4+1;
									value=value+a[ip]*b[ir]*term;
								}
							}
						}
					}
				}
			}
			strad=strad+value*con12;
		}
	}               

	strad=strad*terma;
	return strad;
}


//////////////////////////////////////////////////////////
//	Read the atomic parameters
//////////////////////////////////////////////////////////

void read_atomic_parameters(char * HUCKEL_PARAM)
{
	FILE *Fparameter;
	int i,j;

	//Fparameter=fopen("/Users/nicolasrenaud/Documents/husky.1/SRC/huckel/parameters","r");
	Fparameter=fopen(HUCKEL_PARAM,"r");
  
	if(Fparameter == NULL)
	{
	  //printf("Huckel parameters not found at:  /home/nico/Bureau/husky.1/SRC/huckel/parameters \n");
    printf("Huckel parameters not found at:  %s \n",HUCKEL_PARAM);
	  exit(1);
	}


	for(i=1;i<=TABLESIZE;i++)
	{
		
		period[i].symbol[2]=0;
		
		fscanf(Fparameter,"%c%c",&period[i].symbol[0], &period[i].symbol[1]);
		fscanf(Fparameter,"%i",&period[i].valence_electron);
	
		for(j=0;j<4;j++)
		{
			fscanf(Fparameter,"%i ",&period[i].orb[j].orb_of_e);
		}
	
		for(j=0;j<4;j++)
		{
			fscanf(Fparameter,"%lf %lf %lf %lf %lf %lf",&period[i].orb[j].VSIP, &period[i].orb[j].exp, &period[i].orb[j].exp2[0], &period[i].orb[j].exp2[1], &period[i].orb[j].coef2[0], &period[i].orb[j].coef2[1]);
	
			period[6].orb[j].coef2[1]=0;

			if((period[i].orb[j].coef2[1]==0)||(j<2))
			{
				period[i].orb[j].exp2[0]=period[i].orb[j].exp;
				period[i].orb[j].coef2[0]=1;
				period[i].orb[j].coef2[1]=0;
			}
		}
		
	
		fscanf(Fparameter,"\n");

	}
	
	fclose(Fparameter);
}


//////////////////////////////////////////////////////////
//	Compute the nb or orbitals
//////////////////////////////////////////////////////////

int compute_nb_orb(char *file_name)
{
  int indexi,i,atom_row,noorb_i;
  indexi=0;
  int nb =0;
	
	double KEHT;
	atom *molecule;
	int no_atoms;
	char HKL_PARAM[250];
	
	
	
	// read the input file of the double cell
	read_molecule_input_file(file_name, &molecule, &no_atoms, &KEHT, HKL_PARAM);
  read_atomic_parameters(HKL_PARAM);
  
  for(i=0;i<no_atoms;i++)
  {
      atom_row=molecule[i].atomtype;
      
      // derermine the # of orbitals
      noorb_i=0;		
      
      if(period[atom_row].orb[0].orb_of_e!=0) noorb_i=noorb_i+1;
      if(period[atom_row].orb[1].orb_of_e!=0) noorb_i=noorb_i+3;
      if(period[atom_row].orb[2].orb_of_e!=0) noorb_i=noorb_i+5;
      if(period[atom_row].orb[3].orb_of_e!=0) noorb_i=noorb_i+7;
      
      nb = nb+noorb_i;
  }

  return(nb);
   
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      COMPUTE THE HUCKEL HAMILTONIAN OF THE TWO INTERACTING FRAGMENTS
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void compute_huckel_hamiltonian_general(double *hmat, double *smat, int nb_orb_precomp, char *file_name)
{

	// index for calculation
	int ii,jj,i,j,k,atom_row,atom_col,noorb_i,noorb_j,indexi,indexj;
	double delx,dely,delz;
	
	//temp matrices for the hamiltonian/overlap
	double S[16][16],H[16][16];
	
	// store the moleclar info
	atom *molecule;
	int no_atoms;
	
	// number of orbital for check
	int norb_tot = 0;
	
	// the HKL parameter
	double KEHT;
	
	// the adress of the parameter
	char HKL_PARAM[250];


	// read the input file of the double cell
	read_molecule_input_file(file_name, &molecule, &no_atoms, &KEHT, HKL_PARAM);

		
	// read the EHMO parameter provided by the user
	read_atomic_parameters(HKL_PARAM);


  // debug a bit
  if(0)
  {
	 //check data
	 printf("no_atoms=%d\n",no_atoms);
	 for(i=0;i<no_atoms;i++)
	   printf("%s\t%d\t%f\t%f\t%f\n",molecule[i].atomTypeChar,molecule[i].atomtype,molecule[i].x,molecule[i].y,molecule[i].z);
  }

  // compute the Hamoltonian and overlap
	indexi=0;

	for(i=0;i<no_atoms;i++)	
	{
		indexj=0;
		atom_row=molecule[i].atomtype;
		
		// derermine the # of orbitals
		noorb_i=0;		
		if(period[atom_row].orb[0].orb_of_e!=0) noorb_i=noorb_i+1;
		if(period[atom_row].orb[1].orb_of_e!=0) noorb_i=noorb_i+3;
		if(period[atom_row].orb[2].orb_of_e!=0) noorb_i=noorb_i+5;
		if(period[atom_row].orb[3].orb_of_e!=0) noorb_i=noorb_i+7;
		
		
		// increment the # oftotal orbitals  
		norb_tot += noorb_i;  		  
		
		for(j=0;j<no_atoms;j++)
		{
			
			delx=molecule[j].x-molecule[i].x;
			dely=molecule[j].y-molecule[i].y;
			delz=molecule[j].z-molecule[i].z;
			atom_col=molecule[j].atomtype;
			noorb_j=0;
			if(period[atom_col].orb[0].orb_of_e!=0) noorb_j=noorb_j+1;
			if(period[atom_col].orb[1].orb_of_e!=0) noorb_j=noorb_j+3;
			if(period[atom_col].orb[2].orb_of_e!=0) noorb_j=noorb_j+5;
			if(period[atom_col].orb[3].orb_of_e!=0) noorb_j=noorb_j+7;

			for(ii=0;ii<16;ii++)
			{
				for(jj=0;jj<16;jj++)
				{
					S[ii][jj]=0;
					H[ii][jj]=0;
				}
			}

			if(i==j)
			{
				for(ii=0;ii<16;ii++)
				{
					S[ii][ii]=1;
					if(ii==0) H[ii][ii]=period[atom_col].orb[0].VSIP;
					if((3>=ii)&&(ii>0)) H[ii][ii]=period[atom_col].orb[1].VSIP;
					if((8>=ii)&&(ii>3)) H[ii][ii]=period[atom_col].orb[2].VSIP;
					if((15>=ii)&&(ii>8)) H[ii][ii]=period[atom_col].orb[3].VSIP;
				}
			}
			else
			{
				overlap(atom_row,atom_col,delx,dely,delz,S,H,KEHT);
			}


			for(ii=0;ii<noorb_i;ii++)
			{
				for(jj=0;jj<noorb_j;jj++)
				{
					//smat[indexi+ii][indexj+jj]=S[ii][jj];
					//hmat[indexi+ii][indexj+jj]=H[ii][jj];
					
					// replace that with 1D array
					smat[ (indexi+ii)*nb_orb_precomp + (indexj+jj) ]=S[ii][jj];
					hmat[ (indexi+ii)*nb_orb_precomp + (indexj+jj) ]=H[ii][jj];
					
				}	
			}

			indexj=indexj+noorb_j;

		}

		indexi=indexi+noorb_i;

	}

	no_orbitals=indexi;

/************************************************/
	
	
	// little check
	if(no_orbitals != nb_orb_precomp)
		printf(" warning precomputed %d orbitals and found %d\n",nb_orb_precomp,no_orbitals);
	
/************************************************/

}







