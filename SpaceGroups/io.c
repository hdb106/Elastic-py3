#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sto.h"
#include "type_sg.h"
#include "math_sg.h"
#include "pgrp_dat.h"

void find_atom_pgroup( int lat, int nop_s, int npgrp_s, double sym_op[][4][3], double Rsh[3],
		       int *npgrp_p, double op[][4][3], /*  number and operations of point group  */
		       double T[3][3]); /*  orientation of new elementary basis  */

#define ERROR { fprintf(stderr,"Error: not completed input!"); fclose(fl); return(1); }

/* read data from WIEN input file */
int read_wien_data(char *fname, t_cell *cl_in, char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
	      int *lat, int *is_wien_rhomb)
{
   FILE *fl;
   double d,d1,Z_W;
   char s[256],st[64];
   int i,j,k,l,nat_old=0, NATO_W, MULT_W;
   
   if( (fl=fopen(fname,"rt")) == NULL ) {
      fprintf(stderr,"Can't open file: %s\n.",fname);
      return(0);
   }
   
  *is_wien_rhomb=0;
  /*  read title  */
  if( fgets(s,256,fl)==NULL ) ERROR
  /*  read lattice type  */
  if( fgets(s,256,fl)==NULL ) ERROR
  strncpy(st,s,4); st[4]='\0';
  if( strstr(st,"P") != NULL ) *lat=ORTHOROMBIC_P;
  else
  if( strstr(st,"F") != NULL ) *lat=ORTHOROMBIC_F;
  else
  if( strstr(st,"B") != NULL ) *lat=ORTHOROMBIC_I;
  else
  if( strstr(st,"CXY") != NULL ) *lat=ORTHOROMBIC_C;
  else
  if( strstr(st,"CYZ") != NULL ) *lat=ORTHOROMBIC_A;
  else
  if( strstr(st,"CXZ") != NULL ) *lat=ORTHOROMBIC_B;
     /* WIEN uses the reverse setting for rhombohedral */
  if( strstr(st,"R") != NULL ) { *lat=RHOMBOHEDRAL; *is_wien_rhomb=1; }
  if( strstr(st,"H") != NULL ) *lat=HEXAGONAL;
  /*  NATO param. of WIEN  */
  strncpy(st,s+4+24,2); st[2]='\0';
  NATO_W=atoi(st);
  /*  read a,b,c alpha, beta, gamma  */
  if( fgets(s,256,fl)==NULL ) ERROR   /*   skip 1  */
  if( fgets(s,256,fl)==NULL ) ERROR
  
  j=strlen(s);
  strncpy(st,s,10);    st[10]='\0';  a[0][0]=atof(st);
  strncpy(st,s+10,10); st[10]='\0';  a[0][1]=atof(st);
  strncpy(st,s+20,10); st[10]='\0';  a[0][2]=atof(st);
  if( j > 30 ) strncpy(st,s+30,10); st[10]='\0';  if( stod(st,&(a[1][0])) != 1 ) a[1][0]=90.;
  if( j > 40 ) strncpy(st,s+40,10); st[10]='\0';  if( stod(st,&(a[1][1])) != 1 ) a[1][1]=90.;
  if( j > 50 ) strncpy(st,s+50,10); st[10]='\0';  if( stod(st,&(a[1][2])) != 1 ) a[1][2]=90.;
   
  switch( *lat ) {
   case ORTHOROMBIC_I:
      d=a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2];
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos(( a[0][0]-a[0][1]-a[0][2])/d);
      a[1][1]=180./PI*acos((-a[0][0]+a[0][1]-a[0][2])/d);
      a[1][2]=180./PI*acos((-a[0][0]-a[0][1]+a[0][2])/d);
      a[0][0]=a[0][1]=a[0][2]=sqrt(d)/2;
      break;
   case ORTHOROMBIC_F:
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos( a[0][0]/sqrt((a[0][0]+a[0][1])*(a[0][0]+a[0][2])) );
      a[1][1]=180./PI*acos( a[0][1]/sqrt((a[0][0]+a[0][1])*(a[0][1]+a[0][2])) );
      a[1][2]=180./PI*acos( a[0][2]/sqrt((a[0][0]+a[0][2])*(a[0][1]+a[0][2])) );
      d      =sqrt(a[0][1]+a[0][2])/2;
      d1     =sqrt(a[0][0]+a[0][2])/2;
      a[0][2]=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=d; a[0][1]=d1;
      break;
   case ORTHOROMBIC_C:
/*        a[1][0]=a[1][1]=a[1][2]=90.; WIEN rule for input  */
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1];
      a[1][0]=a[1][1]=90.;
      a[1][2]=180./PI*acos((a[0][0]-a[0][1])/(a[0][0]+a[0][1]));
      d=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=a[0][1]=d;
      break;
   case ORTHOROMBIC_A:
/*      a[1][0]=a[1][1]=a[1][2]=90.; WIEN rule for input  */
      a[1][0]=180./PI*acos((-a[0][1]*a[0][1]+a[0][2]*a[0][2])/
			   ( a[0][1]*a[0][1]+a[0][2]*a[0][2]));
      a[1][1]=180./PI*acos(-a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      a[1][2]=180./PI*acos( a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      d=sqrt(a[0][1]*a[0][1] + a[0][2]*a[0][2])/2;
      a[0][1]=a[0][2]=d;
      break;
   case ORTHOROMBIC_B:
/*        a[1][0]=a[1][1]=90.; WIEN rule for input  */
      d=180./PI*acos(a[0][0]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][0]*a[0][0]+a[0][2]*a[0][2]) ) );
      d1=180./PI*acos((a[0][0]*a[0][0] - a[0][2]*a[0][2])/
			   ( a[0][0]*a[0][0]+a[0][2]*a[0][2]));
      a[1][0]=a[1][2]=d;
      a[1][1]=d1;
      d=sqrt(a[0][0]*a[0][0] + a[0][2]*a[0][2])/2;
      a[0][0]=a[0][2]=d;
      break;
   case RHOMBOHEDRAL:
      /*  A_rhomb^2 = A_hex^2 /3 + C_hex^2 /9  */
      d=sqrt( a[0][0]*a[0][0]/3 + a[0][2]*a[0][2]/9 );
     /*  Cos[alpha]_rhomb = (-A_hex^2 /6 + C_hex^2 /9)/A_rhomb^2;  */
      d1=180./PI*acos((-a[0][0]*a[0][0]/6 + a[0][2]*a[0][2]/9)/(d*d));
      a[0][0]=a[0][1]=a[0][2]=d;
      a[1][0]=a[1][1]=a[1][2]=d1;
     /*  WIEN input supposes the primitive (rhombohedral) coordinates of atoms  */
      *lat=ORTHOROMBIC_P;
      break;
   case HEXAGONAL:
      a[1][0]=a[1][1]=90.; a[1][2]=120.;
      *lat=ORTHOROMBIC_P;
   }
     
  cl_in->nat=cl_in->nsrt=0;
  for(i=0; i<NATO_W; i++) {
    if( fgets(s,256,fl)==NULL ) ERROR
    strncpy(st,s+5+3+4,10);           st[10]='\0';  cl_in->r[cl_in->nat][0]=atof(st);
    strncpy(st,s+5+3+4+10+3,10);      st[10]='\0';  cl_in->r[cl_in->nat][1]=atof(st);
    strncpy(st,s+5+3+4+10+3+10+3,10); st[10]='\0';  cl_in->r[cl_in->nat][2]=atof(st);
    cl_in->nat++;
    if( fgets(s,256,fl)==NULL ) ERROR /*  read MULT  */
    strncpy(st,s+15,2);  st[2]='\0';  MULT_W=atoi(st);
    for(j=0; j<MULT_W-1; j++) {
      if( fgets(s,256,fl)==NULL ) ERROR
      strncpy(st,s+5+3+4,10);           st[10]='\0';  cl_in->r[cl_in->nat][0]=atof(st);
      strncpy(st,s+5+3+4+10+3,10);      st[10]='\0';  cl_in->r[cl_in->nat][1]=atof(st);
      strncpy(st,s+5+3+4+10+3+10+3,10); st[10]='\0';  cl_in->r[cl_in->nat][2]=atof(st);
      cl_in->nat++;
    }
    if( fgets(s,256,fl)==NULL ) ERROR /*  read Z  */
    strncpy(st,s+10+5+5+5+10+5+10+5,5);  st[5]='\0';  Z_W=atof(st);
    /*  find this Z value  */
    sprintf(name_srt[cl_in->nsrt],"z=%.2f",Z_W);
    for(k=0; k<cl_in->nsrt; k++)
    if( !strcmp(name_srt[cl_in->nsrt],name_srt[k]) ) { /*  this sort already exists  */
	 for(l=nat_old; l < cl_in->nat; l++ ) cl_in->srt[l]=k;
         goto NEXT;
    }
    /*  new sort  */
    for(l=nat_old; l < cl_in->nat; l++ ) cl_in->srt[l]=cl_in->nsrt;
    cl_in->nsrt++;
    NEXT:
    nat_old=cl_in->nat;
    if( fgets(s,256,fl)==NULL ) ERROR /*  read ROT. matrix  */
    if( fgets(s,256,fl)==NULL ) ERROR /*  read ROT. matrix  */
    if( fgets(s,256,fl)==NULL ) ERROR /*  read ROT. matrix  */
  } /*  for(i=0; i<NATO_W; i++) {  */
  return(0);
}

int read_data(char *fname, t_cell *cl_in, char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
	      int *lat)
/*  "lat" here means only the type of centring mode (P,F,I,C,A),
 therefore "lat" will take only values: ORTHOROMBIC_P, ORTHOROMBIC_F, ORTHOROMBIC_I,
 * ORTHOROMBIC_A, ORTHOROMBIC_C
 */
{
   FILE *fl;
   double d,d1;
   char s[256],st[64];
   int i,j;

   if( (fl=fopen(fname,"rt")) == NULL ) return 0;
   *lat=ORTHOROMBIC_P;
/*  read type of cell  */
   do
     if( fgets(s,256,fl)==NULL ) ERROR
   while( stos(s,st)==0 );
/*  read translation vectors  */
   do
     if( fgets(s,256,fl)==NULL ) ERROR
   while( stod(s,&a[0][0])==0 );
   if( st[0]=='i' || st[0]=='I') {
      *lat=ORTHOROMBIC_I;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
	 fprintf(stderr,"WARNING: for I-centred type angles must be: "
		"alpha=beta=gamma=90.\nthis relation will be assumed.\n" );
      d=a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2];
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos(( a[0][0]-a[0][1]-a[0][2])/d);
      a[1][1]=180./PI*acos((-a[0][0]+a[0][1]-a[0][2])/d);
      a[1][2]=180./PI*acos((-a[0][0]-a[0][1]+a[0][2])/d);
      a[0][0]=a[0][1]=a[0][2]=sqrt(d)/2;
   } else
   if( st[0]=='f'|| st[0]=='F') {
      *lat=ORTHOROMBIC_F;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
	 fprintf(stderr,"WARNING: for F-centred type angles must be: "
		"alpha=beta=gamma=90.\nthis relation will be assumed.\n" );
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos( a[0][0]/sqrt((a[0][0]+a[0][1])*(a[0][0]+a[0][2])) );
      a[1][1]=180./PI*acos( a[0][1]/sqrt((a[0][1]+a[0][0])*(a[0][1]+a[0][2])) );
      a[1][2]=180./PI*acos( a[0][2]/sqrt((a[0][2]+a[0][0])*(a[0][2]+a[0][1])) );
      d      =sqrt(a[0][1]+a[0][2])/2;
      d1     =sqrt(a[0][0]+a[0][2])/2;
      a[0][2]=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=d; a[0][1]=d1;
   } else
   if( st[0]=='c'|| st[0]=='C') { /*  C-centred: only orthorombic, cubic, tetragonal  */
      *lat=ORTHOROMBIC_C;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
	 fprintf(stderr,"WARNING: for C-centred type angles must be: "
		"alpha=beta=gamma=90.\nthis relation will be assumed.\n" );
      a[1][0]=a[1][1]=90.;
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1];
      a[1][2]=180./PI*acos((a[0][0]-a[0][1])/(a[0][0]+a[0][1]));
      d=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=a[0][1]=d;
   } else
   if( st[0]=='a'|| st[0]=='A') { /*  A-centred: orthorombic and monoclinic  */
      *lat=ORTHOROMBIC_A;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL)
          fprintf(stderr,"WARNING: for A-centred type angles must be: "
                  "alpha=beta=90.\nthis relation will be assumed.\n");
      a[1][0]=180./PI*acos((-a[0][1]*a[0][1]+a[0][2]*a[0][2])/
			   ( a[0][1]*a[0][1]+a[0][2]*a[0][2]));
      a[1][1]=180./PI*acos(-a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      a[1][2]=180./PI*acos( a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      d=sqrt(a[0][1]*a[0][1] + a[0][2]*a[0][2])/2;
      a[0][1]=a[0][2]=d;
   }
/*  read number of atoms  */
    do
     if( fgets(s,256,fl)==NULL ) ERROR
    while( stoi(s,&cl_in->nat)==0 );
    if( cl_in->nat > MAX_ATOMS ) {
       fprintf(stderr,"Parameter MAX_ATOMS=%d is too small! Increase it.\n", MAX_ATOMS);
       fclose(fl);
       return(1);
    }
/*  read coordinates and types  */
    cl_in->nsrt=0;
    for(i=0; i < cl_in->nat; i++)
    {
       do
         if( fgets(s,256,fl)==NULL ) ERROR     	
      while( stod(s,cl_in->r[i])==0 );

       do
         if( fgets(s,256,fl)==NULL ) ERROR     	
       while( stos(s,st)==0 );
       for(j=0; j <= cl_in->nsrt; j++){
	  if(j == cl_in->nsrt) {
	     if( cl_in->nsrt == MAX_ATOMS ) {
		fprintf(stderr,"Parameter MAX_ATOMS=%d is too small! Increase it.\n", MAX_ATOMS);
		fclose(fl);
		return(1);
	     }
	     strcpy(name_srt[j],st);
	     cl_in->nsrt++;
	     cl_in->srt[i]=j;
	     break;
	  } else
          if( !strcmp(st,name_srt[j]) ) { /*  this sort already exists  */
	     cl_in->srt[i]=j;
	     break;
	  }
       }
    }
    fclose(fl);
    return(0);
}

int write_res(char *fname, t_cell *cl, char name_srt[MAX_ATOMS][MAX_CHARS],
	      t_cell *cl_new, char name_srt_new[MAX_ATOMS][MAX_CHARS],
	      int is_noeq, int is_prim, int fix_org, int is_wien_rhomb,
	      int nsh, double Rsh[][3],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int lat, char lat_name[14][32], char * sgrp_name,
              int npgrp, int nop, double sym_op[][4][3],
              char name_pgrp[][32],
	      double Tin[3][3])
{
  FILE *fl;

  double r[3],Tel[3][3], pgrp_at[48][4][3];
  int i, j, npgrp_at, nat_srt;

  if( fname == NULL) fl=stdout;  else
  if( (fl=fopen(fname,"wt")) == NULL ) {
     fprintf(stderr,"Can not open file:%s\n", fname);
     return(1);
  }

  if( is_prim ) {
     fprintf(fl,"NOTE: atom positions and space group operations "
	         "are given in primitive basis.\n\n");
     fprintf(fl,"Bravais lattice : %s\n\nParameters of primitive cell:\n",lat_name[lat]);
  }
  else {
     if( lat == RHOMBOHEDRAL )
       fprintf(fl,"Bravais lattice: %s [hexagonal setting]\n\n",lat_name[lat]);
     else
       fprintf(fl,"Bravais lattice: %s\n\n",lat_name[lat]);
  }
  fprintf(fl,"     a             b            c\n");
  fprintf(fl,"% .8f   % .8f  % .8f\n",a[0][0],a[0][1],a[0][2]);
  fprintf(fl,"    alpha          beta          gamma\n");
  fprintf(fl,"% .8f   % .8f  % .8f\n\n",a[1][0],a[1][1],a[1][2]);
  if( !is_wien_rhomb)
    fprintf(fl,"\n===== Decomposition of new basis vectors over input basis =====\n");
  else
    fprintf(fl,"\n===== Decomposition of new basis vectors over input RHOMBOHEDRAL"
	       "\n===== (not HEXAGONAL!) basis\n");
  fprintf(fl,"% .6f  % .6f % .6f  <--- 1\n",Tin[0][0],Tin[1][0],Tin[2][0]);
  fprintf(fl,"% .6f  % .6f % .6f  <--- 2\n",Tin[0][1],Tin[1][1],Tin[2][1]);
  fprintf(fl,"% .6f  % .6f % .6f  <--- 3\n",Tin[0][2],Tin[1][2],Tin[2][2]);
/*    fprintf(fl,"\n=======================================\n");  */
   
  if( is_noeq && !is_prim ) {
    fprintf(fl,"\n==== Number of atoms in cell (only atoms of primitive cell): %d \n",cl->nat);
    fprintf(fl,"==== Atom positions (only atoms of primitive cell):\n\n");
  } else {
    fprintf(fl,"\n==== Number of atoms in cell: %d\n",cl->nat);
    fprintf(fl,"==== Atom positions:\n\n");
  }
  for(i=0; i<cl->nat; i++) {
    if(i &&  cl->srt[i-1] != cl->srt[i] ) fprintf(fl,"\n");
    fprintf(fl,"% .8f  % .8f % .8f\n",cl->r[i][0],cl->r[i][1],cl->r[i][2]);
    fprintf(fl," %s\n",name_srt[cl->srt[i]]);
  }
/*    fprintf(fl,"\n======================================================================\n");  */
  if( !is_prim)
/*    fprintf(fl,"Nonequivalent atoms, point group and rotation matrix for each sort:\n");  */
  fprintf(fl,"\n==== Nonequivalent atoms, point group for each sort: ====\n");
  else
  fprintf(fl,"\n==== Nonequivalent atoms for each sort: ====\n");
  for(i=0; i<cl_new->nat; i++) {
      /*  find point group for nonequivalent atoms  */
/*     if( ! is_prim )
     if( cl_new->srt[i] == isrt ) {
      isrt++;
      find_atom_pgroup(lat,nop,npgrp,sym_op,cl_new->r[i],&npgrp_at,pgrp_at,Tel);
      fprintf(fl,"Names of point group are: %s\n", name_pgrp[npgrp_at]);
      fprintf(fl,"Rotation matrix for this point group:\n");
      for(j=0; j<3; j++)
         fprintf(fl,"% .4f  % .4f % .4f\n",Tel[j][0],Tel[j][1],Tel[j][2]);
     }*/
     
    if( i==0 || cl_new->srt[i-1] != cl_new->srt[i] ) {
       fprintf(fl, "\nSort number: %d\n",cl_new->srt[i]+1);
       if( ! is_prim ) {
          find_atom_pgroup(lat,nop,npgrp,sym_op,cl_new->r[i],&npgrp_at,pgrp_at,Tel);
          fprintf(fl,"  Names of point group: %s\n", name_pgrp[npgrp_at]);
          fprintf(fl,"  New basis vectors for this point group:\n");
	  fprintf(fl,"   % .4f  % .4f % .4f  <--- 1\n",Tel[0][0],Tel[1][0],Tel[2][0]);
	  fprintf(fl,"   % .4f  % .4f % .4f  <--- 2\n",Tel[0][1],Tel[1][1],Tel[2][1]);
	  fprintf(fl,"   % .4f  % .4f % .4f  <--- 3\n",Tel[0][2],Tel[1][2],Tel[2][2]);
          fprintf(fl,"\n");
       }
       /* count number of atoms of a given sort */
       for(j=nat_srt=0; j<cl_new->nat; j++) 
	 if( cl_new->srt[j] == cl_new->srt[i] ) nat_srt++;

       fprintf(fl,"  Atom positions: %d\n", nat_srt);
    }
     
    fprintf(fl,"  % .8f  % .8f % .8f\n",cl_new->r[i][0],
	   cl_new->r[i][1],cl_new->r[i][2]);
    fprintf(fl,"   %s\n",name_srt_new[cl_new->srt[i]]);
 
/*    fprintf(fl,"          X=%.8f Y=%.8f Z=%.8f\n",cl_new->r[i][0],
	   cl_new->r[i][1],cl_new->r[i][2]);
    fprintf(fl," %s\n",name_srt_new[cl_new->srt[i]]);*/

  }
  fprintf(fl,"\n======================================================================\n");
  if( lat==RHOMBOHEDRAL && !is_prim )
     fprintf(fl,"\nNumber and name of space group: %s [h axes]\n",sgrp_name);
   else
     fprintf(fl,"\nNumber and name of space group: %s\n",sgrp_name);
   
  fprintf(fl,"- Short - Full - Schoenflies - names of point group:\n %s\n",
	  name_pgrp[npgrp]);
  fprintf(fl,"\nNumber of symmetry operations: %d\n",nop);
  for(i=0; i < nop; i++) {
       fprintf(fl,"Operation: %d\n",i+1);
       for(j=0; j<3; j++)
           if( fabs(sym_op[i][3][j]-1./6) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  1/6\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
	   else
           if( fabs(sym_op[i][3][j]-1./3) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  1/3\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
	   else
           if( fabs(sym_op[i][3][j]-2./3) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  2/3\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
	   else
           if( fabs(sym_op[i][3][j]-5./6) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  5/6\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
	   else
           fprintf(fl,"% .1f  % .1f  % .1f % .3f\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2],sym_op[i][3][j]);
       fprintf(fl,"\n");
   }

     fprintf(fl,"============================================================");
     if( fix_org == 001 ) {
	fprintf(fl,"\nNote that shift vectors for this space group are defined\n"
	       "only up to the vector ");
	if( is_prim ) {
	   if( lat==CUBIC_P || lat==HEXAGONAL ||
	       lat==TETRAGONAL_P || lat==ORTHOROMBIC_P || lat==MONOCLINIC_P)
	       fprintf(fl,"{ 0, 0, Z }."); else
	   if( lat==RHOMBOHEDRAL )
	       fprintf(fl,"{ Z, Z, Z }."); else	
           if( lat==CUBIC_I || lat==TETRAGONAL_I || lat== ORTHOROMBIC_I)
	       fprintf(fl,"{ Z, Z, 0 }."); else
           if( lat==CUBIC_F || lat==ORTHOROMBIC_F )
	       fprintf(fl,"{ Z, Z, -Z }."); else
           if( lat==ORTHOROMBIC_C )
	       fprintf(fl,"{ 0, 0, Z }."); else
           if( lat==ORTHOROMBIC_A || lat==MONOCLINIC_A )
	       fprintf(fl,"{ 0, Z, Z }.");
	}
	else fprintf(fl,"{ 0, 0, Z }.");	
	fprintf(fl,"\nHere Z can take any value.");	
	
     } else
     if( fix_org == 110 ) {
	fprintf(fl,"\nNote that shift vectors for this space group are defined\n"
	       "only up to the vector ");
	if( is_prim ) {
	   if( lat==CUBIC_P || lat==HEXAGONAL ||
	       lat==TETRAGONAL_P || lat==ORTHOROMBIC_P || lat==MONOCLINIC_P)
	       fprintf(fl,"{ X, Y, 0 }."); else
	   if( lat==RHOMBOHEDRAL )
	       fprintf(fl,"{ -Y, X, -X+Y }."); else	
           if( lat==CUBIC_I || lat==TETRAGONAL_I || lat== ORTHOROMBIC_I)
	       fprintf(fl,"{ Y, X, X+Y }."); else
           if( lat==CUBIC_F || lat==ORTHOROMBIC_F )
	       fprintf(fl,"{ -X+Y, X-Y, X+Y }."); else
           if( lat==ORTHOROMBIC_C )
	       fprintf(fl,"{ X-Y, X+Y, 0 }."); else
           if( lat==ORTHOROMBIC_A || lat==MONOCLINIC_A )
	       fprintf(fl,"{ X, Y, -Y }.");
	}
	else fprintf(fl,"{ X, Y, 0 }.");	
	fprintf(fl,"\nHere X and Y can take any values.");
     } else
     if( fix_org == 111 )
	fprintf(fl,"\nNote that shift vector for this space group is defined\n"
	            "only up to the vector { X, Y, Z }.\nHere X, Y and Z can take any values.");

     fprintf(fl,"\n\n===== Number of possible shift vectors: %d =====\n",nsh);
     fprintf(fl,"===== List of shift vectors:\n");
     for(i=0; i<nsh; i++)
        fprintf(fl,"% .4f  % .4f  % .4f\n",Rsh[i][0], Rsh[i][1], Rsh[i][2]);
   if( nsh > 1 ) {
      fprintf(fl,"\nList of shifted cells:\n");
      for(j=1; j<nsh; j++) {
	 fprintf(fl,"\nCell #%d \n",j);
	 fprintf(fl,"  Shift vector: % .4f  % .4f  % .4f\n",Rsh[j][0], Rsh[j][1], Rsh[j][2]);
	 fprintf(fl,"  Nonequivalent atoms :\n");
	 for(i=0; i<cl_new->nat; i++) {
            if(i &&  cl_new->srt[i-1] != cl_new->srt[i] ) fprintf(fl,"\n");
            sub_vv(r,cl_new->r[i],Rsh[j]); reduce(r);
	    fprintf(fl,"% .8f  % .8f % .8f\n",r[0], r[1],r[2]);
	    fprintf(fl," %s\n",name_srt_new[cl_new->srt[i]]);
	 }
      }
   }
   return(0);
}

