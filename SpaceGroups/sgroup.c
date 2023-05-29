#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "type_sg.h"
#include "sto.h"
#include "io.h"
#include "math_sg.h"
#include "lat.h"

#include "pgrp_dat.h"
#include "pgrp_op.h"
#include "pgrp.h"
#include "sgrp.h"

/*  #define  DEBUG 0  */
#include "rotb.h"

char lat_name[15][32]={
     "Triclinic",
     "Monoclinic primitive","Monoclinic A bace-centred",
     "Orthorombic primitive",   "Orthorombic body-centred",
     "Orthorombic C bace-centred","Orthorombic A bace-centred",
     "Orthorombic all-face centred",
     "Tetragonal primitive", "Tetragonal body-centred",
     "Rhombohedral", "Hexagonal",
     "Cubic primitive", "Cubic body centred", "Cubic face-centred"
};

char *help_str=
"\nDESCRIPTION\n"
"=============\n"
"     sgroup [options] input_file [output_file]\n"
"\n"
"  OPTIONS:\n"
"     -noeq  Output contains only atoms not connected by lattice translations,\n"
"            i.e. only atoms of a primitive cell. For example, in the case of pure\n"
"            fcc lattice only (0,0,0) position will be printed, but not the\n"
"            positions (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5).\n"
"\n"	
"     -prim  Use the basis of a primitive cell for output.\n"
"\n"
"     -wien  Read data from a file written to be input for the WIEN package.\n"
"            That's WIEN's case.struct file.\n"
"\n"
"     -help  Print this message and exit.\n";
  
/**********************************************************************
 * generate T matrix; from double  a[2][3] - a,b,c, alpha, beta, gama
 * Ti-column is decomposition of Ti vector in cubic basis
 * ********************************************************************/

void T_matr(double T[3][3], double a[2][3])
{
   int i,j;

   for(i=0; i<2; i++)
   for(j=0; j<3; j++) {
      if( (i==0 && fabs(a[i][j]) < 1e-5) ||
          (i==1 && fabs(sin(PI/180.*a[i][j])) < 1e-5) )  {
           fprintf(stderr,"error: degenerated basis\n");
           exit(1);
      }
   }
   T[0][0]=a[0][0]; T[0][1]=0.; T[0][2]=0.;
   
   T[1][0]=a[0][1]*cos(PI/180.*a[1][2]); /*  b*cos(gamma)  */
   T[1][1]=a[0][1]*sin(PI/180.*a[1][2]); /*  b*sin(gamma)  */
   T[1][2]=0.;
   
   T[2][0]=a[0][2]*cos(PI/180.*a[1][1]); /*  c*cos(beta)  */
   /*  cy = c*(cos(al)-cos(bet)*cos(gam))/sin(gam)  */
   T[2][1]=a[0][2]*(cos(PI/180.*a[1][0])-cos(PI/180.*a[1][1])*
                    cos(PI/180.*a[1][2]))/sin(PI/180.*a[1][2]);
   
   T[2][2]=sqrt(a[0][2]*a[0][2] - T[2][0]*T[2][0] -
                T[2][1]*T[2][1]);
   for(i=0; i<9; i++) if( fabs(T[0][i])<TOL ) T[0][i]=0.;
   transp(T,T);
}

/********************************************************
 * transform atom coordinates from old basis to new one
 * Ti column is decompostion of i new basis vector
 * T_1=T^-1; r'=A^-1*r
 *******************************************************/
void transform_coor(t_cell *cl_in, double T_1[3][3])
{
   int i;
   for(i=0; i < cl_in->nat; i++)
     mul_mv(cl_in->r[i],T_1,cl_in->r[i]);
   return;
}
     
/******************************************************
 * find sort corresponding to minimal number of atoms
******************************************************/
int find_min_sort( double Rshft[3], t_cell *cl_in)
{
      int i, min_srt, isrt, nat, nat_prev=-1, at_shft;

      for(isrt=0; isrt < cl_in->nsrt; isrt++) {
	 /*  count atoms of sort isrt  */
	 for(i=nat=0; i < cl_in->nat; i++)
	 if( cl_in->srt[i] == isrt ) {
	    nat++;
	    at_shft=i;
	 }
	
	 if( nat_prev == -1 || nat_prev > nat) {
	    nat_prev=nat;
	    min_srt=isrt;
	    asgn(Rshft,cl_in->r[at_shft]);
	 }
      }
      return(min_srt);
}

/*************************************
 * shift origin of coordinate system
 *************************************/
void  shift_cell(double R[3], t_cell *cl)
{
   int i;
   double Rshft[3];
   
   asgn(Rshft,R);
   for(i=0; i < cl->nat; i++) {
      sub_vv(cl->r[i],cl->r[i],Rshft);
      reduce(cl->r[i]);
   }
}
     
/***************************
 * form cell for one sort
 ***************************/
void make_cell_for_sort(int srt, t_cell *cell, t_cell *cell_in)
{
   int atc;

   cell->nat=0;
   for(atc=0; atc < cell_in->nat; atc++)
   if( cell_in->srt[atc] == srt ) {
      asgn(cell->r[cell->nat],cell_in->r[atc]);
      cell->nat++;
   }
   return;
}

/*** build new cell, all input coordinates in kartesian basis ***/
void build_cell(t_cell *cl, t_cell *cln, double Tp[3][3])
{
   int i,at;
   double Tp_1[3][3];
   
   cln->nsrt=cl->nsrt;
   cln->nat=0;
   inv_matr(Tp,Tp_1);
   for(at=0; at < cl->nat; at++) {
      asgn(cln->r[cln->nat],cl->r[at]);
      cln->srt[cln->nat]=cl->srt[at];
   /*  transform from cubic into prim. basis and move to cell  */
      mul_mv(cln->r[cln->nat],Tp_1,cln->r[cln->nat]);
      reduce(cln->r[cln->nat]);
      /*  find this atom  */
      for(i=0; i < cln->nat; i++)
        if( cln->srt[i] == cl->srt[at] &&
            ddvec(cln->r[cln->nat],cln->r[i]) < TOL
          ) goto FOUND;
/*        printf("% .3f % .3f % .3f\n",cln->r[cln->nat][0],  */
/*               cln->r[cln->nat][1],cln->r[cln->nat][2]);  */
      cln->nat++;
      FOUND: ;
   }
}

/* coor. of cl_in are given in  prim. basis,
 * T[3][3] decomposition  of new trans. vectors (elem. basis)
 * rold=T*rnew ==> rnew=T^-1*rold
 * is_noeq=1  construct elementary cell containing only atoms of primitive one
 * is_noeq=0  construct full elementary cell */
void build_large_cell(t_cell *cl, t_cell *cl_out, double T[3][3], int is_noeq)
{
   int i,i1,i2,i3,at,n;
   double T_1[3][3];
   
   if( is_noeq==0 ) n=3; else n=1;
   cl_out->nsrt=cl->nsrt;
   cl_out->nat=0;
   inv_matr(T,T_1);
   for(at=0; at < cl->nat; at++) {
      for(i1=0; i1<n; i1++)
      for(i2=0; i2<n; i2++)
      for(i3=0; i3<n; i3++) {
         cl_out->srt[cl_out->nat]=cl->srt[at];
         cl_out->r[cl_out->nat][0]=cl->r[at][0]+i1;
         cl_out->r[cl_out->nat][1]=cl->r[at][1]+i2;
         cl_out->r[cl_out->nat][2]=cl->r[at][2]+i3;
         mul_mv(cl_out->r[cl_out->nat], T_1, cl_out->r[cl_out->nat]);
         reduce(cl_out->r[cl_out->nat]);
        /*  find this atom  */
         for(i=0; i < cl_out->nat; i++)
           if( cl_out->srt[i] == cl->srt[at] &&
               ddvec(cl_out->r[cl_out->nat],cl_out->r[i]) < TOL
           ) goto FOUND;
	 if( cl_out->nat == MAX_ATOMS ) {
	    fprintf(stderr,"Parameter MAX_ATOMS=%d is too small! Increase it.\n", MAX_ATOMS);
	    exit(1);
	 }
	 cl_out->nat++;
	 FOUND: ;
      }
   }
   return;
}

/*****************************************************
 * find valid operations for space group;
 * op[][4][3] - point group operations are given;
 * cl - cell for one sort (minimal number of atoms)
 * clf - full cell (all atoms)
 *****************************************************/
void find_operations(int nop, int npgrp, int *nop_out, double op[][4][3],
                     int *ind, t_cell *cl, t_cell *clf)
{
    int nopp=0,i,j,io,k;
    double r[3],r1[3],t[3], op_out[48][4][3];

    for(io=0; io < nop; io++) {
        mul_mv(r,op[io],cl->r[0]);
        for(i=0; i < cl->nat; i++) { /*  probe all partial translations  */
            sub_vv(t,cl->r[i],r);
            reduce(t);
           /*  check for this partial translation  */
           for(j=0; j < clf->nat; j++) {
               mul_mv(r1,op[io],clf->r[j]);
               add_vv(r1,r1,t);
               reduce(r1);
              /*  find this atom  */
              for(k=0; k < clf->nat; k++)
              if( clf->srt[j] ==  clf->srt[k] &&
                  ddvec(r1,clf->r[k]) < TOL ) goto ATOM_FOUND;
              goto SKIP_TRANSLATION;
              ATOM_FOUND:;
           }
           asgn(op_out[nopp][0],op[io][0]);
           asgn(op_out[nopp][1],op[io][1]);
           asgn(op_out[nopp][2],op[io][2]);
           asgn(op_out[nopp][3],t);
           ind[nopp]=ind_pgrp[npgrp][io];
           nopp++;
           break; /*  next operation  */
      SKIP_TRANSLATION:;
        }
    }
    *nop_out=nopp;
    for(io=0;io < nopp; io++)
    for(i=0; i<4; i++)
    for(j=0; j<3; j++) op[io][i][j]=op_out[io][i][j];
   
    return;
}

/***************************************************
 * transform group operations into new basis
 * Uprim=(A^T)^-1*U*(A^T), rprim=r*A_1
 * *************************************************/
void transform_group_opr(int nop, double op[][4][3], double A[3][3])
{
   double A_1[3][3];
   int i;
   
   inv_matr(A,A_1);
   for(i=0; i<nop; i++) {
      mul_mm(op[i],A_1,op[i]);
      mul_mm(op[i],op[i],A);
   }
}

int is_inversion(int nop,double sym_op[][4][3],double Rsh[3])
{
   int i;
   for(i=0; i<nop; i++)
   if( fabs(sym_op[i][0][0]+1.) < TOL &&  fabs(sym_op[i][1][1]+1.) < TOL &&
       fabs(sym_op[i][2][2]+1.) < TOL &&
       fabs(sym_op[i][0][1]) < TOL &&  fabs(sym_op[i][0][2]) < TOL &&
       fabs(sym_op[i][1][0]) < TOL &&  fabs(sym_op[i][1][2]) < TOL &&
       fabs(sym_op[i][2][0]) < TOL &&  fabs(sym_op[i][2][1]) < TOL ) {
      /*  Rsh=-((Inv-ONE)^-1)*Rsh  */
       Rsh[0]=0.5*sym_op[i][3][0];
       Rsh[1]=0.5*sym_op[i][3][1];
       Rsh[2]=0.5*sym_op[i][3][2];
       return(1);
   }
   return(0);
}
   
void rotate_group_transl(int nop, double op[][4][3], double A[3][3])
{
   double A_1[3][3];
   int i;

   inv_matr(A,A_1);
   for(i=0; i<nop; i++) {
      mul_mv(op[i][3],A_1,op[i][3]);
      reduce(op[i][3]);
   }
}

/** shift origin on R vector **/
void shift_group_transl(int nop, double op[][4][3],double R[3])
{
   int i;
   double A[3][3],r1[3];
/*   r'=r+(op-E)*R  */
   for(i=0; i<nop; i++) {
      asgn_n(A[0],op[i][0],9);
      A[0][0]-=1.; A[1][1]-=1.; A[2][2]-=1.; /*  A=op-E  */
      mul_mv(r1,A,R);
      add_vv(op[i][3],op[i][3],r1);
      reduce(op[i][3]);
   }
}

void holohedral_pgrp(int lat, int *npgrp, double sym_op[][4][3])
{
   int i;
   if( lat==CUBIC_P || lat==CUBIC_I || lat==CUBIC_F ) {
      *npgrp=NOh;
      for(i=0;i<nop_pgrp[NOh];i++) asgn_n(sym_op[i][0],Oh_pgrp[i][0],9);
      return;
   }
   if( lat==HEXAGONAL ) {
      *npgrp=ND6h;
      for(i=0;i<nop_pgrp[ND6h];i++) asgn_n(sym_op[i][0],D6h_pgrp[i][0],9);
      return;
   }
   if( lat==RHOMBOHEDRAL ) {
      *npgrp=ND3d2;
      for(i=0;i<nop_pgrp[ND3d2];i++) asgn_n(sym_op[i][0],D3d2_pgrp[i][0],9);
      return;
   }
   if( lat==TETRAGONAL_P || lat==TETRAGONAL_I) {
      *npgrp=ND4h;
      for(i=0;i<nop_pgrp[ND4h];i++) asgn_n(sym_op[i][0],D4h_pgrp[i][0],9);
      return;
   }
   if( lat==ORTHOROMBIC_P || lat==ORTHOROMBIC_I ||
       lat==ORTHOROMBIC_C || lat==ORTHOROMBIC_A ||
       lat==ORTHOROMBIC_F ) {
      *npgrp=ND2h;
      for(i=0;i<nop_pgrp[ND2h];i++) asgn_n(sym_op[i][0],D2h_pgrp[i][0],9);
      return;
   }
   if( lat==MONOCLINIC_P || lat==MONOCLINIC_A ) {
      *npgrp=NC2h;
      for(i=0;i<nop_pgrp[NC2h];i++) asgn_n(sym_op[i][0],C2h_pgrp[i][0],9);
      return;
   }
   if( lat==TRICLINIC ) {
      *npgrp=NCi;
      for(i=0;i<nop_pgrp[NCi];i++) asgn_n(sym_op[i][0],Ci_pgrp[i][0],9);
      return;
   }
}

/* decompositions of primitive vectors over vectors of elem. cell */
/* 1 vector is 1 column, ... */
int get_Telpr(int lat, double Telpr[3][3])
{
   switch( lat ) {
    case CUBIC_P:
    case HEXAGONAL:
    case TETRAGONAL_P:
    case ORTHOROMBIC_P:
    case MONOCLINIC_P:
    case TRICLINIC:
                    Telpr[0][0]= 1.; Telpr[0][1]= 0.; Telpr[0][2]= 0.;
                    Telpr[1][0]= 0.; Telpr[1][1]= 1.; Telpr[1][2]= 0.;
                    Telpr[2][0]= 0.; Telpr[2][1]= 0.; Telpr[2][2]= 1.;
                    break;
    case CUBIC_I:
    case TETRAGONAL_I:
    case ORTHOROMBIC_I:
                    Telpr[0][0]=-0.5;  Telpr[0][1]= 0.5; Telpr[0][2]= 0.5;
                    Telpr[1][0]= 0.5;  Telpr[1][1]=-0.5; Telpr[1][2]= 0.5;
                    Telpr[2][0]= 0.5;  Telpr[2][1]= 0.5; Telpr[2][2]=-0.5;
                    break;
    case CUBIC_F:
    case ORTHOROMBIC_F:
                    Telpr[0][0]= 0.0; Telpr[0][1]= 0.5; Telpr[0][2]= 0.5;
                    Telpr[1][0]= 0.5; Telpr[1][1]= 0.0; Telpr[1][2]= 0.5;
                    Telpr[2][0]= 0.5; Telpr[2][1]= 0.5; Telpr[2][2]= 0.0;
                    break;
    case ORTHOROMBIC_C:
                    Telpr[0][0]= 0.5; Telpr[0][1]= 0.5; Telpr[0][2]= 0.0;
                    Telpr[1][0]=-0.5; Telpr[1][1]= 0.5; Telpr[1][2]= 0.0;
                    Telpr[2][0]= 0.0; Telpr[2][1]= 0.0; Telpr[2][2]= 1.0;
                    break;
    case ORTHOROMBIC_A:
    case MONOCLINIC_A:
                    Telpr[0][0]= 1.0; Telpr[0][1]= 0.0; Telpr[0][2]= 0.0;
                    Telpr[1][0]= 0.0; Telpr[1][1]= 0.5; Telpr[1][2]=-0.5;
                    Telpr[2][0]= 0.0; Telpr[2][1]= 0.5; Telpr[2][2]= 0.5;
                    break;
    case RHOMBOHEDRAL:
                    Telpr[0][0]=-1./3; Telpr[0][1]= 2./3; Telpr[0][2]=-1./3;
                    Telpr[1][0]=-2./3; Telpr[1][1]= 1./3; Telpr[1][2]= 1./3;
                    Telpr[2][0]= 1./3; Telpr[2][1]= 1./3; Telpr[2][2]= 1./3;
                    break;
/*   not implemented  */
    case ORTHOROMBIC_B:
                    Telpr[0][0]= 0.5; Telpr[0][1]= 0.0; Telpr[0][2]= 0.5;
                    Telpr[1][0]= 0.0; Telpr[1][1]= 1.0; Telpr[1][2]= 0.0;
                    Telpr[2][0]=-0.5; Telpr[2][1]= 0.0; Telpr[2][2]= 0.5;
                    break;
      
    default:        return(0);
   }
   return(1);
}
   
int find_space_group(int lat, int nop, int npgrp, int *nshft,
                     double sym_op[][4][3], char **sgrp_name,
		     double Rshft[3], double Rsh[][3], double Tnew[3][3],
		     int *fix_org, int is_prim)
{
   char **sname;
   int i,i1,i2,j1,j2,j3,k1,k2,k3,j,isgrp,ib,num_sgrp,nsh;
   int *lat_ptr,Nb=24;
   double (*r_ptr)[3], (*A)[3][3];
   double op[48][4][3],r[48][3],rop[48][3],vec[4][3],ge1[3][3];
   double A_1[3][3],Telpr[3][3],Tprel[3][3];

   switch( lat ) {
    case CUBIC_P: A=Rot_cub_p; Nb=2; break;
    case CUBIC_I: A=Rot_one; Nb=1; break;
    case CUBIC_F: A=Rot_one; Nb=1; break;

    case HEXAGONAL: A=Rot_one; Nb=1; break;
    case RHOMBOHEDRAL:  A=Rot_one; Nb=1; break;
      
    case TETRAGONAL_P: A=Rot_one; Nb=1; break;
    case TETRAGONAL_I: A=Rot_one; Nb=1; break;
      
    case ORTHOROMBIC_P: A=Rot_ort_p; Nb=6; break;
    case ORTHOROMBIC_I: A=Rot_ort_i; Nb=4; break;
    case ORTHOROMBIC_C: A=Rot_ort_c; Nb=2; break;
    case ORTHOROMBIC_A: A=Rot_one; Nb=1; break;
    case ORTHOROMBIC_F: A=Rot_one; Nb=1; break;
      
    case MONOCLINIC_P: A=Rot_mon_p; Nb=3; break;
    case MONOCLINIC_A: A=Rot_one; Nb=1; break;

    case TRICLINIC: A=Rot_one; Nb=1; break;
   }
   get_Telpr(lat,Telpr);  inv_matr(Telpr,Tprel);

   /*** TRICLINIC ***/
   if( npgrp==NC1 ) {
      r_ptr=rC1_sgrp[0]; lat_ptr=lat_C1_sgrp; num_sgrp=NC1_sgrp; sname=comnt_C1_sgrp;
      i1=0; i2=0; *fix_org=111;/*  x,y,z not fixed  */
   } else
   if( npgrp==NCi ) {
      r_ptr=rCi_sgrp[0]; lat_ptr=lat_Ci_sgrp; num_sgrp=NCi_sgrp; sname=comnt_Ci_sgrp;
      i1=1; i2=1; *fix_org=000;
   } else
   /*** MONOCLIINIC ***/
   if( npgrp==NC2 ) {
      r_ptr=rC2_sgrp[0]; lat_ptr=lat_C2_sgrp; num_sgrp=NC2_sgrp; sname=comnt_C2_sgrp;
      i1=1; i2=0; *fix_org=001;/*  z not fixed (shift (0,0,z) with arbitrary z value is allowed)  */
   } else
   if( npgrp==NCs ) {
      r_ptr=rCs_sgrp[0]; lat_ptr=lat_Cs_sgrp; num_sgrp=NCs_sgrp; sname=comnt_Cs_sgrp;
      i1=0; i2=1; *fix_org=110; /*   x,y not fixed  */
   } else
   if( npgrp==NC2h ) {
      r_ptr=rC2h_sgrp[0]; lat_ptr=lat_C2h_sgrp; num_sgrp=NC2h_sgrp; sname=comnt_C2h_sgrp;
      i1=2; i2=2; *fix_org=000;
   } else
 /*** ORTHOROMBIC ***/
   if( npgrp==ND2 ) {
      r_ptr=rD2_sgrp[0]; lat_ptr=lat_D2_sgrp; num_sgrp=ND2_sgrp; sname=comnt_D2_sgrp;
      i1=1; i2=2; *fix_org=0;
   } else
   if( npgrp==NC2v ) {
      r_ptr=rC2v_sgrp[0]; lat_ptr=lat_C2v_sgrp; num_sgrp=NC2v_sgrp; sname=comnt_C2v_sgrp;
      i1=1; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==ND2h ) {
      r_ptr=rD2h_sgrp[0]; lat_ptr=lat_D2h_sgrp; num_sgrp=ND2h_sgrp; sname=comnt_D2h_sgrp;
      i1=4; i2=4; *fix_org=0;
   } else
 /*** TETRAGONAL ***/
   if( npgrp==NC4 ) {
      r_ptr=rC4_sgrp[0]; lat_ptr=lat_C4_sgrp; num_sgrp=NC4_sgrp; sname=comnt_C4_sgrp;
      i1=1; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==NS4 ) {
      r_ptr=rS4_sgrp[0]; lat_ptr=lat_S4_sgrp; num_sgrp=NS4_sgrp; sname=comnt_S4_sgrp;
      i1=2; i2=2; *fix_org=0;
   } else
   if( npgrp==NC4h ) {
      r_ptr=rC4h_sgrp[0]; lat_ptr=lat_C4h_sgrp; num_sgrp=NC4h_sgrp; sname=comnt_C4h_sgrp;
      i1=4; i2=4; *fix_org=0;
   } else
   if( npgrp==ND4 ) {
      r_ptr=rD4_sgrp[0]; lat_ptr=lat_D4_sgrp; num_sgrp=ND4_sgrp; sname=comnt_D4_sgrp;
      i1=1; i2=4; *fix_org=0;
   } else
   if( npgrp==NC4v ) {
      r_ptr=rC4v_sgrp[0]; lat_ptr=lat_C4v_sgrp; num_sgrp=NC4v_sgrp; sname=comnt_C4v_sgrp;
      i1=1; i2=0; *fix_org=001;
   } else
   if( npgrp==ND2d ) {
      r_ptr=rD2d_sgrp[0]; lat_ptr=lat_D2d_sgrp; num_sgrp=ND2d_sgrp; sname=comnt_D2d_sgrp;
      i1=1; i2=4; *fix_org=0;
   } else
   if( npgrp==ND2d2 ) {
      r_ptr=rD2d2_sgrp[0]; lat_ptr=lat_D2d2_sgrp; num_sgrp=ND2d2_sgrp; sname=comnt_D2d2_sgrp;
      i1=2; i2=2; *fix_org=0;
   } else
   if( npgrp==ND4h ) {
      r_ptr=rD4h_sgrp[0]; lat_ptr=lat_D4h_sgrp; num_sgrp=ND4h_sgrp; sname=comnt_D4h_sgrp;
      i1=8; i2=8; *fix_org=0;
   } else
 /*** TRIGONAL ***/
   if( npgrp==NC3 ) {
      r_ptr=rC3_sgrp[0]; lat_ptr=lat_C3_sgrp; num_sgrp=NC3_sgrp; sname=comnt_C3_sgrp;
      i1=1; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==NC3i ) {
      r_ptr=rC3i_sgrp[0]; lat_ptr=lat_C3i_sgrp; num_sgrp=NC3i_sgrp; sname=comnt_C3i_sgrp;
      i1=3; i2=3; *fix_org=0;
   }
   if( npgrp==ND3 ) {
      r_ptr=rD3_sgrp[0]; lat_ptr=lat_D3_sgrp; num_sgrp=ND3_sgrp; sname=comnt_D3_sgrp;
      i1=1; i2=3; *fix_org=0;
   } else
   if( npgrp==ND32 ) {
      r_ptr=rD32_sgrp[0]; lat_ptr=lat_D32_sgrp; num_sgrp=ND32_sgrp; sname=comnt_D32_sgrp;
      i1=1; i2=3; *fix_org=0;
   } else
   if( npgrp==NC3v ) {
      r_ptr=rC3v_sgrp[0]; lat_ptr=lat_C3v_sgrp; num_sgrp=NC3v_sgrp; sname=comnt_C3v_sgrp;
      i1=1; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==NC3v2 ) {
      r_ptr=rC3v2_sgrp[0]; lat_ptr=lat_C3v2_sgrp; num_sgrp=NC3v2_sgrp; sname=comnt_C3v2_sgrp;
      i1=1; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==ND3d ) {
      r_ptr=rD3d_sgrp[0]; lat_ptr=lat_D3d_sgrp; num_sgrp=ND3d_sgrp; sname=comnt_D3d_sgrp;
      i1=6; i2=6; *fix_org=0;
   } else
   if( npgrp==ND3d2 ) {
      r_ptr=rD3d2_sgrp[0]; lat_ptr=lat_D3d2_sgrp; num_sgrp=ND3d2_sgrp; sname=comnt_D3d2_sgrp;
      i1=6; i2=6; *fix_org=0;
   } else
 /*** HEXAGONAL ***/
   if( npgrp==NC6 ) {
      r_ptr=rC6_sgrp[0]; lat_ptr=lat_C6_sgrp; num_sgrp=NC6_sgrp; sname=comnt_C6_sgrp;
      i1=4; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==NC3h ) {
      r_ptr=rC3h_sgrp[0]; lat_ptr=lat_C3h_sgrp; num_sgrp=NC3h_sgrp; sname=comnt_C3h_sgrp;
      i1=4; i2=4; *fix_org=0;
   } else
   if( npgrp==NC6h ) {
      r_ptr=rC6h_sgrp[0]; lat_ptr=lat_C6h_sgrp; num_sgrp=NC6h_sgrp; sname=comnt_C6h_sgrp;
      i1=10; i2=10; *fix_org=0;
   } else
   if( npgrp==ND6 ) {
      r_ptr=rD6_sgrp[0]; lat_ptr=lat_D6_sgrp; num_sgrp=ND6_sgrp; sname=comnt_D6_sgrp;
      i1=3; i2=6; *fix_org=0; /*  C2(001) C2(110)  */
   } else
   if( npgrp==NC6v ) {
      r_ptr=rC6v_sgrp[0]; lat_ptr=lat_C6v_sgrp; num_sgrp=NC6v_sgrp; sname=comnt_C6v_sgrp;
      i1=4; i2=0; *fix_org=001; /*  z direction not fixed  */
   } else
   if( npgrp==ND3h ) {
      r_ptr=rD3h_sgrp[0]; lat_ptr=lat_D3h_sgrp; num_sgrp=ND3h_sgrp; sname=comnt_D3h_sgrp;
      i1=4; i2=4; *fix_org=0;
   } else
   if( npgrp==ND3h2 ) {
      r_ptr=rD3h2_sgrp[0]; lat_ptr=lat_D3h2_sgrp; num_sgrp=ND3h2_sgrp; sname=comnt_D3h2_sgrp;
      i1=4; i2=4; *fix_org=0;
   } else
   if( npgrp==ND6h ) {
      r_ptr=rD6h_sgrp[0]; lat_ptr=lat_D6h_sgrp; num_sgrp=ND6h_sgrp; sname=comnt_D6h_sgrp;
      i1=16; i2=16; *fix_org=0;
   } else
 /*** CUBIC ***/
   if( npgrp==NT )  {
      r_ptr=rT_sgrp[0];  lat_ptr=lat_T_sgrp;  num_sgrp=NT_sgrp; sname=comnt_T_sgrp;
      i1=1; i2=2; *fix_org=0;
   } else
   if( npgrp==NTh ) {
      r_ptr=rTh_sgrp[0]; lat_ptr=lat_Th_sgrp; num_sgrp=NTh_sgrp; sname=comnt_Th_sgrp;
      i1=1; i2=2; *fix_org=0;
   } else
   if( npgrp==NO )  {
      r_ptr=rO_sgrp[0];  lat_ptr=lat_O_sgrp;  num_sgrp=NO_sgrp; sname=comnt_O_sgrp;
      i1=1; i2=2; *fix_org=0;
   } else
   if( npgrp==NTd ) {
      r_ptr=rTd_sgrp[0]; lat_ptr=lat_Td_sgrp; num_sgrp=NTd_sgrp; sname=comnt_Td_sgrp;
      i1=1; i2=2; *fix_org=0;
   } else
   if( npgrp==NOh ) {
      r_ptr=rOh_sgrp[0]; lat_ptr=lat_Oh_sgrp; num_sgrp=NOh_sgrp; sname=comnt_Oh_sgrp;
      i1=1; i2=2; *fix_org=0;
/*        i1=24; i2=24; *fix_org=0;  */
   };
   
   /*  construct (g-E)^-1  */
   /*  transform group operations to prim. basis  */
   for(i=0; i<nop; i++) {
       mul_mm(sym_op[i],sym_op[i],Telpr);
       mul_mm(sym_op[i],Tprel,sym_op[i]);
       mul_mv(sym_op[i][3],Tprel,sym_op[i][3]);  /*  change translation part  */
       reduce(sym_op[i][3]);
   }
   for( ib=0; ib<Nb; ib++ ) {
        mul_mm(A[ib],A[ib],Telpr);
        mul_mm(A[ib],Tprel,A[ib]);
   }

   for(i=0;i<9;i++) ge1[0][i]=0.;
   asgn_n(op[0][0],sym_op[i1][0],9); op[0][0][0]--; op[0][1][1]--; op[0][2][2]--;
   asgn_n(op[1][0],sym_op[i2][0],9); op[1][0][0]--; op[1][1][1]--; op[1][2][2]--;
   
   k1=0; k2=1; k3=2;
   if( *fix_org == 110 ) { /*   Cs group only  */
      if( lat==MONOCLINIC_P) ge1[2][2]=-0.5;
      else ge1[2][2]=-1.0;
   } else
   if( *fix_org == 001 ) {
     for(k1=0; k1<3; k1++)
     for(k2=k1+1; k2<3; k2++) {
        vec[0][0]=op[0][k1][k1]*op[0][k2][k2]-op[0][k1][k2]*op[0][k2][k1];
        if( fabs(vec[0][0])>TOL ) {
	   ge1[k1][k1]= op[0][k2][k2]/vec[0][0];
	   ge1[k2][k2]= op[0][k1][k1]/vec[0][0];
	   ge1[k1][k2]=-op[0][k1][k2]/vec[0][0];
	   ge1[k2][k1]=-op[0][k2][k1]/vec[0][0];
           if( (k3=(k1+1)%3) == k2 ) k3=(k1+2)%3;
	   goto GE1;
	}
     }
     fprintf(stderr,"Error 101 in find_space_group().\n");
     exit(1);
   } else 	
   if( *fix_org == 0 ) {
     for(k1=0; k1<3; k1++)
     for(k2=k1+1; k2<3; k2++) {
        if( (k3=(k1+1)%3) == k2 ) k3=(k1+2)%3;
	asgn(ge1[k1],op[0][k1]);  asgn(ge1[k2],op[0][k2]);  asgn(ge1[k3],op[1][k3]);
        if( fabs(detr(ge1)) > TOL ) {
	   inv_matr(ge1,ge1);
	   goto GE1;
	}
     }
     fprintf(stderr,"Error 102 in find_space_group().\n");
     exit(1);
   }
   GE1:

   /* find now all shift vectors "r" which hold relation: (g-E)r=R
      where R - is vector with integer components,
      g - operation given by indexes i1,i2 */
   if( lat==RHOMBOHEDRAL || lat==HEXAGONAL ) {
      rop[0][0]=0.; rop[0][1]=1./2; rop[0][2]=1./3; rop[0][3]=2./3;
   } else {
      rop[0][0]=0.; rop[0][1]=1./4; rop[0][2]=1./2; rop[0][3]=3./4;
   }
   nsh=0;
   asgn(vec[k1],op[0][k1]); asgn(vec[k2],op[0][k2]); asgn(vec[k3],op[1][k3]);
   for(j1=0;j1<=3;j1++)
   for(j2=0;j2<=3;j2++)
   for(j3=0;j3<=3;j3++) {
      Rsh[nsh][0]=rop[0][j1];
      Rsh[nsh][1]=rop[0][j2];
      Rsh[nsh][2]=rop[0][j3];
      mul_mv(Rsh[nsh],Telpr,Rsh[nsh]);
      if( *fix_org==111 ) Rsh[nsh][0]=Rsh[nsh][1]=Rsh[nsh][2]=0; else
      if( *fix_org==110 ) Rsh[nsh][0]=Rsh[nsh][1]=0; else
      if( *fix_org==001 ) Rsh[nsh][2]=0;
      mul_mv(Rsh[nsh],Tprel,Rsh[nsh]);
      reduce(Rsh[nsh]);
      mul_mv(rop[2],vec,Rsh[nsh]); reduce(rop[2]);
      if( dvec(rop[2]) > TOL ) goto NEXT_SHIFT1;
      for(i=0; i<nsh; i++)
	if( ddvec(Rsh[nsh],Rsh[i]) < TOL ) goto NEXT_SHIFT1; /*  already exist  */
#if DEBUG==1
    printf("% .4f % .4f % .4f\n ",Rsh[nsh][0],Rsh[nsh][1],Rsh[nsh][2]);
#endif
      
      nsh++;
      if( nsh > 8 ) {
         fprintf(stderr,"Error103 in find_space_group()!\n");	
	 exit(1);
      }
      NEXT_SHIFT1:;
   }
/*     printf("%d ",nsh);  */
   for( isgrp=0; isgrp < num_sgrp; isgrp++ ) {
      if( lat != lat_ptr[isgrp] ) {
          r_ptr+=nop;
          continue;
      }
      asgn_n(r[0],&(*r_ptr)[0],nop*3);
      for(i=0; i<nop; i++) { mul_mv(r[i],Tprel,r[i]); reduce(r[i]); }
      for( ib=0; ib<Nb; ib++ ) { /*  probe all orientation of basis  */
        inv_matr(A[ib],A_1);
        asgn_n(op[0][0],sym_op[0][0],nop*4*3);
	/*  probe possible orientations of the elementary basis  */
        for(i=0; i<nop; i++) { /*  change point operations  */
          mul_mm(op[i],A_1,op[i]);
          mul_mm(op[i],op[i],A[ib]);
          mul_mv(op[i][3],A_1,op[i][3]);  /*  change translation part  */
        }
       /*  sort group  */
        for(i=0; i<nop; i++) {
          for(j=i; j<nop; j++)
          if( ddvec_n(sym_op[i][0],op[j][0],9) < TOL ) { /*  swap op[i] op[j]  */
              asgn_n(vec[0],op[i][0],4*3);
              asgn_n(op[i][0],op[j][0],4*3);
              asgn_n(op[j][0],vec[0],4*3);
              goto FOUND;
          }
	  goto NEXT_BASIS;
          FOUND:;
        }
     /*  shift coordinate system; r' = r + (g-E)*Rsh; Rsh= (g-E)^(-1)*(r'-r)  */
	vec[3][k1]=r[i1][k1]-op[i1][3][k1];
	vec[3][k2]=r[i1][k2]-op[i1][3][k2];
	vec[3][k3]=r[i2][k3]-op[i2][3][k3];
	mul_mv(vec[3],ge1,vec[3]);
/*   r'=r+(op-E)*Rsh  */
        for(j1=0; j1<nsh; j1++) {
#if DEBUG==1
    printf("\n--------------------------------------------\n");
#endif
 	  for(i=0; i<nop; i++) {
            asgn_n(vec[0],op[i][0],9);
            vec[0][0]--; vec[1][1]--; vec[2][2]--; /*  A=op-E  */
	    add_vv(rop[i],vec[3],Rsh[j1]);
	    mul_mv(rop[i],vec,rop[i]);
            add_vv(rop[i],op[i][3],rop[i]);
            reduce(rop[i]);
#if DEBUG==1
    printf("(% .4f % .4f % .4f) ",rop[i][0],rop[i][1],rop[i][2]);
#endif
            if( ddvec(rop[i],r[i]) > TOL ) goto NEXT_SHIFT;	
	  }
     /*   found!  */
          *sgrp_name=sname[isgrp];
	  asgn_n(Tnew[0],A[ib][0],9);
	  mul_mm(Tnew,Tnew,Tprel); mul_mm(Tnew,Telpr,Tnew);
	  add_vv(Rshft,vec[3],Rsh[j1]);
	  if( !is_prim ) mul_mv(Rshft,Telpr,Rshft);
	 /*  find possible shift vectors  */
	  for(i=0; i<nsh; i++) asgn(rop[i],Rsh[i]);
	  j=nsh; nsh=0;
	  for(i=0; i<j; i++) {
	     for(j2=0; j2<nop; j2++) {
		asgn_n(vec[0],sym_op[j2][0],9); vec[0][0]--; vec[1][1]--; vec[2][2]--;
		mul_mv(Rsh[nsh],vec,rop[i]);
		reduce(Rsh[nsh]);
		if( dvec(Rsh[nsh]) > TOL ) goto NEXT_SHIFT2; /*  next shift vector  */
	     }
	     asgn(Rsh[nsh],rop[i]);
	     if( !is_prim ) mul_mv(Rsh[nsh],Telpr,Rsh[nsh]);
	     nsh++;
	     NEXT_SHIFT2:;
	  }
	  *nshft=nsh;	
	  for(i=0; i<nop; i++) asgn(sym_op[i][3],*(r_ptr+i));
	  if( is_prim ) /*  transform to prim. basis  */
	  for(i=0; i<nop; i++) {
	      mul_mv(sym_op[i][3],Tprel,sym_op[i][3]);  /*  change translation part  */
	      reduce(sym_op[i][3]);
	  }
	  if( !is_prim )  /*  transform to elem. basis  */
          for(i=0; i<nop; i++) {
              mul_mm(sym_op[i],sym_op[i],Tprel);
              mul_mm(sym_op[i],Telpr,sym_op[i]);
	  }
      /*  transform shifts to elem. basis  */
/*	   for(i=0; i<nsh; i++) {
	      mul_mv(Rsh[i],Telpr,Rsh[i]);
	      if( *fix_org==111 ) Rsh[i][0]=Rsh[i][1]=Rsh[i][2]=0; else
	      if( *fix_org==110 ) Rsh[i][0]=Rsh[i][1]=0; else
	      if( *fix_org==001 ) Rsh[i][2]=0;
	      reduce(Rsh[i]);
	   }*/
          return(0);
          NEXT_SHIFT:;
	}
	NEXT_BASIS:;
      } /*   ib=0  */
      r_ptr+=nop;
   } /*   isgrp  */
   fprintf(stderr,"Error2 in find_space_group()!\n");
   return(1);
}

void split_cell(int is_prim, t_cell *cl_in, t_cell *cl_new, int lat,
		int nop, double op[][4][3],
		char name_srt[MAX_ATOMS][MAX_CHARS],char name_srt_new[MAX_ATOMS][MAX_CHARS])
{
   int i1,i2,j1,j2,j,iop,n,nsub_srt;
   double Telpr[3][3],Tprel[3][3],r[3],rj[3];
   
   /*  sort cell_in  */
   for( j=0; j<cl_in->nat; j++ )
   for( n=j+1; n<cl_in->nat; n++ )
      if( cl_in->srt[n] < cl_in->srt[j] ) {
          swap(cl_in->r[j],cl_in->r[n]);
          i1=cl_in->srt[n]; cl_in->srt[n]=cl_in->srt[j]; cl_in->srt[j]=i1;
      }
   if( is_prim )
     get_Telpr(CUBIC_P,Telpr); /*  if cl_new is given in prim. basis Telpr=unity matrix  */
   else
     get_Telpr(lat,Telpr);
   inv_matr(Telpr,Tprel);
   
   cl_new->nat=cl_in->nat; cl_new->nsrt=0;
   for(j=0; j<cl_in->nat; j++ ) asgn(cl_new->r[j],cl_in->r[j]);
   /*  i1 - left boundary; i2 - right boundary  */
   i1=i2=0;
   for(n=0; n<cl_in->nsrt; n++) {
      while( (i2+1)<cl_in->nat && cl_in->srt[i2+1]==n ) i2++;
      nsub_srt=0; /*  find all subsorts  */
      for(j1=i1; j1<=i2; j1=j2+1 ) {
	 j2=j1;
	 cl_new->srt[j1]=cl_new->nsrt;
         for(iop=0; iop<nop; iop++ ) {	
           mul_mv(r,op[iop],cl_new->r[j1]);
           add_vv(r,r,op[iop][3]);
	   mul_mv(r,Tprel,r); reduce(r);
           for(j=j2+1; j<=i2; j++) {
      	      mul_mv(rj,Tprel,cl_new->r[j]); reduce(rj);
              if( ddvec(rj,r)<TOL ) {
	        j2++;
  	        swap(cl_new->r[j2],cl_new->r[j]);
 	        cl_new->srt[j2]=cl_new->nsrt;
	    }
	   }
	 }
         sprintf(name_srt_new[cl_new->nsrt],"%s (%i)",name_srt[n],nsub_srt+1);
	 if( cl_new->nsrt == MAX_ATOMS ) {
	    fprintf(stderr,"Parameter MAX_ATOMS=%d too small! Increase it.\n", MAX_ATOMS);
	    exit(1);
	 }
	 cl_new->nsrt++;
         nsub_srt++;
      }
      i1=i2+1; i2=i1;
   }
} /*  end of  split_cell()  */

/********************************************************************
 * find point group for each atom.  shift origin to atom at position Rsh,
 * and select only operations with zero translation part
 * Rnew= (g-E)Rsh + Rold;  Rnew=0;
 * nops - number of operations for space group
 *********************************************************************/
void find_atom_pgroup( int lat, int nop_s, int npgrp_s, double sym_op[][4][3], double Rsh[3],
		       int *npgrp_p, double op[][4][3], /*  number and operations of point group  */
		       double T[3][3]) /*  orientation of new elementary basis  */
{
   int i, lat_new, nop_p=0, ind[48];
   double opE[4][3];
   
   get_Telpr(lat,T);  inv_matr(T,T);
   for(i=0; i<nop_s; i++) {
       asgn_n(opE[0],sym_op[i][0],9);
       opE[0][0]--; opE[1][1]--; opE[2][2]--; /*  A=op-E  */
       mul_mv(opE[3],opE,Rsh);
       add_vv(opE[3],opE[3],sym_op[i][3]);
       mul_mv(opE[3],T,opE[3]);
       reduce(opE[3]);
       if( dvec(opE[3]) < TOL ) {
	  asgn_n(op[nop_p][0],sym_op[i][0],12);
          ind[nop_p++]=ind_pgrp[npgrp_s][i];
       }
   }
/*  for triclinic lattice axes preserve old orentaton; here Ttrcl=opE  */
   for(i=0; i<9; i++) opE[0][i]=0.; opE[0][0]=opE[1][1]=opE[2][2]=1.;
  /*  determine new orientation of basis for this point group  */
   new_transl(lat,nop_p,ind,opE,&lat_new,T);
  /*  transform group opr. from old elementary basis into new one  */
   transform_group_opr(nop_p,op,T);
   *npgrp_p=det_pgrp(nop_p,op);
}

double     a[2][3], a_pr[2][3],
           Rsh[9][3], /*  possible shifts for space group  */
           Tin[3][3],Tin_old[3][3],        /*  input basis vectors  */
           Tp[3][3],                       /*  primitive basis vectors; Tp_1=Tp^-1  */
           Tel[3][3],                      /*  nonprimitive basis vectors  */
           Telpr[3][3],Tprel[3][3], /*  transf. matrix from el to prim. basis=(Tel^-1)*Tp  */
           Ttrcl[3][3];    /*  for triclinic lattice only, chtobi ne mudoxat'sia s tricl.  */
                           /*  bazisom, on naxoditsia v basis_transl odin raz i vse  */
           
t_cell     cl_in, cl_new, cl_pr;
int        lat, lat_new, nop, npgrp, npgrp_at, ind[48], fix_org, nsh,
           is_noeq=0, /* print only atoms of primitive cell in elem. basis*/
           is_prim=0, /* print atoms of pritive cell in prim. basis */
           is_wien_input=0, /* read data from WIEN input file*/
           is_wien_rhomb=0; /* for RHOMBOHEDRAL WIEN uses the hexagonal lattice         */
			    /* parameters but coord. of atoms are given in RHOMB. basis */
double     sym_op[48][4][3], /*  operations of space group  */
           pgrp_at[48][4][3]; /*  point group for atom; determined for each atom and printed  */
char       name_srt[MAX_ATOMS][MAX_CHARS],name_srt_new[MAX_ATOMS][MAX_CHARS];
char      *file_in=NULL, *file_out = NULL;
int        min_srt;
char      *sgrp_name;
 
int main(int argc, char *argv[])
{
   int i,j;
   /*  parse command line  */
   j=0;
   for(i=1; i<argc; i++) {
      if( strcmp("-help",argv[i])==0 ) {
	 printf("%s",help_str);
	 return(0);
      }	else
      if( strcmp("-prim",argv[i])==0 ) is_prim=1;       else
      if( strcmp("-noeq",argv[i])==0 ) is_noeq=1;       else
      if( strcmp("-wien",argv[i])==0 ) is_wien_input=1; else
      if( file_in  == NULL ) file_in=argv[i];           else
      if( file_out == NULL ) file_out=argv[i];
   }
   if( file_in == NULL ) {
      printf("No input file specified.\n");
      return(1);
   }
   if( is_prim && is_noeq ) is_noeq=0;
      
/******************************************************************************
 * Tin - decompos. of the input  basis translation vectors over cubic basis
 * basis_transl() Tel - find out  transl. vector (elem. cell)
 * new_transl()   Tel - find out  transl. vector consistent with orientation of
 *                      point group
 * find_space_group() - find out space group and transl. vectors
 ******************************************************************************/
/* section of geting data */
   if( is_wien_input ) {
       if( read_wien_data(file_in,&cl_in,name_srt,a,&lat,&is_wien_rhomb) ) exit(1);
   }
   else {
       if( read_data(file_in,&cl_in,name_srt,a,&lat) ) exit(1);
   }
   T_matr(Tin,a);
   get_Telpr(lat,Telpr); inv_matr(Telpr,Tprel);
   mul_mm(Tin_old,Tin,Tprel); /*  decompos. of input translations over cubic basis  */
   transform_coor(&cl_in,Tprel);
   
   for(i=0; i < cl_in.nat; i++) reduce(cl_in.r[i]);
   for(i=0; i < cl_in.nat; i++)
   for(j=i+1; j<cl_in.nat; j++)
     if( ddvec(cl_in.r[i],cl_in.r[j]) < TOL ) {
	fprintf(stderr,"Error: duplicated atoms found! Atoms #%d and #%d coincide.\n",i+1,j+1);
	exit(1);
     }
/* end of read section */
   
   min_srt=find_min_sort(Rsh[0],&cl_in);
/* atom of sort min_srt must occupy (0,0,0) position to find translations vectors*/
   shift_cell(Rsh[0],&cl_in);
   transform_coor(&cl_in,Tin); /*  transform into cubic basis  */
   make_cell_for_sort(min_srt,&cl_new,&cl_in); /*  cl_new here contains only atoms of single sort  */
   basis_transl(&cl_in,&cl_new,&lat,Tin,Tp,Tel,Ttrcl);

#if DEBUG==1
   inv_matr(Tin,Rsh);
   asgn_n(Rsh[3],Rsh[0],9);
   mul_mm(Rsh,Rsh,Tel);
   mul_mm(Rsh[3],Rsh[3],Tp);
#endif
   
   asgn_n(Tin[0],Tel[0],9);
   build_cell(&cl_in,&cl_new,Tp); /*  cl_in in cart. coord., cl_new in prim.  */
   make_cell_for_sort(min_srt,&cl_in,&cl_new);/*  cl_in here will contain only atoms of one sort  */
   inv_matr(Tel,Telpr);
   mul_mm(Ttrcl,Telpr,Ttrcl); /*  razl. tricl. bazisa po elem.(ne cubicheskomu)  */
   mul_mm(Telpr,Telpr,Tp); inv_matr(Telpr,Tprel);
   holohedral_pgrp(lat,&npgrp,sym_op); /*  npgrp = numb. of hol. point group  */
/*  transform group opr. from elementary into prim. basis  */
   transform_group_opr(nop_pgrp[npgrp],sym_op,Telpr);
   find_operations(nop_pgrp[npgrp],npgrp,&nop,sym_op,ind,&cl_in,&cl_new);
/*  transform from prim. to elementary  basis  */
   transform_group_opr(nop,sym_op,Tprel);
   
   /*  get new basis translation Tel for point group  */
   new_transl(lat,nop,ind,Ttrcl,&lat_new,Tel);
   mul_mm(Tprel,Tprel,Tel);
   mul_mm(Tin,Tin,Tel);
   if( lat_new==TRICLINIC || lat_new==MONOCLINIC_P || lat_new==ORTHOROMBIC_P ||
       lat_new==TETRAGONAL_P || lat_new==HEXAGONAL || lat_new==CUBIC_P )
       is_prim=is_noeq=0;
/*  transform group opr. from old elementary basis into new one  */
   transform_group_opr(nop,sym_op,Tel);
   rotate_group_transl(nop,sym_op,Tprel);
   npgrp=det_pgrp(nop,sym_op); /*  not only determines but also sorts sym_op in correct order  */
 
/*  find space group and new orientation of the basis  */
   find_space_group(lat_new,nop,npgrp,&nsh,sym_op,&sgrp_name,
		    Rsh[8],Rsh,Tel,&fix_org,is_prim);

   mul_mm(Tin,Tin,Tel);
   mul_mm(Tprel,Tprel,Tel); /*  decomp. of new elem. basis into old primitive  */
   if( !is_prim ) {
      build_large_cell(&cl_new,&cl_in,Tprel,is_noeq);
   } else  { /*  primitive lattice for output  */
      get_Telpr(lat_new,Telpr);
      mul_mm(Tp,Tprel,Telpr); /*  get new primitive basis from orientation of new elem.  */
      build_large_cell(&cl_new,&cl_in,Tp,1);
      mul_mm(Tin,Tin,Telpr);
   }
/*  shift origin of the cell  */
   shift_cell(Rsh[8],&cl_in);
   split_cell(is_prim,&cl_in,&cl_new,lat_new,nop,sym_op,name_srt,name_srt_new);

   transp(Tin,Tin);
   for(i=0; i<3; i++) a[0][i]=dvec(Tin[i]);
   a[1][0]=180./PI*acos(mul_vv(Tin[1],Tin[2])/(a[0][1]*a[0][2]));
   a[1][1]=180./PI*acos(mul_vv(Tin[0],Tin[2])/(a[0][0]*a[0][2]));
   a[1][2]=180./PI*acos(mul_vv(Tin[0],Tin[1])/(a[0][0]*a[0][1]));

   transp(Tin,Tin);
/*  decomposition of new basis vectors over oldest ones  */
   inv_matr(Tin_old,Tin_old);
   mul_mm(Tin_old,Tin_old,Tin);
   
   round_to_zero(sym_op[0][0], nop*3*4);
   round_to_zero(Rsh[0], nsh*3);
   round_to_zero(cl_in.r[0], cl_in.nat*3);
   round_to_zero(cl_new.r[0], cl_new.nat*3);
   round_to_zero(Tin_old[0],9);
     
   write_res(file_out,&cl_in,name_srt,&cl_new,name_srt_new,
	     is_noeq,is_prim,fix_org,is_wien_rhomb,
	     nsh,Rsh,a,lat_new,lat_name,&sgrp_name[0],
             npgrp,nop,sym_op,name_pgrp,Tin_old);
   
/*  find point group for nonequivalent atoms (all sorts)  */
/*   for(i=j=0; i<cl_new.nat; i++)
     if( cl_new.srt[i] == j ) {
      j++;
      find_atom_pgroup(lat_new, nop, npgrp, sym_op, cl_new.r[i], &npgrp_at, pgrp_at, Tel);
     }*/

   return(0);
}
   

