#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "type_sg.h"
#include "math_sg.h"
#include "pgrp_dat.h"
#include "sto.h"

int get_Telpr(int lat, double Telpr[3][3]);
void det_lat_NSM(int *lat, double T1[3], double T2[3], double T3[3],
                 double T1p[3], double T2p[3], double T3p[3],
        	 double T1t[3], double T2t[3], double T3t[3] );

/***************************************************
 * check T to be some translation vector
 * *************************************************/
int is_transl( t_cell *cl, double tr[3],
              double T[3][3], double T_1[3][3])
{
   int i,at;
   double r[3];

   for(at=0; at < cl->nat; at++)   {
      add_vv(r,cl->r[at],tr);
      put_in_cell(r,T,T_1);
      /*  find this atom  */
      for(i=0; i < cl->nat; i++)
        if( cl->srt[i] == cl->srt[at] &&
            ddvec(r,cl->r[i]) < TOL
          ) goto FOUND;
      return(0); /*  not found  */
      FOUND:;
   }
   return(1);
}

/**************************************************************
* determine Bravais lattice and translation vectors T1,T2,T3
* T1p,T2p,T3p - primitive lattice vectors
* on input T1,T2,T3 have minimal length
***************************************************************/
void det_lat(int *lat, double T1[3], double T2[3], double T3[3],
             double T1p[3], double T2p[3], double T3p[3],
	     double T1t[3], double T2t[3], double T3t[3] )
{
   int nat=0,i1,i2,i3,j1,j2,j3;
   double d1,d2,d3,t12,t13,t23,dmin,v;
   double r[37][3],p1[3],p2[3],p3[3];

   for(i1= 0; i1<=2; i1++)
   for(i2=-2; i2<=2; i2++)
   for(i3=-2; i3<=2; i3++) {
       if( (i1==0 && i2==0 && i3==0) ||
	   (abs(i1)==2 && i2==0 && i3==0) ||
	   (i1==0 && abs(i2)==2 && i3==0) ||
           (i1==0 && i2==0 && abs(i3)==2) ||
	   abs(i1)+abs(i2)==4 || abs(i1)+abs(i3)==4 || abs(i2)+abs(i3)==4
	 ) continue;
       r[nat][0]=i1*T1[0]+i2*T2[0]+i3*T3[0];
       r[nat][1]=i1*T1[1]+i2*T2[1]+i3*T3[1];
       r[nat][2]=i1*T1[2]+i2*T2[2]+i3*T3[2];
       if( dvec(r[nat]) < TOL) {
	  fprintf(stderr,"Error1 in det_lat().\n");
	  exit(1);
       }
       mul_cv(p1,-1.,r[nat]);
       for(j1=0; j1<nat; j1++) /*  find inverted vector  */
 	 if( ddvec(p1,r[j1]) < TOL ) goto NEXT_VEC;
       if( nat==37 ) {
	  fprintf(stderr,"Error2 in det_lat(): nat==37\n");
	  exit(1);
       }
       nat++;
      NEXT_VEC:;
   }
   /*** form TRICLINIC basis ***/
   asgn(p1,T1); asgn(p2,T2); asgn(p3,T3);
   d1=dvec(p1); d2=dvec(p2); if( d1>d2 ) swap(p1,p2);
   d1=dvec(p1); d3=dvec(p3); if( d1>d3 ) swap(p1,p3);
   d2=dvec(p2); d3=dvec(p3); if( d2>d3 ) swap(p2,p3);
   mulv_vv(T1t,p1,p2); /*   right basis  */
   if( mul_vv(T1t,p3) < 0 )  mul_cv(p1,-1.,p1);
   for(i1=-1; i1<=1; i1++)
   for(i2=-1; i2<=1; i2++)
   for(i3=-1; i3<=1; i3++) {
      if( i1*i2*i3<=0 ) continue;
      mul_cv(T1t,i1,asgn(T1t,p1));
      mul_cv(T2t,i2,asgn(T2t,p2));
      mul_cv(T3t,i3,asgn(T3t,p3));
      t12=mul_vv(T1t,T2t); t13=mul_vv(T1t,T3t); t23=mul_vv(T2t,T3t);
      if( ( t12 > TOL  && t13 > TOL && t23 > TOL ) ||
	  ( (t12 <= 0. || fabs(t12<TOL)) &&
	    (t13 <= 0. || fabs(t13<TOL)) &&
	    (t23 <= 0. || fabs(t23<TOL))) ) goto TRICL_FOUND;
   }
   fprintf(stderr,"Error in det_lat: TRICL_FOUND\n");
   exit(1);
   TRICL_FOUND:
   /**** check for CUBIC lattice ****/
   d1=dvec(T1); d2=dvec(T2); d3=dvec(T3);
   t12= fabs(mul_vv(T1,T2));
   t13= fabs(mul_vv(T1,T3));
   t23= fabs(mul_vv(T2,T3));
   v=vol(T1,T2,T3);
   
   if( fabs(d1-d2)<TOL && fabs(d1-d3)<TOL ) {
       if( t12<TOL && t13<TOL && t23<TOL) {
          *lat=CUBIC_P;
          right_basis(T1,T2,T3);
          asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
          return;
       }
      /*  find 3 orthogonal vectors of equal length  */
      for(j1=0; j1<nat; j1++) {
          d1=dvec(r[j1]);
          if( fabs(d1*d1*d1-2.*v) > TOL &&
              fabs(d1*d1*d1-4.*v) > TOL ) continue;
          for(j2=j1+1; j2<nat; j2++) {
              d2=dvec(r[j2]);
              if( fabs(d1-d2) >TOL ||
                  fabs(mul_vv(r[j1],r[j2])) > TOL ) continue;
              for(j3=j2+1; j3<nat; j3++) {
                  d3=dvec(r[j3]);
                  if( fabs(d1-d3) >TOL ||
                      fabs(mul_vv(r[j1],r[j3])) > TOL ||
                      fabs(mul_vv(r[j2],r[j3])) > TOL
                    ) continue;
                  if( fabs(vol(r[j1],r[j2],r[j3])-2.*v) < TOL ) {
                     *lat=CUBIC_I;
                     asgn(T1,r[j1]); asgn(T2,r[j2]); asgn(T3,r[j3]);
                     right_basis(T1,T2,T3);
                     for(i1=0;i1<3; i1++) {
                         T1p[i1]=0.5*( -T1[i1] + T2[i1] + T3[i1] );
                         T2p[i1]=0.5*(  T1[i1] - T2[i1] + T3[i1] );
                         T3p[i1]=0.5*(  T1[i1] + T2[i1] - T3[i1] );
                     }
                     return;
                  }
                  if( fabs(vol(r[j1],r[j2],r[j3])-4.*v) < TOL ) {
                     *lat=CUBIC_F;
                     asgn(T1,r[j1]); asgn(T2,r[j2]); asgn(T3,r[j3]);
                     right_basis(T1,T2,T3);
                     for(i1=0;i1<3; i1++) {
                         T1p[i1]=0.5*( T2[i1] + T3[i1] );
                         T2p[i1]=0.5*( T1[i1] + T3[i1] );
                         T3p[i1]=0.5*( T1[i1] + T2[i1] );
                     }
                     return;
                  }
              }
          }
      }
   }
   /***  check for HEXAGONAL lattice ****/
   /*  find 2 vectors with angle 2PI/3 and equal length  */
   dmin=0.;
   for(i1=0; i1<nat; i1++)
   for(i2=i1+1; i2<nat; i2++) {
      d1=dvec(r[i1]); d2=dvec(r[i2]);
      if( fabs(d1-d2) > TOL ||
          fabs(mul_vv(r[i1],r[i2])/(d1*d1)+0.5) > TOL ) continue;
      if( (d1 < dmin && fabs(d1-dmin) > TOL) || dmin < TOL ) {
         j1=i1; j2=i2; dmin=d1;
      }
   }
   if( dmin < TOL ) goto CHK_TETRAGONAL;
   /*  find c vector to be orthogonal to (a,a)  */
   for(j3=0; j3<nat; j3++) {
     if(fabs(mul_vv(r[j1],r[j3])) < TOL &&
        fabs(mul_vv(r[j2],r[j3])) < TOL &&
        fabs(vol(r[j1],r[j2],r[j3])-v) < TOL )
     {
        *lat=HEXAGONAL;
        asgn(T1,r[j1]); asgn(T2,r[j2]); asgn(T3,r[j3]);
        right_basis(T1,T2,T3);
        asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
        return;
     }
   }
   
 CHK_TETRAGONAL:
   /***  check for TETRAGONAL lattice ****/
   d1=dvec(T1); d2=dvec(T2); d3=dvec(T3);
   if( fabs(mul_vv(T1,T2)) <TOL &&
       fabs(mul_vv(T1,T3)) <TOL &&
       fabs(mul_vv(T2,T3)) <TOL &&
       (fabs(d1-d2)<TOL || fabs(d1-d3)<TOL || fabs(d2-d3)<TOL)
     ) { /*  primitive  */
         *lat=TETRAGONAL_P;
          if( fabs(d1-d3)<TOL ) swap(T2,T3);
          else
          if( fabs(d2-d3)<TOL ) swap(T1,T3);
          right_basis(T1,T2,T3);
          asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
          return;
       }
   /*  find 2 orthogonal vectors of equal and minimal length (a=b)  */
   dmin=0.;
   for(i1=0; i1<nat; i1++)
   for(i2=i1+1; i2<nat; i2++) {
      d1=dvec(r[i1]); d2=dvec(r[i2]);
      if( fabs(mul_vv(r[i1],r[i2])) > TOL ||
          fabs(d1-d2) > TOL ) continue;
      if( (d1 < dmin && fabs(d1-dmin) > TOL) || dmin < TOL ) {
         j1=i1; j2=i2; dmin=d1;
      }
   }
   if( dmin < TOL ) goto CHK_ORTHOROMBIC;
   /*  find c vector to be orthogonal to (a,a)  */
   for(j3=0; j3<nat; j3++)
   if(fabs(mul_vv(r[j1],r[j3])) < TOL &&
      fabs(mul_vv(r[j2],r[j3])) < TOL &&
      fabs(vol(r[j1],r[j2],r[j3])-2.*v) < TOL /*  vol_I=2*vol_P  */
     ) { /*  base centered ?  */
      add_vv(p1,r[j1],r[j2]); add_vv(p1,p1,r[j3]); mul_cv(p1,0.5,p1);
      if( find_vec_ivec(p1,r,nat) != -1 ) {
            *lat=TETRAGONAL_I;
             asgn(T1,r[j1]); asgn(T2,r[j2]); asgn(T3,r[j3]);
             right_basis(T1,T2,T3);
             for(i1=0;i1<3; i1++) {
                 T1p[i1]=0.5*( -T1[i1] + T2[i1] + T3[i1] );
                 T2p[i1]=0.5*(  T1[i1] - T2[i1] + T3[i1] );
                 T3p[i1]=0.5*(  T1[i1] + T2[i1] - T3[i1] );
             }
             return;
	 }
   }
   
 CHK_ORTHOROMBIC:
   /**** check for ORTHOROMBIC lattice ****/
   /*  find 3 orthogonal vectors  */
   dmin=0.;
   for(i1=0; i1<nat; i1++)
   for(i2=i1+1; i2<nat; i2++) {
      if( fabs(mul_vv(r[i1],r[i2])) > TOL ) continue;
      for(i3=i2+1; i3<nat; i3++) {
        if(fabs(mul_vv(r[i1],r[i3])) > TOL || fabs(mul_vv(r[i2],r[i3])) > TOL ) continue;
        dmin=vol(r[i1],r[i2],r[i3]);
        if( fabs(dmin-v)>TOL && fabs(dmin-2.*v)>TOL && fabs(dmin-4.*v)>TOL ) continue;
        asgn(T1,r[i1]); asgn(T2,r[i2]); asgn(T3,r[i3]);
        right_basis(T1,T2,T3);
        /*  check for primitive  */
        if( fabs(dmin-v) < TOL ) {
           *lat=ORTHOROMBIC_P;
           asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
           return;
        }
        /*  check for all faces centered  */
        if( fabs(dmin-4.*v) < TOL ) {
            for(j1=0;j1<3; j1++) {
               T1p[j1]=0.5*( T2[j1] + T3[j1] );
               T2p[j1]=0.5*( T1[j1] + T3[j1] );
               T3p[j1]=0.5*( T1[j1] + T2[j1] );
            }
	    if( find_vec_ivec(T1p,r,nat) != -1  &&
	        find_vec_ivec(T2p,r,nat) != -1  &&
	        find_vec_ivec(T3p,r,nat) != -1     )
	    {
	       *lat=ORTHOROMBIC_F;
	       return;
	    }
        }
        /*  check for body and base centred  */
        if( fabs(dmin-2.*v) < TOL ) {
	    add_vv(p1,T1,T2); add_vv(p1,p1,T3); mul_cv(p1,0.5,p1);
	    if( find_vec_ivec(p1,r,nat) != -1 ) {
	       *lat=ORTHOROMBIC_I;
	       for(j2=0;j2<3; j2++) {
		  T1p[j2]=0.5*( -T1[j2] + T2[j2] + T3[j2] );
		  T2p[j2]=0.5*(  T1[j2] - T2[j2] + T3[j2] );
		  T3p[j2]=0.5*(  T1[j2] + T2[j2] - T3[j2] );
	       }
              return;
            }
        /*  check for base centred  */
            for(j1=0; j1<3; j1++) {
                p1[j1]=0.5*(T2[j1]+T3[j1]); /*  A-centred  */
                p2[j1]=0.5*(T1[j1]+T3[j1]); /*  B-centred  */
                p3[j1]=0.5*(T1[j1]+T2[j1]); /*  C-centred  */
            }
 	    if( find_vec_ivec(p1,r,nat) != -1 ) { /*  A-centred  */
                swap(T1,T2); swap(T2,T3); goto ORTH_C;
            } else
            if( find_vec_ivec(p2,r,nat) != -1 ) { /*  B-centred  */
                swap(T1,T3); swap(T2,T3); goto ORTH_C;
            } else
            if( find_vec_ivec(p3,r,nat) != -1 ) goto ORTH_C;
	    continue;

           ORTH_C:
           *lat=ORTHOROMBIC_C;
	   sub_vv(T1p,T1,T2); mul_cv(T1p,0.5,T1p); /*  T1p=1/2*(T1-T2)  */
	   add_vv(T2p,T1,T2); mul_cv(T2p,0.5,T2p); /*  T2p=1/2*(T1+T2)  */
	   asgn(T3p,T3);                           /*  T3p=T3  */
	   return;
	}
      }
   }
/*     CHK_RHOMBOHEDRAL:  */
   /**** check for RHOMBOHEDRAL lattice ****/
   /*  find 3 vectors of equal length and equal angles  */
   for(j1=0; j1<nat; j1++) {
       d1=dvec(r[j1]);
       for(j2=j1+1; j2<nat; j2++) {
           d2=dvec(r[j2]);
           if( fabs(d1-d2) >TOL ) continue;
	   mul_cv(p2,-1.,r[j2]);
           for(j3=j2+1; j3<nat; j3++) {
               d3=dvec(r[j3]);
               if( fabs(d1-d3) > TOL || fabs(vol(r[j1],r[j2],r[j3])-v) > TOL ) continue;
	       mul_cv(p3,-1.,r[j3]);
	
	       t12=mul_vv(r[j1],r[j2]); t13=mul_vv(r[j1],r[j3]); t23=mul_vv(r[j2],r[j3]);
               if( fabs(t12-t13) < TOL  &&  fabs(t12-t23) < TOL ) {
		  asgn(T1p,r[j1]);
		  asgn(T2p,r[j2]);
		  asgn(T3p,r[j3]);
		  goto RHOMB;
	       }
	       t12=mul_vv(r[j1],p2); t13=mul_vv(r[j1],r[j3]); t23=mul_vv(p2,r[j3]);
               if( fabs(t12-t13) < TOL  &&  fabs(t12-t23) < TOL ) {
		  asgn(T1p,r[j1]);
		  asgn(T2p,p2);
		  asgn(T3p,r[j3]);
		  goto RHOMB;
	       }
	       t12=mul_vv(r[j1],r[j2]); t13=mul_vv(r[j1],p3); t23=mul_vv(r[j2],p3);
               if( fabs(t12-t13) < TOL  &&  fabs(t12-t23) < TOL ) {
		  asgn(T1p,r[j1]);
		  asgn(T2p,r[j2]);
		  asgn(T3p,p3);
		  goto RHOMB;		
	       }
	       t12=mul_vv(r[j1],p2); t13=mul_vv(r[j1],p3); t23=mul_vv(p2,p3);
               if( fabs(t12-t13) < TOL  &&  fabs(t12-t23) < TOL ) {
		  asgn(T1p,r[j1]);
		  asgn(T2p,p2);
		  asgn(T3p,p3);
		  goto RHOMB;
	       }
	       continue;
	      RHOMB:
                  *lat=RHOMBOHEDRAL;
                  right_basis(T1p,T2p,T3p);
                  for(i1=0;i1<3; i1++) { /*   hexagonal setting  */
                     T1[i1]= T1p[i1] - T2p[i1];
                     T2[i1]= T2p[i1] - T3p[i1];
                     T3[i1]= T1p[i1] + T2p[i1] + T3p[i1] ;
		  }
                  return;
	   }
       }
   }

/*     CHK_MONOCLINIC:  only c-type setting: (c,a)=(c,b)=0  */
   /**** check for MONOCLINIC lattice ****/
   for(j3=0; j3<nat; j3++)
   for(j2=j3+1; j2<nat; j2++) {
     if( fabs(mul_vv(r[j3],r[j2])) > TOL ) continue;
     for(j1=j3+1; j1<nat; j1++) {
       if( fabs(mul_vv(r[j3],r[j1])) < TOL && vol(r[j1],r[j2],r[j3]) > TOL ) {
           asgn(T1,r[j1]); asgn(T2,r[j2]); asgn(T3,r[j3]);
       } else
       if( fabs(mul_vv(r[j2],r[j1])) < TOL && vol(r[j1],r[j2],r[j3]) > TOL ) {
           asgn(T1,r[j1]); asgn(T2,r[j3]); asgn(T3,r[j2]);
       } else continue;
	
       right_basis(T1,T2,T3); /*  right_basis function swaps only 1 and 2 vector  */
       if( fabs(vol(T1,T2,T3)-v) < TOL ) { /*  monoclinic primitive  */
          *lat=MONOCLINIC_P;
           asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
           return;
       }
       if( fabs(vol(T1,T2,T3)-2.*v) < TOL ) { /*  monoclinic base centred  */
          *lat=MONOCLINIC_A;
           for(i1=0; i1<3; i1++) {
             p1[i1]=0.5*(T2[i1]+T3[i1]); /*  A-centred  */
             p2[i1]=0.5*(T1[i1]+T3[i1]); /*  B-centred  */
           }
	   if( find_vec_ivec(p1,r,nat) != -1 ) { /*  A-centred  */
      	      asgn(T1p,T1);                           /*  T1p=T1  */
	      add_vv(T2p,T2,T3); mul_cv(T2p,0.5,T2p); /*  T2p=1/2*( T2+T3)  */
	      sub_vv(T3p,T3,T2); mul_cv(T3p,0.5,T3p); /*  T3p=1/2*(-T2+T3)  */
	      return;
	   } else
	   if( find_vec_ivec(p2,r,nat) != -1 ) { /*  B-centred  */
              mul_cv(T2,-1.,T2); /*  transform to A-centred  */
              swap(T1,T2);
      	      asgn(T1p,T1);                           /*  T1p=T1  */
	      add_vv(T2p,T2,T3); mul_cv(T2p,0.5,T2p); /*  T2p=1/2*( T2+T3)  */
	      sub_vv(T3p,T3,T2); mul_cv(T3p,0.5,T3p); /*  T3p=1/2*(-T2+T3)  */
              return;
	   }
       }/*  end of monoclinic base centred  */
     }
   } /*   END OF CHK_MONOCLINIC section  */
      
   /**** lattice is TRICLINIC ****/
  *lat=TRICLINIC;
   asgn(T1,T1t); asgn(T2,T2t); asgn(T3,T3t);
   asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
   return;

}

/***************************************************
 * find basis translations
 * Ttrcl -  triclinic basis
 ***************************************************/
void basis_transl(t_cell *cl_in, t_cell *cl, int *lat, double Tin[3][3],
                  double Tp[3][3], double T[3][3], double Ttrcl[3][3] )
{
   int at,fl_ch,i1,i2,i3;
   double d,dmin,d1,d2,d3,v,p[3],r[3];
   double T1[3],T2[3],T3[3],
          T1p[3],T2p[3],T3p[3],Tin_1[3][3];
   
   inv_matr(Tin,Tin_1);
   /*  find minimal transl. vector to be parallel cl->T1  */
   colmn(T1,Tin,0);
   p[0]=-T1[1];   p[1]= T1[0];   p[2]= T1[2];
   norm_vec(p);
   for(at=0; at < cl->nat; at++) {
     d=dvec(cl->r[at]);
     if( d < TOL                          || /*  zero length  */
         fabs(cl->r[at][2])  > TOL        || /*   z != 0  */
         fabs(mul_vv(p,cl->r[at])) > TOL  || /*   r is not parallel T1  */
         dvec(T1) <= d
       ) continue;
     if( is_transl(cl_in,cl->r[at],Tin,Tin_1) ) asgn(T1,cl->r[at]);
   }
   /*   find minimal transl. vector T2  */
   /*    PI/2 rotation of T1  */
   colmn(T2,Tin,1);
   dmin=dvec(T2);
   for(at=0; at < cl->nat; at++) {
     d=fabs(mul_vv(p,cl->r[at]));
     if( d < TOL || dmin <= d ||
         fabs(cl->r[at][2]) > TOL
        ) continue;
     if( is_transl(cl_in,cl->r[at],Tin,Tin_1) ) {
         asgn(T2,cl->r[at]);
         dmin=d;
     }
   }
   /*   find minimal transl. vector T3  */
   p[0]=0.;
   p[1]=0.;
   p[2]=1.;
   colmn(T3,Tin,2);
   dmin=dvec(T3);
   for(at=0; at < cl->nat; at++) {
     d=fabs(mul_vv(p,cl->r[at]));
     if( d < TOL || dmin <= d ) continue;
     if( is_transl(cl_in,cl->r[at],Tin,Tin_1) ) {
         asgn(T3,cl->r[at]);
         dmin=d;
     }
   }
#if DEBUG==1
   mul_mv(T1p,Tin_1,T1);
   mul_mv(T2p,Tin_1,T2);
   mul_mv(T3p,Tin_1,T3);
#endif
   /* reduce T1,T2,T3 to minimal length */
   d1=dvec(T1); d2=dvec(T2); d3=dvec(T3);
   v=vol(T1,T2,T3);
   do {
      fl_ch=0;
      for(i1=-1; i1<=1; i1++)
      for(i2=-1; i2<=1; i2++)
      for(i3=-1; i3<=1; i3++) {
         r[0]=i1*T1[0]+i2*T2[0]+i3*T3[0];
         r[1]=i1*T1[1]+i2*T2[1]+i3*T3[1];
         r[2]=i1*T1[2]+i2*T2[2]+i3*T3[2];
         d=dvec(r);
         if( d<d1 && fabs(d-d1) > TOL &&
             fabs(vol(r,T2,T3)-v) < TOL ) {
             asgn(T1,r);
             d1=d;
             fl_ch=1;
         } else
         if( d<d2 && fabs(d-d2) > TOL &&
             fabs(vol(T1,r,T3)-v) < TOL ) {
             asgn(T2,r);
             d2=d;
             fl_ch=1;
         } else
         if( d<d3 && fabs(d-d3) > TOL &&
             fabs(vol(T1,T2,r)-v) < TOL ) {
             asgn(T3,r);
             d3=d;
             fl_ch=1;
         }
      }
   } while(fl_ch);
#if DEBUG==1
/*     printf("LAT= %s\n",lat_name[lat]);  */
{
   int i;
   print_dbl("Tin=\n",3,Tin[0]);
   print_dbl(NULL,3,Tin[1]);
   print_dbl(NULL,3,Tin[2]);
   
   print_dbl("T1=",3,T1);
   print_dbl("T2=",3,T2);
   print_dbl("T3=",3,T3);
   
   mul_mv(T1p,Tin_1,T1);
   mul_mv(T2p,Tin_1,T2);
   mul_mv(T3p,Tin_1,T3);
}
#endif
   /* determine Brave's lattice */
/*     det_lat(lat,T1,T2,T3,T1p,T2p,T3p,Ttrcl[0],Ttrcl[1],Ttrcl[2]);  */
   det_lat_NSM(lat,T1,T2,T3,T1p,T2p,T3p,Ttrcl[0],Ttrcl[1],Ttrcl[2]);
   asgn(T[0],T1); asgn(Tp[0],T1p);
   asgn(T[1],T2); asgn(Tp[1],T2p);
   asgn(T[2],T3); asgn(Tp[2],T3p);
   transp(T,T);   transp(Tp,Tp);  transp(Ttrcl,Ttrcl);

#if DEBUG==1
/*     printf("LAT= %s\n",lat_name[lat]);  */
{
	int i;
   printf("Tin=\n");
   for(i=0; i<3; i++)
     printf(" % .8f, % .8f, % .8f\n",Tin[i][0],Tin[i][1],Tin[i][2]);
   
   printf("\nT1= % .8f, % .8f, % .8f\n",T1[0],T1[1],T1[2]);
   printf("T2= % .8f, % .8f, % .8f\n",T2[0],T2[1],T2[2]);
   printf("T3= % .8f, % .8f, % .8f\n",T3[0],T3[1],T3[2]);
   printf("\nT1p= % .8f, % .8f, % .8f\n",T1p[0],T1p[1],T1p[2]);
   printf("T2p= % .8f, % .8f, % .8f\n",T2p[0],T2p[1],T2p[2]);
   printf("T3p= % .8f, % .8f, % .8f\n",T3p[0],T3p[1],T3p[2]);
}
#endif
/*     exit(0);  */
} /*    reduce_trans()  */

void new_basis( double a11, double a12, double a13,
	        double a21, double a22, double a23,
 	        double a31, double a32, double a33,
	        double T1[3], double T2[3], double T3[3] )
{
   double T_on[3][3],T_xo[3][3];

   T_on[0][0]=a11;  T_on[0][1]=a12;  T_on[0][2]=a13;
   T_on[1][0]=a21;  T_on[1][1]=a22;  T_on[1][2]=a23;
   T_on[2][0]=a31;  T_on[2][1]=a32;  T_on[2][2]=a33;
   transp(T_on,T_on);
   
   asgn(T_xo[0],T1); asgn(T_xo[1],T2); asgn(T_xo[2],T3);
   transp(T_xo,T_xo);
   
   mul_mm(T_xo, T_xo, T_on);
   transp(T_xo,T_xo);
   asgn(T1,T_xo[0]); asgn(T2,T_xo[1]); asgn(T3,T_xo[2]);
   return;
}
   
#define ERR_main { fprintf(stderr,"Error: the main conditions not satisfied.\n"); exit(1); }
#define ERR_spec { fprintf(stderr,"Error: the special conditions not satisfied.\n"); exit(1); }

/**************************************************************
* determine Bravais lattice and translation vectors T1,T2,T3
* T1p,T2p,T3p - primitive lattice vectors
* on input T1,T2,T3 have minimal lengths,
* using table given by A.Santoro, A.D. Mighel, Determination of reduced cell,
* Acta Cryst. A26, p.124., to satisfy the special conditions
*
* and table for 44 Niggli's transformations of reduced basis vectors
* Int. Tables (1992), table 9.3.1
* to conventional basis
***************************************************************/
void det_lat_NSM(int *lat, double T1[3], double T2[3], double T3[3],
/*  	     double t1[3], double t2[3], double t3[3],  */
             double T1p[3], double T2p[3], double T3p[3],
	     double T1t[3], double T2t[3], double T3t[3] )
{
   int i,i1,i2,i3,itype,updated;
   double S11,S22,S33,S12,S13,S23, A,B,C,D,E,F;
   double T_on[3][3];
/*     double T1[3]; double T2[3]; double T3[3];  */

/*     asgn(T1,t1);    asgn(T2,t2);     asgn(T3,t3);  */
   if( gt(dvec(T1),dvec(T2)) ) swap(T1,T2);
   if( gt(dvec(T2),dvec(T3)) ) swap(T2,T3);
   if( gt(dvec(T1),dvec(T2)) ) swap(T1,T2);
   asgn(T_on[0],T1); asgn(T_on[1],T2); asgn(T_on[2],T3);
   if( lt(detr(T_on),0) ) mul_cv(T1,-1,T1); /*  right handed basis  */
   
/*** See, A.Santoro, A.D. Mighel, Determination of reduced cell, ***/
/*** Acta Cryst. A26, p.124.                                     ***/
   for(i=0; i<100; i++) {
      updated=itype=0;
      for(i1=1; i1 >= -1; i1-=2)
      for(i2=1; i2 >= -1; i2-=2)
      for(i3=1; i3 >= -1; i3-=2) {
	 if( i1*i2*i3 < 0 ) continue;
	 mul_cv(T_on[0],i1,T1);
	 mul_cv(T_on[1],i2,T2);
	 mul_cv(T_on[2],i3,T3);
	 S12=mul_vv(T_on[0],T_on[1]);
	 S13=mul_vv(T_on[0],T_on[2]);
	 S23=mul_vv(T_on[1],T_on[2]);
	 if( gt(S12,0) && gt(S13,0) && gt(S23,0) ) {
	    itype=1;
	    asgn(T1,T_on[0]); asgn(T2,T_on[1]); asgn(T3,T_on[2]);
	    goto FOUND;
	 }
	 if( le(S12,0) && le(S13,0) && le(S23,0) ) {
	    itype=2;
	    asgn(T1,T_on[0]); asgn(T2,T_on[1]); asgn(T3,T_on[2]);
	    goto FOUND;
	 }
      }
      FOUND:;
      if( itype==0 ) {
	 fprintf(stderr,"Error: basis can not satisfy the main conditions.\n");
	 exit(1);
      }
      
      S11=mul_vv(T1,T1); S22=mul_vv(T2,T2); S33=mul_vv(T3,T3);
      if( eq(S11,S22) && lt(fabs(S13),fabs(S23)) ) {
	  new_basis( 0, -1, 0,  -1, 0, 0,   0, 0, -1, T1,T2,T3);
	  updated=1;
      } else
      if( eq(S22,S33) && lt(fabs(S12),fabs(S13)) ) {
	  new_basis( -1, 0, 0,   0, 0,-1,   0,-1, 0, T1,T2,T3);
	  updated=1;
      } else /*  type I  */
      if( itype==1 && eq(S23,S22/2.) &&  lt(2.*S13, S12)) {
	  new_basis( -1, 0, 0,   0,-1, 0,   0,-1, 1, T1,T2,T3);
	  updated=1;
      } else
      if( itype==1 && eq(S13,S11/2.) && lt(2.*S23, S12)) {
 	  new_basis( -1, 0, 0,   0,-1, 0,  -1, 0, 1, T1,T2,T3);
	  updated=1;
      } else
      if( itype==1 && eq(S12,S11/2.) && lt(2.*S23, S13)) {
	  new_basis( -1, 0, 0,  -1, 1, 0,   0, 0,-1, T1,T2,T3);
	  updated=1;
      } else /*  type II  */
      if( itype==2 && eq(S23, -S22/2.) && ne(S12,0) ) {
	  new_basis(  1, 0, 0,   0,-1, 0,   0,-1,-1, T1,T2,T3);
	  updated=1;
      } else
      if( itype==2 && eq(S13, -S11/2.) && ne(S12,0) ) {
	  new_basis( -1, 0, 0,   0, 1, 0,  -1, 0,-1, T1,T2,T3);
	  updated=1;
      } else
      if( itype==2 && eq(S12, -S11/2.) && ne(S13,0) ) {
	  new_basis( -1, 0, 0,  -1,-1, 0,   0, 0, 1, T1,T2,T3);
	  updated=1;
      } else
      if( itype==2 && eq(S23,(S11+S22-2*fabs(S13)-2*fabs(S12))/2.) &&
	  lt(2*fabs(S13)+fabs(S12),S11) ) {
	  new_basis( -1, 0, 0,   0,-1, 0,   1, 1, 1, T1,T2,T3);
	  updated=1;
     }
     if( updated==0 ) break;
   }
   if( i==100 ) {
      fprintf(stderr,"Error: basis can not satisfy the special"
	              "conditions within 100 iterations.\n");
      exit(1);
   }
   
   /*  check again the main and special conditions  */
   S11=mul_vv(T1,T1); S22=mul_vv(T2,T2); S33=mul_vv(T3,T3);
   S12=mul_vv(T1,T2); S13=mul_vv(T1,T3); S23=mul_vv(T2,T3);
   if( gt(S12,0) && gt(S13,0) && gt(S23,0) ) itype=1;
   else itype=2;
   if( itype==1 ) {
     /*  main conditions  */
     if(!( le(S11,S22) && le(S22,S33) &&
	   le(S23,S22/2) && le(S13,S11/2) && le(S12,S11/2) )) ERR_main;
     /*  special conditions  */
      if( eq(S11,S22)   ) if( ! le(S23, S13) ) ERR_spec;
      if( eq(S22,S33)   ) if( ! le(S13, S12) ) ERR_spec;
      if( eq(S23,S22/2) ) if( ! le(S12,2*S13)) ERR_spec;
      if( eq(S13,S11/2) ) if( ! le(S12,2*S23)) ERR_spec;
      if( eq(S12,S11/2) ) if( ! le(S13,2*S23)) ERR_spec;
   } else { /*  type==2  */
     /*  main conditions  */
     if(!( le(S11,S22) && le(S22,S33) &&
	   le(fabs(S23),S22/2) && le(fabs(S13),S11/2) && le(fabs(S12),S11/2) &&
	   le(fabs(S23)+fabs(S13)+fabs(S12),(S11+S22)/2) )) ERR_main;
     /*  special conditions  */
      if( eq(S11,S22)   )       if( ! le(fabs(S23), fabs(S13)) ) ERR_spec;
      if( eq(S22,S33)   )       if( ! le(fabs(S13), fabs(S12)) ) ERR_spec;
      if( eq(fabs(S23),S22/2) ) if( ! le(S12,0) )  ERR_spec;
      if( eq(fabs(S13),S11/2) ) if( ! le(S12,0) )  ERR_spec;
      if( eq(fabs(S12),S11/2) ) if( ! le(S13,0) )  ERR_spec;
      if( eq(fabs(S23)+fabs(S13)+fabs(S12),(S11+S22)/2) )
 	  if( ! le(S11,2*fabs(S13)+fabs(S12)) ) ERR_spec;
   }
   
   /* table 9.3.1 from IT(1992) - 44 cases of Niggli's reduction */
   A=mul_vv(T1,T1); B=mul_vv(T2,T2); C=mul_vv(T3,T3);
   D=mul_vv(T2,T3); E=mul_vv(T1,T3); F=mul_vv(T1,T2);
 /*  A=B=C  */
   if( eq(A,B) && eq(A,C) ) {
      if( itype==1 && eq(D,A/2) && eq(E,A/2) && eq(F,A/2) ) {   /*  1 CUBIC_F  */
	 new_basis( 1,-1, 1,  1, 1,-1, -1, 1, 1, T1,T2,T3);
	 *lat=CUBIC_F; goto F_44;
      } else
      if( itype==1 && eq(E,D) && eq(F,D) ) {   /*  2 RHOMBOHEDRAL  */
	 new_basis( 1,-1, 0, -1, 0, 1, -1,-1,-1, T1,T2,T3);
	 *lat=RHOMBOHEDRAL; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) && eq(F,0) ) {   /*  3 CUBIC_P  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=CUBIC_P; goto F_44;
      } else
      if( itype==2 && eq(D,-A/3.) && eq(E,-A/3.) && eq(F,-A/3.) ) { /*  5 CUBIC_I  */
	 new_basis( 1, 0, 1,  1, 1, 0,  0, 1, 1, T1,T2,T3);
	 *lat=CUBIC_I; goto F_44;
      } else
      if( itype==2 && eq(E,D) && eq(F,D) ) {   /*  4 RHOMBOHEDRAL  */
	 new_basis( 1,-1, 0, -1, 0, 1, -1,-1,-1, T1,T2,T3);
	 *lat=RHOMBOHEDRAL; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) && eq(E,D) ) { /*  6 TETRAGONAL_I  */
	 new_basis( 0, 1, 1,  1, 0, 1,  1, 1, 0, T1,T2,T3);
	 *lat=TETRAGONAL_I; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) && eq(F,E) ) { /*  7 TETRAGONAL_I  */
	 new_basis( 1, 0, 1,  1, 1, 0,  0, 1, 1, T1,T2,T3);
	 *lat=TETRAGONAL_I; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) ) { /*  8 ORTHOROMBIC_I  */
	  new_basis(-1,-1, 0, -1, 0,-1,  0,-1,-1, T1,T2,T3);
	 *lat=ORTHOROMBIC_I; goto F_44;
      }
   }
 /*  A=B  */
   if( eq(A,B) ) {
      if( itype==1 && eq(D,A/2) && eq(E,A/2) && eq(F,A/2) ) { /*  9 RHOMBOHEDRAL  */
	 new_basis( 1, 0, 0, -1, 1, 0, -1,-1, 3, T1,T2,T3);
	 *lat=RHOMBOHEDRAL; goto F_44;
      } else
      if( itype==1 && eq(E,D) ) { /*  10 MONOCLINIC_C  */
	 new_basis( 1, 1, 0,  1,-1, 0,  0, 0,-1, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) && eq(F,0) ) {   /*  11 TETRAGONAL_P  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=TETRAGONAL_P; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) && eq(F,-A/2.) ) {   /*  12 HEXAGONAL  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=HEXAGONAL; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) ) {   /*  13 ORTHOROMBIC_C  */
	 new_basis( 1, 1, 0, -1, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=ORTHOROMBIC_C; goto F_44;
      } else
      if( itype==2 && eq(D,-A/2) && eq(E,-A/2) && eq(F,0) ) { /*  15 TETRAGONAL_I  */
	 new_basis( 1, 0, 0,  0, 1, 0,  1, 1, 2, T1,T2,T3);
	 *lat=TETRAGONAL_I; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) && eq(E,D) ) { /*  16 ORTHOROMBIC_F  */
	  new_basis(-1,-1, 0,  1,-1, 0,  1, 1, 2, T1,T2,T3);
	 *lat=ORTHOROMBIC_F; goto F_44;
      } else
      if( itype==2 && eq(E,D) ) { /*  14 MONOCLINIC_C  */
	 new_basis( 1, 1, 0, -1, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) ) { /*  17 MONOCLINIC_C  */
	  new_basis( 1,-1, 0,  1, 1, 0, -1, 0,-1, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      }
   }
 /*  B=C  */
   if( eq(B,C) ) {
      if( itype==1 && eq(D,A/4) && eq(E,A/2) && eq(F,A/2) ) { /*  18 TETRAGONAL_I  */
	 new_basis( 0,-1, 1,  1,-1,-1,  1, 0, 0, T1,T2,T3);
	 *lat=TETRAGONAL_I; goto F_44;
      } else
      if( itype==1 && eq(E,A/2) && eq(F,A/2) ) { /*  19 ORTHOROMBIC_I  */
	 new_basis(-1, 0, 0,  0,-1, 1, -1, 1, 1, T1,T2,T3);
	 *lat=ORTHOROMBIC_I; goto F_44;
      } else
      if( itype==1 && eq(F,E) ) { /*  20 MONOCLINIC_C  */
	 new_basis( 0, 1, 1,  0, 1,-1, -1, 0, 0, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) && eq(F,0) ) {   /*  21 TETRAGONAL_P  */
	 new_basis( 0, 1, 0,  0, 0, 1,  1, 0, 0, T1,T2,T3);
	 *lat=TETRAGONAL_P; goto F_44;
      } else
      if( itype==2 && eq(D,-B/2) && eq(E,0) && eq(F,0) ) {   /*  22 HEXAGONAL  */
	 new_basis( 0, 1, 0,  0, 0, 1,  1, 0, 0, T1,T2,T3);
	 *lat=HEXAGONAL; goto F_44;
      } else
      if( itype==2 && eq(E,0) && eq(F,0) ) {   /*  23 ORTHOROMBIC_C  */
	 new_basis( 0, 1, 1,  0,-1, 1,  1, 0, 0, T1,T2,T3);
	 *lat=ORTHOROMBIC_C; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) && eq(E,-A/3) && eq(F,-A/3)  ) { /*  24 RHOMBOHEDRAL  */
	  new_basis( 1, 2, 1,  0,-1, 1,  1, 0, 0, T1,T2,T3);
	 *lat=RHOMBOHEDRAL; goto F_44;
      } else
      if( itype==2 && eq(F,E) ) { /*  25 MONOCLINIC_C  */
	 new_basis( 0, 1, 1,  0,-1, 1,  1, 0, 0, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      }
   }
/*   no conditions on A,B,C  */
      if( itype==1 && eq(D,A/4) && eq(E,A/2) && eq(F,A/2) ) { /*  26 ORTHOROMBIC_F  */
	 new_basis( 1, 0, 0, -1, 2, 0, -1, 0, 2, T1,T2,T3);
	 *lat=ORTHOROMBIC_F; goto F_44;
      } else
      if( itype==1 && eq(E,A/2) && eq(F,A/2) ) { /*  27 MONOCLINIC_C  */
	 new_basis(-1, 2, 0, -1, 0, 0,  0,-1, 1, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==1 && eq(E,A/2) && eq(F,2*D) ) { /*  28 MONOCLINIC_C  */
	 new_basis(-1, 0, 0, -1, 0, 2,  0, 1, 0, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==1 && eq(E,2*D) && eq(F,A/2) ) { /*  29 MONOCLINIC_C  */
	 new_basis( 1, 0, 0,  1,-2, 0,  0, 0,-1, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==1 && eq(D,B/2) && eq(F,2*E) ) { /*  30 MONOCLINIC_C  */
	 new_basis( 0, 1, 0,  0, 1,-2, -1, 0, 0, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==1 ) { /*  31 TRICLINIC  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=TRICLINIC; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) && eq(F,0) ) {   /*  32 ORTHOROMBIC_P  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);	
	 *lat=ORTHOROMBIC_P; goto F_44;
      } else
      if( itype==2 && eq(D,-B/2) && eq(E,0) && eq(F,0) ) {   /*  40 ORTHOROMBIC_C  */
	 new_basis( 0,-1, 0,  0, 1, 2, -1, 0, 0, T1,T2,T3);	
	 *lat=ORTHOROMBIC_C; goto F_44;
      } else
      if( itype==2 && eq(E,0) && eq(F,0) ) {   /*  35 MONOCLINIC_P  */
	 new_basis( 0,-1, 0, -1, 0, 0,  0, 0,-1, T1,T2,T3);	
	 *lat=MONOCLINIC_P; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,-A/2) && eq(F,0) ) {   /*  36 ORTHOROMBIC_C  */
	 new_basis( 1, 0, 0, -1, 0,-2,  0, 1, 0, T1,T2,T3);	
	 *lat=ORTHOROMBIC_C; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(F,0) ) {   /*  33 MONOCLINIC_P  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=MONOCLINIC_P; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) && eq(F,-A/2) ) {   /*  38 ORTHOROMBIC_C  */
	 new_basis(-1, 0, 0,  1, 2, 0,  0, 0,-1, T1,T2,T3);	
	 *lat=ORTHOROMBIC_C; goto F_44;
      } else
      if( itype==2 && eq(D,0) && eq(E,0) ) {   /*  34 MONOCLINIC_P  */
	 new_basis(-1, 0, 0,  0, 0,-1,  0,-1, 0, T1,T2,T3);
	 *lat=MONOCLINIC_P; goto F_44;
      } else
      if( itype==2 && eq(D,-B/2) && eq(E,-A/2) && eq(F,0) ) {   /*  42 ORTHOROMBIC_I  */
	 new_basis(-1, 0, 0,  0,-1, 0,  1, 1, 2, T1,T2,T3);
	 *lat=ORTHOROMBIC_I; goto F_44;
      } else
      if( itype==2 && eq(D,-B/2) && eq(F,0) ) {   /*  41 MONOCLINIC_C  */
	 new_basis( 0,-1,-2,  0,-1, 0, -1, 0, 0, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==2 && eq(E,-A/2) && eq(F,0) ) {   /*  37 MONOCLINIC_C  */
	 new_basis( 1, 0, 2,  1, 0, 0,  0, 1, 0, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==2 && eq(E,0) && eq(F,-A/2) ) {   /*  39 MONOCLINIC_C  */
	 new_basis(-1,-2, 0, -1, 0, 0,  0, 0,-1, T1,T2,T3);
	 *lat=MONOCLINIC_C; goto F_44;
      } else
      if( itype==2 && eq(2*fabs(D+E+F),A+B) && eq(fabs(2*D+F),B) ) { /*  43 MONOCLINIC_I  */
	  new_basis(-1, 0, 0, -1,-1,-2,  0,-1, 0, T1,T2,T3);
	 *lat=MONOCLINIC_I; goto F_44;
      } else
      if( itype==2 ) { /*  44 TRICLINIC  */
	 new_basis( 1, 0, 0,  0, 1, 0,  0, 0, 1, T1,T2,T3);
	 *lat=TRICLINIC; goto F_44;
      } else {
	 fprintf(stderr,"Error: no one choise found for Niggli's transformation.\n");
	 exit(1);
      }

   F_44:;
   /*  transform MONOCLIC b-setting to c-setting  */
   if( *lat==MONOCLINIC_P || *lat==MONOCLINIC_C || *lat==MONOCLINIC_I ) {
      new_basis(0,0,1,  1,0,0,  0,1,0, T1,T2,T3);
      if( *lat==MONOCLINIC_C ) *lat=MONOCLINIC_A;
      if( *lat==MONOCLINIC_I ) { /*  transform MONOCLINIC_I to MONOCLINIC_A  */
	 new_basis(0,-1,0, 1,1,0, 0,0,1, T1,T2,T3);	
	 *lat=MONOCLINIC_A;
      }
   }
   /*  obtain primitive vectors  */
   asgn(T1p,T1); asgn(T2p,T2); asgn(T3p,T3);
   get_Telpr(*lat,T_on);
   new_basis( T_on[0][0],T_on[1][0],T_on[2][0],
	      T_on[0][1],T_on[1][1],T_on[2][1],
	      T_on[0][2],T_on[1][2],T_on[2][2],
	      T1p,T2p,T3p);
   asgn(T1t,T1); asgn(T2t,T2); asgn(T3t,T3);
/*     asgn(t1,T1); asgn(t2,T2); asgn(t3,T3);  */
   return;
}
