#ifndef IO_H
#define IO_H
/* read data from WIEN input file */
extern int read_wien_data(char *fname, t_cell *cl_in, char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
	      int *lat, int *is_wien_rhomb);

extern int read_data(char *fname, t_cell *cl_in, char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
	      int *lat) ;

extern int write_res(char *fname, t_cell *cl, char name_srt[MAX_ATOMS][MAX_CHARS],
	      t_cell *cl_new, char name_srt_new[MAX_ATOMS][MAX_CHARS],
	      int is_noeq, int is_prim, int fix_org, int is_wien_rhomb,
	      int nsh, double Rsh[][3],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int lat, char lat_name[14][32], char * sgrp_name,
              int npgrp, int nop, double sym_op[][4][3],
              char name_pgrp[][32],
	      double Tin[3][3]);

#endif
