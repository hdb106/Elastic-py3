#ifndef    TYPE_SG_H
#define    TYPE_SG_H

#ifndef PI
#define PI 3.14159265358979323844
#endif

#define  MAX_ATOMS    1024
#define  MAX_CHARS    32
#define  TOL          1e-4

/* type definition for cell */
typedef struct {
   int     nat;                    /*  number of atoms  */
   int     nsrt;                   /*  number of species  */
   double  r[MAX_ATOMS][3];        /*  coordinates  */
   int     srt[MAX_ATOMS];         /*  sort for each atom  */
}  t_cell;

#define DEBUG 0

#endif
