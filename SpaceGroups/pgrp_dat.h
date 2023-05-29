#ifndef   PGRP_DAT_H
#define   PGRP_DAT_H

#define  TRICLINIC      0
#define  MONOCLINIC_P   1
#define  MONOCLINIC_A   2
#define  ORTHOROMBIC_P  3
#define  ORTHOROMBIC_I  4
#define  ORTHOROMBIC_C  5
#define  ORTHOROMBIC_A  6
#define  ORTHOROMBIC_F  7
#define  TETRAGONAL_P   8
#define  TETRAGONAL_I   9
#define  RHOMBOHEDRAL   10
#define  HEXAGONAL      11
#define  CUBIC_P        12
#define  CUBIC_I        13
#define  CUBIC_F        14
/*  lattice for input only  */
#define  ORTHOROMBIC_B  100
/*  lattice for reduction only  */
#define  MONOCLINIC_C  101
#define  MONOCLINIC_I  102

/*   numbers of point groups  */
#define   NC1        0
#define   NCi        1

#define   NC2        2
#define   NCs        3
#define   NC2h       4

#define   ND2        5
#define   NC2v       6
#define   ND2h       7

#define   NC4        8
#define   NS4        9
#define   NC4h       10
#define   ND4        11
#define   NC4v       12
#define   ND2d       13
#define   ND2d2      14
#define   ND4h       15


#define   NC3        16
#define   NC3i	     17
#define   ND3	     18
#define   ND32	     19
#define   NC3v	     20
#define   NC3v2	     21
#define   ND3d	     22
#define   ND3d2	     23

#define   NC6        24
#define   NC3h       25
#define   NC6h       26
#define   ND6        27
#define   NC6v       28
#define   ND3h       29
#define   ND3h2      30
#define   ND6h       31

#define   NT         32
#define   NTh        33
#define   NO         34
#define   NTd        35
#define   NOh        36

/*  symmetry operations  */
#define   ONE        0
 
#define   C2001      1
#define   C2010      2
#define   C2100      3

#define   C3p111     4
#define   C3p1_1_1   5
#define   C3p_11_1   6
#define   C3p_1_11   7

#define   C3m111     8
#define   C3m1_1_1   9
#define   C3m_11_1   10
#define   C3m_1_11   11

#define   C2110      12
#define   C2101      13
#define   C2011      14
#define   C21_10     15
#define   C2_101     16
#define   C201_1     17

#define   C4p001     18
#define   C4p010     19
#define   C4p100     20

#define   C4m001     21
#define   C4m010     22
#define   C4m100     23

#define   INV        24

#define   Mxy0       25
#define   Mx0z       26
#define   M0yz       27

#define   S3p111     28
#define   S3p1_1_1   29
#define   S3p_11_1   30
#define   S3p_1_11   31

#define   S3m111     32
#define   S3m1_1_1   33
#define   S3m_11_1   34
#define   S3m_1_11   35

#define   Mx_xz      36
#define   M_xyx      37
#define   Mxy_y      38
#define   Mxxz       39
#define   Mxyx       40
#define   Mxyy       41

#define   S4p001     42
#define   S4p010     43
#define   S4p100     44

#define   S4m001     45
#define   S4m010     46
#define   S4m100     47

/*** for hexagonal coordinate system ***/
#define   ONEh       48
#define   C3p001h    49
#define   C3m001h    50
#define   C2001h     51
#define   C6p001h    52
#define   C6m001h    53
#define   C2110h     54
#define   C2100h     55
#define   C2010h     56
#define   C21_10h    57
#define   C2120h     58
#define   C2210h     59
#define   INVh       60
#define   S3p001h    61
#define   S3m001h    62
#define   Mxy0h      63
#define   S6p001h    64
#define   S6m001h    65
#define   Mx_xzh     66
#define   Mx2xzh     67
#define   M2xxzh     68
#define   Mxxzh      69
#define   Mx0zh      70
#define   M0yzh      71

/*  numbers of space groups for each point group  */
#define   NC1_sgrp    1
#define   NCi_sgrp    1

#define   NC2_sgrp    3
#define   NCs_sgrp    4
#define   NC2h_sgrp   6

#define   ND2_sgrp    9
#define   NC2v_sgrp   22
#define   ND2h_sgrp   28

#define   NC4_sgrp    6
#define   NS4_sgrp    2
#define   NC4h_sgrp   6
#define   ND4_sgrp    10
#define   NC4v_sgrp   12
#define   ND2d_sgrp   6
#define   ND2d2_sgrp  6
#define   ND4h_sgrp   20

#define   NC3_sgrp    4
#define   NC3i_sgrp   2
#define   ND3_sgrp    3
#define   ND32_sgrp   4
#define   NC3v_sgrp   4
#define   NC3v2_sgrp  2
#define   ND3d_sgrp   2
#define   ND3d2_sgrp  4

#define   NC6_sgrp    6
#define   NC3h_sgrp   1
#define   NC6h_sgrp   2
#define   ND6_sgrp    6
#define   NC6v_sgrp   4
#define   ND3h_sgrp   2
#define   ND3h2_sgrp  2
#define   ND6h_sgrp   4

#define   NT_sgrp     5
#define   NTh_sgrp    7
#define   NO_sgrp     8
#define   NTd_sgrp    6
#define   NOh_sgrp    10

#endif
