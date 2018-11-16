//======================================================================================
// Name         : cff_blas.h
// Date         : 1/1/18 (DLR)
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Description  : C-header for ISO-bound cache-friendly optimized Fortran 
//                linear algebra routines
//======================================================================================

#include "gtypes.h"

extern "C" {
//======================================================================================
//   Quad-precision  versions:
//======================================================================================

#pragma acc routine vector
void qmxm     (GQUAD [], GQUAD [], GSIZET *, GSIZET *, GQUAD [], GSIZET *, GSIZET *, GINT *);
//#pragma acc routine vector
//void qmxmcf   (GQUAD [], GQUAD [], GSIZET *, GSIZET *, GQUAD [], GSIZET *, GSIZET *, GINT *);
#pragma acc routine vector
void qmxv     (GQUAD [], GQUAD [], GQUAD [], GSIZET *, GSIZET *, GINT *);
#pragma acc routine vector
void qmxDm    (GQUAD [], GQUAD [], GSIZET *, GSIZET *, GQUAD [], GSIZET *, GINT *);
#pragma acc routine vector
void qDmxm    (GQUAD [], GQUAD [], GSIZET *, GQUAD [], GSIZET *, GSIZET *, GINT *); 
#pragma acc routine vector
void qaApbB   (GQUAD [], GQUAD [], GQUAD [], GSIZET *, GSIZET *,  GQUAD *, GQUAD *, GINT *); 
#pragma acc routine vector
void qzaxpby  (GQUAD [], GQUAD [], GQUAD *, GQUAD [], GQUAD *, GSIZET *, GINT *); 
#pragma acc routine vector
void qxaxpby  (GQUAD [], GQUAD *,  GQUAD [], GQUAD *, GSIZET *, GINT *); 
#pragma acc routine vector
void qzvxvpt  (GQUAD [], GQUAD [], GQUAD [], GSIZET *, GINT *); 
#pragma acc routine vector
void qvvxvpt  (GQUAD [], GQUAD [], GSIZET *, GINT *); 
#pragma acc routine vector
void qvvxvptpv(GQUAD [], GQUAD [], GQUAD [], GQUAD *, GSIZET *, GINT *); 
#pragma acc routine vector
void qvpvxvpt (GQUAD [], GQUAD [], GQUAD [], GQUAD *, GSIZET *, GINT *); 
#pragma acc routine vector
void qdotg    (GQUAD  *, GQUAD [], GQUAD [],  GSIZET *, GINT *); 
#pragma acc routine vector
void qcopy    (GQUAD [], GQUAD [], GSIZET *, GINT *); 
#pragma acc routine vector
void qoop     (GQUAD [], GQUAD [], GSIZET [], GSIZET *, GSIZET *); 
#pragma acc routine vector
void qisassign(GQUAD [], GQUAD [], GSIZET [], GSIZET *, GSIZET *); 
#pragma acc routine vector
void qisum    (GQUAD [], GQUAD [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void qiprod   (GQUAD [], GQUAD [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void qimax    (GQUAD [], GQUAD [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void qimin    (GQUAD [], GQUAD [], GSIZET [], GSIZET *); 

//======================================================================================
//   Double-precision  versions:
//======================================================================================

#pragma acc routine vector
void dmxm     (GDOUBLE [], GDOUBLE [], GSIZET *, GSIZET *, GDOUBLE [], GSIZET *, GSIZET *, GINT *);
//#pragma acc routine vector
//void dmxmcf   (GDOUBLE [], GDOUBLE [], GSIZET *, GSIZET *, GDOUBLE [], GSIZET *, GSIZET *, GINT *);
#pragma acc routine vector
void dmxv     (GDOUBLE [], GDOUBLE [], GDOUBLE [], GSIZET *, GSIZET *, GINT *);
#pragma acc routine vector
void dmxDm    (GDOUBLE [], GDOUBLE [], GSIZET *, GSIZET *, GDOUBLE [], GSIZET *, GINT *);
#pragma acc routine vector
void dDmxm    (GDOUBLE [], GDOUBLE [], GSIZET *, GDOUBLE [], GSIZET *, GSIZET *, GINT *); 
#pragma acc routine vector
void daApbB   (GDOUBLE [], GDOUBLE [], GDOUBLE [], GSIZET *, GSIZET *,  GDOUBLE *, GDOUBLE *, GINT *); 
#pragma acc routine vector
void dzaxpby  (GDOUBLE [], GDOUBLE [], GDOUBLE *, GDOUBLE [], GDOUBLE *, GSIZET *, GINT *); 
#pragma acc routine vector
void dxaxpby  (GDOUBLE [], GDOUBLE *,  GDOUBLE [], GDOUBLE *, GSIZET *, GINT *); 
#pragma acc routine vector
void dzvxvpt  (GDOUBLE [], GDOUBLE [], GDOUBLE [], GSIZET *, GINT *); 
#pragma acc routine vector
void dvvxvpt  (GDOUBLE [], GDOUBLE [], GSIZET *, GINT *); 
#pragma acc routine vector
void dvvxvptpv(GDOUBLE [], GDOUBLE [], GDOUBLE [], GDOUBLE *, GSIZET *, GINT *); 
#pragma acc routine vector
void dvpvxvpt (GDOUBLE [], GDOUBLE [], GDOUBLE [], GDOUBLE *, GSIZET *, GINT *); 
#pragma acc routine vector
void ddotg    (GDOUBLE  *, GDOUBLE [], GDOUBLE [],  GSIZET *, GINT *); 
#pragma acc routine vector
void dcopy    (GDOUBLE [], GDOUBLE [], GSIZET *, GINT *); 
#pragma acc routine vector
void doop     (GDOUBLE [], GDOUBLE [], GSIZET [], GSIZET *, GSIZET *); 
#pragma acc routine vector
void disassign(GDOUBLE [], GDOUBLE [], GSIZET [], GSIZET *, GSIZET *); 
#pragma acc routine vector
void disum    (GDOUBLE [], GDOUBLE [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void diprod   (GDOUBLE [], GDOUBLE [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void dimax    (GDOUBLE [], GDOUBLE [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void dimin    (GDOUBLE [], GDOUBLE [], GSIZET [], GSIZET *); 

//======================================================================================
//   Float/single-precision  versions:
//======================================================================================


#pragma acc routine vector
void fmxm     (GFLOAT [], GFLOAT [], GSIZET *, GSIZET *, GFLOAT [], GSIZET *, GSIZET *, GINT *);
//#pragma acc routine vector
//void fmxmcf   (GFLOAT [], GFLOAT [], GSIZET *, GSIZET *, GFLOAT [], GSIZET *, GSIZET *, GINT *);
#pragma acc routine vector
void fmxv     (GFLOAT [], GFLOAT [], GFLOAT [], GSIZET *, GSIZET *, GINT *);
#pragma acc routine vector
void fmxDm    (GFLOAT [], GFLOAT [], GSIZET *, GSIZET *, GFLOAT [], GSIZET *, GINT *);
#pragma acc routine vector
void fDmxm    (GFLOAT [], GFLOAT [], GSIZET *, GFLOAT [], GSIZET *, GSIZET *, GINT *); 
#pragma acc routine vector
void faApbB   (GFLOAT [], GFLOAT [], GFLOAT [], GSIZET *, GSIZET *,  GFLOAT *, GFLOAT *, GINT *); 
#pragma acc routine vector
void fzaxpby  (GFLOAT [], GFLOAT [], GFLOAT *, GFLOAT [], GFLOAT *, GSIZET *, GINT *); 
#pragma acc routine vector
void fxaxpby  (GFLOAT [], GFLOAT *,  GFLOAT [], GFLOAT *, GSIZET *, GINT *); 
#pragma acc routine vector
void fzvxvpt  (GFLOAT [], GFLOAT [], GFLOAT [], GSIZET *, GINT *); 
#pragma acc routine vector
void fvvxvpt  (GFLOAT [], GFLOAT [], GSIZET *, GINT *); 
#pragma acc routine vector
void fvvxvptpv(GFLOAT [], GFLOAT [], GFLOAT [], GFLOAT *, GSIZET *, GINT *); 
#pragma acc routine vector
void fvpvxvpt (GFLOAT [], GFLOAT [], GFLOAT [], GFLOAT *, GSIZET *, GINT *); 
#pragma acc routine vector
void fdotg    (GFLOAT  *, GFLOAT [], GFLOAT [],  GSIZET *, GINT *); 
#pragma acc routine vector
void fcopy    (GFLOAT [], GFLOAT [], GSIZET *, GINT *); 
#pragma acc routine vector
void foop     (GFLOAT [], GFLOAT [], GSIZET [], GSIZET *, GSIZET *); 
#pragma acc routine vector
void fisassign(GFLOAT [], GFLOAT [], GSIZET [], GSIZET *, GSIZET *); 
#pragma acc routine vector
void fisum    (GFLOAT [], GFLOAT [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void fiprod   (GFLOAT [], GFLOAT [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void fimax    (GFLOAT [], GFLOAT [], GSIZET [], GSIZET *); 
#pragma acc routine vector
void fimin    (GFLOAT [], GFLOAT [], GSIZET [], GSIZET *); 

}
