//==================================================================================
// Module       : gtest_gen_type.cpp
// Date         : 2/1/18 (DLR)
// Description  : GeoFLOW test of C-calling-Fortran binding
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gtypes.h"
#include <iostream>

#if defined(__cplusplus)
extern "C" {
void dofprod(GFLOAT  [], GFLOAT  [], GFLOAT  [], GINT *);
void dodprod(GDOUBLE [], GDOUBLE [], GDOUBLE [], GINT *);
}
#else
extern void dofprod(GFLOAT  [], GFLOAT  [], GFLOAT  [], GINT *);
extern void dodprod(GDOUBLE [], GDOUBLE [], GDOUBLE [], GINT *);
#endif

#define N 10

int main(int argc, char* argv[])
{
  GString    serr = "GeoFLOW::main: ";
  GINT       n=N;
  GFLOAT     gf[N], fa[N], fb[N];
  GDOUBLE    gd[N], da[N], db[N];

  for ( GINT j=0; j<n; j++ ) {
    fa[j] = 10.0; 
    fb[j] = 20.0+j; 
    da[j] = 30.0; 
    db[j] = 40.0+j; 
  }

  dofprod(gf, fa, fb, &n);
  dodprod(gd, da, db, &n);

  std::cout << "float array check:" << std::endl;
  for ( GINT j=0; j<n; j++ ) {
    std::cout << fa[j] << "*" << fb[j] << "=" << gf[j] << std::endl;
  }

  std::cout << "double array check:" << std::endl;
  for ( GINT j=0; j<n; j++ ) {
    std::cout << da[j] << "*" << db[j] << "=" << gd[j] << std::endl;
  }

} // end of function main


