//==================================================================================
// Module       : gtest_gmm.cpp
// Date         : 2/28/18 (DLR)
// Description  : GeoFLOW test of GMM, memory manager
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gmm.h"
#include "gtvector.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
  GString    serr = "GeoFLOW::main: ";
  GFLOAT     *gf[10];
  GSIZET      n3, n7, nn,nsz;
  GDOUBLE    *gd;
  GTVector<GDOUBLE> *gvec;

  GMM::bverbose(0);

  nsz = 10*sizeof(GFLOAT); // Initial size of GFLOAT pointer array
  for ( GSIZET j=0; j<10; j++ ) {
    std::cout << "main: allocating region  GFLOAT[" << j << "]:" << std::endl;
    nn = 10*j+3;
    if ( j == 3 ) n3 = nn;
    if ( j == 7 ) n7 = nn;
    gf[j] = new GFLOAT [nn];
    nsz += nn*sizeof(GFLOAT);
    std::cout << "main: region  GFLOAT[" << j << "]=" << gf[j] << std::endl;
  }
  std::cout << "main: after GFLOAT new, bytesused=" << GMM::bytesused() << " expected=" << nsz << std::endl;

  nn    = n3*sizeof(GFLOAT);
  std::cout << "main: deleting GFLOAT[" << 3 << "]:" << std::endl;
  delete [] gf[3];
  gf[3] = NULLPTR;
  nsz  -= nn;
  std::cout << "main: after GFLOAT[3] delete, bytesused=" << GMM::bytesused() << " expected=" << nsz << std::endl;
  
  std::cout << "main: allocating region GDOUBLE:" << std::endl;
  gd   = new GDOUBLE[10];
  nn   = 10*sizeof(GDOUBLE);
  std::cout << "main: region  GDOUBEL=" << gd << std::endl;
  nsz += nn;
  std::cout << "main: after GDOUBLE new , bytesused=" << GMM::bytesused() << " expected=" << nsz << std::endl;

  nn    = n7*sizeof(GFLOAT);
  std::cout << "main: deleting GFLOAT[" << 7 << "]:" << std::endl;
  delete [] gf[7];
  gf[7] = NULLPTR;
  nsz  -= nn;
  std::cout << "main: after GFLOAT[7] delete, bytesused=" << GMM::bytesused() << " expected=" << nsz << std::endl;
  

  std::cout << "main: allocating region GTVector<GDOUBLE>: " << std::endl;
  nn = GMM::bytesused();
  gvec = new GTVector<GDOUBLE>(10);
  *gvec = 13.0e-8;
  std::cout << "main: region  GTVector<GDOUBEL>=" << *gvec << std::endl;
  std::cout << "main: after GTVector<GDOUBLE> new , bytesused=" << GMM::bytesused() << std::endl;
  std::cout << "main: size of GTVector<GDOUBLE>=" << GMM::bytesused()-nn << std::endl;

  std::cout << std::endl;
  for ( GINT j=0; j<10; j++ ) {
    if ( gf[j] != NULLPTR ) {
      std::cout << "main: deleting GFLOAT[" << j << "]:" << std::endl;
      delete [] gf[j];
    }
  }

  std::cout << "main: deleting GDOUBLE:" << std::endl;
  delete [] gd;

  std::cout << "main: deleting GTVector<GDOUBLE>:" << std::endl;
  delete gvec;

  GMM::prstats("gmm_stats.txt");
  GMM::prstats(NULLPTR);

  exit(0);
} // end of function main


