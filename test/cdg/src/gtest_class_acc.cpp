//==================================================================================
// Module       : gtest_class_acc.cpp
// Date         : 2/1/18 (DLR)
// Description  : GeoFLOW test of class usage with OpenACC
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gtypes.h"
#include "gtvector.hpp"
#include "gtstlvec.hpp"
#include <unistd.h>
#include <gptl.h>

int main(int argc, char* argv[])
{
  std::cout << "main: enter..." << std::endl;
  const char   *serr="GeoFLOW::main: ";
  GINT       iopt;
  GINT       N=10, nrpt=1;
  GDVector            gvec1 , gvec2 , gvec3 , tmp;
//GTSTLVec<GDOUBLE>   gsvec1, gsvec2, gsvec3;

  std::cout << "main: do option parsing..." << std::endl;

  // Parse command line. ':' after char
  // option indicates that it takes an argument:
  while ((iopt = getopt(argc, argv, "i:n:h")) != -1) {
    switch (iopt) {
    case 'i': // get iteration count
        nrpt = atoi(optarg);
        break;
    case 'n': // get array/vector size 
        N = atoi(optarg);
        break;
    case 'h': // help
        std::cout << "usage: " << std::endl <<
        argv[0] << " [-h] [-n #Size] [-i #Repeat] " << std::endl;
        exit(1); 
        break;
    case ':': // missing option argument
        std::cout << argv[0] << ": option " << optopt << " requires an argument" << std::endl;
        exit(1); 
        break;
    case '?':
    default: // invalid option
        std::cout << argv[0] << ": option " << optopt << " invalid" << std::endl;
        exit(1);
        break;
    }
  }

  std::cout << "main: initializing GPTL..." << std::endl;

  // Initialize GPTL:
  GPTLsetoption(GPTL_DCMRT,1);   // L1 miss rate
  GPTLsetoption(GPTL_L2MRT,1);   // L2 miss rate
  GPTLsetoption(GPTL_L3MRT,1);   // L3 miss rate
  GPTLsetoption(GPTL_CI   ,1);   // Comput intensity
  GPTLsetoption(GPTL_LSTPDCM,1); // Load-stores per L1 miss 
  GPTLsetoption(GPTL_FPI,1);     // P ops per instruction
  GPTLinitialize();
  
  std::cout << "main: initializing vectors..." << std::endl;
  std::cout << "main: ... resizing..." << std::endl;
  gvec1.resize(N);
  gvec2.resize(N);
  gvec3.resize(N);
//gsvec1.resize(N);
//gsvec2.resize(N);
//gsvec3.resize(N);
  std::cout << "main: ... setting..." << std::endl;
  for ( GSIZET j=0; j<N; j++ )
  {
    gvec1 [j] = 2.0*j;
    gvec2 [j] = 3.0*j;
    gvec3 [j] = gvec2[j];
//  gsvec1[j] = 2.0*j;
//  gsvec2[j] = 3.0*j;
//  gsvec3[j] = gsvec2[j];
  }
  std::cout << "main: updating on device..." << std::endl;
  gvec1.updatedev();
  gvec2.updatedev();
  gvec3.updatedev();
//gsvec1.updatedev();
//gsvec2.updatedev();
//gsvec3.updatedev();
 
  std::cout << "main: enter repeat loop..." << std::endl;

  std::cout << "main: do compute gvec +=..." << std::endl;
  GPTLstart("gvec+=");
  for ( GINT j=0; j<nrpt; j++ )
  {
  #pragma acc parallel 
    gvec2 += gvec1; 
  }
  GPTLstop("gvec+=");

  std::cout << "main: update host..." << std::endl;
  gvec2.updatehost();
  std::cout << "main: do gvec print..." << std::endl;
  std::cout << serr << " gvec2=";
  for ( GSIZET j=0; j< MIN(N,20); j++ ) std::cout << gvec2[j] << " ";
  std::cout << std::endl;

/*
  std::cout << "main: do compute gsvec +=..." << std::endl;
  GPTLstart("gsvec+=");
  for ( GINT j=0; j<nrpt; j++ )
  {
  #pragma acc parallel 
    gsvec2 += gsvec1; 
  }
  GPTLstop("gsvec+=");
*/

  std::cout << "main: do compute gvec = gv1 + gv2..." << std::endl;
  GPTLstart("gvec=gv1+gv2");
  for ( GINT j=0; j<nrpt; j++ )
  {
  #pragma acc parallel 
    gvec1.add(gvec3, gvec2, 1.0, 1.0);
  }
  GPTLstop("gvec=gv1+gv2");

/*
  std::cout << "main: do compute gsvec = gsv1 + gsv2..." << std::endl;
  GPTLstart("gsvec=gsv1+gsv2");
  for ( GINT j=0; j<nrpt; j++ )
  {
  #pragma acc parallel 
    gsvec1.add(gsvec3, gsvec2, 1.0, 1.0);
  }
  GPTLstop("gsvec=gsv1+gsv2");
*/

  std::cout << "main: update host..." << std::endl;
  gvec3.updatehost();
  std::cout << "main: do print..." << std::endl;
  std::cout << serr << " gvec3=";
  for ( GSIZET j=0; j< MIN(N,20); j++ ) std::cout << gvec3[j] << " ";
  std::cout << std::endl;

  std::cout << "main: GPTL print..." << std::endl;
  GPTLpr_file("gptl.txt");
  GPTLfinalize();

  std::cout << "main: done." << std::endl;

  exit(0);

} // end of function main
