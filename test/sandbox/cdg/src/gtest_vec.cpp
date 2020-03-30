//==================================================================================
// Module       : gtest_vec.cpp
// Date         : 7/24/18 (DLR)
// Description  : GeoFLOW test of misc vector or matrix functionality (non-BLAS)
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "gptl.h"
#include <random>
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gtpoint.hpp"


int main(int argc, char **argv)
{
    // Initialize comm:
    GComm::InitComm(&argc, &argv);

    GString serr ="main: ";
    GBOOL  bret;
    GINT   nrpt=1;
    GINT   iopt;
    GINT   errcode=0, gerrcode;
    GSIZET ne=4;    
    GSIZET np=1;    // elem 'order'
    GC_COMM  comm = GC_COMM_WORLD;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:n:p:h")) != -1) {
      switch (iopt) {
      case 'i': // get iteration count
          nrpt = atoi(optarg);
          break;
      case 'p': // get array/vector size 
          np = atoi(optarg);
          break;
      case 'n': // # elems per task
          ne = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-n #Elems per task] [-p expansion order] [-i #Repeat] " << std::endl;
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

#endif

    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();

    // Check push_back, reserve:
    GTPoint<GFTYPE>   pt;
    GTVector<GFLOAT>  fvec;
    GTVector<GTPoint<GFTYPE>>  gpush;      // check push_back with point
    GTVector<GTPoint<GFTYPE>>  gassign(ne);
    GTVector<GTVector<GFLOAT>> gpushv     ; // check push_back with vector

    gpush .reserve(2);
    gpushv.reserve(2);
    for ( GSIZET k=0; k<ne; k++) {
      pt.x1 = static_cast<GFTYPE>(k);
      pt.x2 = static_cast<GFTYPE>(2*k);
      pt.x3 = -1;
      gassign[k] = pt;
      gpush.push_back(pt);

      fvec.resize(2*k+1);
      fvec = static_cast<GFLOAT>(2*k+1);
      std::cout << "main: fvec[" << k << "]=" << fvec << std::endl;
      gpushv.push_back(fvec);
    }

    GPP(comm, "main: gpush (pre-reserve)  =" << gpush);
    GPP(comm, "main: gassign=" << gassign);
    GPP(comm, "main: gpushv =" << gpushv);

    std::cout << "main: print individual elements of gpushv: " << std::endl;
    for ( GSIZET k=0; k<ne; k++) {
      for ( GSIZET j=0; j<gpushv[k].size(); j++) {
        std::cout << "gpushv[" << k << "][" << j << "]=" << gpushv[k][j] << " ";
      }
    }
    std::cout << std::endl;;

    gpush.reserve(2*ne);
    GPP(comm, "main: gpush (post-reserve) =" << gpush);


    for ( GSIZET k=ne; k<2*ne; k++) {
      pt.x1 = static_cast<GFTYPE>(2*k);
      pt.x2 = static_cast<GFTYPE>(6*k);
      pt.x3 = -1;
      gpush.push_back(pt);
    }
    GPP(comm, "main: gpush (post-post-resize)  =" << gpush);

    // Check matrix as a container for vectors:
    GTMatrix<GTVector<GFLOAT>> vmat       ; // matrix's container access

    ne = 2;
    vmat.resize(ne,ne);
    for ( GSIZET k=0; k<ne; k++) {
      for ( GSIZET j=0; j<ne; j++) {
        fvec.resize(2*k+1);
        fvec = static_cast<GFLOAT>(2*k+1 + j);
        std::cout << "main: fvec[" << k << "]=" << fvec << std::endl;
        vmat(j,k) = fvec;
      }
    }
    GPP(comm, "main: vmat=" << vmat);

    GTVector<GINT> v1(10);
    GTVector<GINT> *v2;

    v1 = 13;
    v2 = &v1;
    std::cout << "v2=" << *v2 << std::endl;
    for ( GSIZET k=0; k<15; k++ ) v1.push_back(k);
    std::cout << "v1=" << v1 << std::endl;

    std::default_random_engine rgen;
    std::uniform_int_distribution<GINT> dist(1,100);

    v1.resize(10);
    GTVector<GSIZET> isort;
    for ( GSIZET k=0; k<v1.size(); k++ ) v1[k] = dist(rgen);
    std::cout << "v1=" << v1 << std::endl;
    v1.sortincreasing(isort);
    std::cout << "isort.increasing=" << isort << std::endl; 
    std::cout << "v1(isort)= ";
    for ( GSIZET k=0; k<v1.size()-1; k++ )
      std::cout << v1[isort[k]] << " ";
    std::cout << v1[isort[v1.size()-1]] << std::endl;

    v1.sortincreasing();
    std::cout << "v1.sortincreasing=" << v1 << std::endl; 

    for ( GSIZET k=0; k<v1.size(); k++ ) v1[k] = dist(rgen);
    std::cout << "v1=" << v1 << std::endl;
    v1.sortdecreasing(isort);
    std::cout << "isort.decreasing=" << isort << std::endl; 
    std::cout << "v1(isort)= ";
    for ( GSIZET k=0; k<v1.size()-1; k++ )
      std::cout << v1[isort[k]] << " ";
    std::cout << v1[isort[v1.size()-1]] << std::endl;

    v1.sortdecreasing();
    std::cout << "v1.sortdecreasing=" << v1 << std::endl; 
 
    GPTLpr_file("timing.txt");
    GPTLfinalize();


term:
    if ( gerrcode != 0 ) {
      GPP(comm,serr << " Error: code=" << errcode);
    }
    else {
      GPP(comm,serr << "     Success!");
    }

    GComm::TermComm();
    return(gerrcode);
} // end, main
