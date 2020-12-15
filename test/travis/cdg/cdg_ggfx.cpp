//==================================================================================
// Module       : gtest_ggfx.cpp
// Date         : 5/30/18 (DLR)
// Description  : GeoFLOW test of GGFX operator
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <vector>
#if defined(GEOFLOW_USE_GPTL)
#include "gptl.h"
#endif
#include "gexec.h"
#include "ggfx.hpp"

int main(int argc, char **argv)
{
    // Initialize comm:
    GComm::InitComm(&argc, &argv);

    GString serr ="main: ";
    GBOOL  bret;
    GINT   nrpt=1;
    GINT   iopt;
    GINT   errcode=0, gerrcode;
    GSIZET ne=16; // no. elements per task (# 'y' elems)
    GSIZET np=4;    // elem 'order'
    GFLOAT eps = std::numeric_limits<GFLOAT>::epsilon();
    GC_COMM comm = GC_COMM_WORLD;
    GGFX<GFLOAT>   ggfx;


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

    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    GPP(comm, "*********************************************************");
    GPP(comm,serr << ": no. elems per task: " << ne << " np: " << np << " nprocs: " << nprocs );
    GPP(comm, "*********************************************************" << std::endl );

#if defined(GEOFLOW_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();
#endif


    // Build 'grid':

/*
   . ne
   .
   .               ...
    |         |         |         |         |
 -- o -- o -- o -- o -- o -- o -- o -- o -- o
    |    |    |    |    |    |    |    |    |
 -- o -- o -- o -- o -- o -- o -- o -- o -- o
    |    |    |    |    |    |    |    |    |
 -- o -- o -- o -- o -- o -- o -- o -- o -- o
    |   T0    |   T1    |  T2     |   T3    |
   
    0    1    2    3    4    5    6    7    8


NOTE: global ids are labeled starting from left on bottom-most
      row, advancing along the row, continuing to the next 
      row, etc. T0, T1, ... represent task ids.
*/
    GTVector<GNODEID>    glob_indices;
    glob_indices.resize(ne*(np+1)*(np+1));
    glob_indices = 0;

    GSIZET  n;
    GSIZET  jstart;

    GSIZET rowlen = nprocs*np + 1; // # global dof in a row of nodes
    GSIZET collen = ne    *np + 1; // # global dof in a col of nodes
    GSIZET c, r, icol, irow;
    

    // Label each node point of local elements with
    // their global index values:
    n = 0;
    for ( GSIZET k=0; k<ne; k++) {
      for ( GSIZET i=0; i<np+1; i++) {
        irow   = k*np + i; // global row id
        jstart = irow*rowlen + myrank*np ; // global start of col index
        for ( GSIZET j=0; j<np+1; j++) {
          glob_indices[n] = jstart + j;
          n++;
        }
      }
    }

    GPP(comm,serr << " glob_indices=" << glob_indices);
 
    // Set field, and analytic field buffers:
    GTVector<GFLOAT>  u (glob_indices.size()); // trial field
    GTVector<GFLOAT>  ua(glob_indices.size()); // analytic solution
    GTVector<GFLOAT>  utmp(glob_indices.size()); // tmp buffer


    // Set analytic solution of GGFX_OP_SUM:
    ua = 1; 

    GPP(comm,serr << " ..........................start analytic loop:");
    GBOOL isLeftTask = nprocs > 1 && myrank > 0 ? TRUE: FALSE;
    GBOOL isRghtTask = nprocs > 1 && myrank < nprocs - 1 ? TRUE: FALSE;
    GBOOL isTopElem, isBotElem;
    for ( GSIZET k=0,n=0; k<ne; k++) { // loop over # elements on this task
      isTopElem = ne > 1 && k < ne-1 ? TRUE : FALSE;
      isBotElem = ne > 1 && k > 0 ? TRUE : FALSE;
      for ( GSIZET i=0; i<np+1; i++) {
        irow  = k*np + i;    // global row id
        jstart = irow*rowlen + myrank*np ; // global start of col index
        for ( GSIZET j=0; j<np+1; j++) {
          icol = jstart + j; // global col id
#if 1
          if ( isTopElem  && i == np ) ua[n] = 2;
          if ( isBotElem  && i ==  0 ) ua[n] = 2;
          if ( isLeftTask && j ==  0 ) ua[n] = 2;
          if ( isRghtTask && j == np ) ua[n] = 2;
          if ( isLeftTask && isTopElem && i == np && j == 0  ) ua[n] = 4;
          if ( isLeftTask && isBotElem && i ==  0 && j == 0  ) ua[n] = 4;
          if ( isRghtTask && isTopElem && i == np && j == np ) ua[n] = 4;
          if ( isRghtTask && isBotElem && i ==  0 && j == np ) ua[n] = 4;
#else
          ua[n] = irow+icol;
          r = irow; c = icol;
          if ( isTopElem  && i == np ) ua[n] = irow+icol + irow+icol+1;
          if ( isBotElem  && i ==  0 ) ua[n] = irow+icol + irow+icol-1;
          if ( isLeftTask && j ==  0 ) ua[n] = irow+icol + irow-1+icol;
          if ( isRghtTask && j == np ) ua[n] = irow+icol + irow+1+icol;
          if ( isLeftTask && isTopElem && i == np && j == 0  ) ua[n] = r+c + r+c+1 + r-1+c+1 + r-1+c;
          if ( isLeftTask && isBotElem && i ==  0 && j == 0  ) ua[n] = r+c + r+1+c + r+1+c+1 + r+c+1;
          if ( isRghtTask && isTopElem && i == np && j == np ) ua[n] = r+c + r-1+c + r-1+c-1 + r+c-1;
          if ( isRghtTask && isBotElem && i ==  0 && j == np ) ua[n] = r+c + r+c-1 + r+1+c-1 + r+1+c;
#endif
          n++;
        }
      }
    }
    GPP(comm,serr << " ..........................end analytic loop:");

    // Print ua to task 0:

    GPP(comm, serr << " doing GGFX::Init...");
    u = 1;
#if defined(GEOFLOW_USE_GPTL)
    GPTLstart("ggfx_init");
#endif
    for ( GSIZET k=0; k<nrpt && errcode==0; k++ ) { 
      bret = ggfx.init(glob_indices);
      if ( !bret ) errcode = 2;
    }
    GComm::Synch();
    GPP(comm, serr << " GGFX::Init done.");
#if defined(GEOFLOW_USE_GPTL)
    GPTLstop("ggfx_init");
#endif

    GPP(comm,serr << " doing GGFX::doOp...");
#if defined(GEOFLOW_USE_GPTL)
    GPTLstart("ggfx_doop");
#endif
    for ( GSIZET k=0; k<nrpt && errcode==0; k++ ) { 
      bret = ggfx.doOp(u, GGFX_OP_SUM);
      if ( !bret ) errcode = 3;
    }
#if defined(GEOFLOW_USE_GPTL)
    GPTLstop("ggfx_doop");
#endif
    if ( errcode == 0 ) GPP(comm,serr << " GGFX::doOp done.");


    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm);


    // Check against truth solution:
    // Print error code to task 0
    GFLOAT           ldel, gdel;
    GTVector<GFLOAT> del(ua.size());
    if ( gerrcode == 0 ) {
      del = ua - u; 
      ldel = del.max();
      GComm::Allreduce(&ldel, &gdel, 1, T2GCDatatype<GFLOAT>() , GC_OP_MAX, comm);
      if ( gdel > eps ) {
        GPP(comm,serr << "  ua=" << ua);
        GPP(comm,serr << "  ua=" << ua);
        GPP(comm,serr << "  u =" << u );
        GPP(comm,serr << " del=" << del);
        GPP(comm,serr << "gdel=" << gdel);
        GPP(comm,"main: Solution error in GGFX test");
        errcode = 4;
      }
    }

#if defined(GEOFLOW_USE_GPTL)
    GPTLpr_file("timing.txt");
    GPTLfinalize();
#endif

    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm);

term:
    if ( myrank == 0 ) {
      if ( gerrcode != 0 ) {
        cout << serr << "  Error: code=" << gerrcode << endl;
      }
      else {
        cout << serr << "     Success!" << endl;
      }
    }

    GComm::TermComm();
    return(gerrcode);
} // end, main
