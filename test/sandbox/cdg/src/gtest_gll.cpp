//==================================================================================
// Module       : gtest_ggrid.cpp
// Date         : 7/24/18 (DLR)
// Description  : GeoFLOW tests of GL, GLL classes
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <vector>
#include "gptl.h"
#include "gcomm.hpp"
#include "gllbasis.hpp"


int main(int argc, char **argv)
{
    // Initialize comm:
    GComm::InitComm(&argc, &argv);

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   errcode=0, gerrcode;
    GINT   np=1;    // elem 'order'
    GC_COMM comm = GC_COMM_WORLD;


    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "p:h")) != -1) {
      switch (iopt) {
      case 'p': // get nodal exp order
          np = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-i #Elems in r] [-j #Elems in lat]  [-k #Elems in long] [-l refine level] [-p expansion order] [-q rad_inner] [-r rad_outer]  " << std::endl;
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

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();


    // Create basis:
    GLLBasis<GCTYPE,GFTYPE> gllbasis(np);
    
    GTVector<GFTYPE> gllnodes(np+1);
    GTVector<GFTYPE> gllwghts(np+1);
    GTMatrix<GFTYPE> deriv(np+1,np+1);

    gllbasis.getXiNodes(gllnodes);
    gllbasis.getWeights(gllwghts);
    gllbasis.evalDBasis(gllnodes, deriv);

    std::cout << std::setw(18) << std::setprecision(16) << "gllnodes =" << gllnodes << std::endl; 
    std::cout << std::setw(18) << std::setprecision(16) << "gllwghts =" << gllwghts << std::endl << std::endl; 
    std::cout << std::setw(18) << std::setprecision(16) << "deriv=" << deriv<< std::endl; 

    // Accumulate errors:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm);

 
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
