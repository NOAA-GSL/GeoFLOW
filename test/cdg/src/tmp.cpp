//==================================================================================
// Module       : gtest_mass.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW tests of GMass operator
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include <random>
#include "gcomm.hpp"
#include "gllbasis.hpp"
#include "ggrid_icos.hpp"
#include "gmass.hpp"


int main(int argc, char **argv)
{
    // Initialize comm:
    GComm::InitComm(&argc, &argv);

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   errcode=0, gerrcode;
    GINT   nlat=10, nlong=20;
    GINT   np=1;    // elem 'order'
    GFTYPE radiusi=1, radiuso=2;
    GTVector<GINT> ne[3]; // # elements in each direction in 3d
    GC_COMM comm = GC_COMM_WORLD;

    for ( GSIZET j=0; j<3; j++ ) ne[j] = 10;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:j:k;l:p:q:r:h")) != -1) {
      switch (iopt) {
      case 'i': // get # elements in r
          ne[0] = atoi(optarg);
          break;
      case 'j': // get # elements in lat
          ne[1] = atoi(optarg);
          break;
      case 'k': // get # elements in long
          ne[2] = atoi(optarg);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
          break;
      case 'q': // inner radius for 2d/3d
          radiusi = atoi(optarg);
          break;
      case 'r': // outer radius for 3d
          radiuso = atoi(optarg);
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

#endif

    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();


    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) gbasis[k] = new GLLBasis<GCTYPE,GFTYPE>(np);
    
    GGrid grid;

    #if defined(_G_IS2D)
    GGridIcos gen_icos(radiusi, ilevel, gbasis, nprocs);
    #elif defined(_G_IS3D)
    GGridIcos gen_icos(radiusi, radiuso, ne, gbasis, nprocs);
    #endif


    // Generate grid:
    gen_icos.do_grid(grid, myrank);
    grid.print("grid.txt",TRUE); // print internal dof too

    // Compute function, f=1 over grid:
    GTVector<GFTYPE> f(grid.ndof());
    GTVector<GFTYPE> g(grid.ndof());
    f = 1.0;

    GMass Mass(grid);
    Mass.opVec_prod(f,g);

    #if defined(_G_IS2D)
    std::cout << "integral(f)= " << g.sum() << std::endl;
    std::cout << "Area= " << 4*PI*pow(radiusi,2)<< std::endl;
    #elif defined(_G_IS3D)
    std::cout << "integral(f)= " << g.sum() << std::endl;
    std::cout << "Volume= " << 4*PI*pow(radiusi,3)/3.0 << std::endl;
    #endif
    

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
