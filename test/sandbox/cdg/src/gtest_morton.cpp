//==================================================================================
// Module       : gtest_morton.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW test of Morton indices object
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#if defined(GEOFLOW_USE_GPTL)
#include "gptl.h"
#endif
#include <random>
#include "gcomm.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   errcode=0, gerrcode;
    GINT   np=1;    // elem 'order'
    GFTYPE scalelen;
    GC_COMM comm = GC_COMM_WORLD;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "p:r:h")) != -1) {
      switch (iopt) {
      case 'p': // get nodal exp order
          np = atoi(optarg);
          break;
      case 'r': // scale length
          scalelen = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] -p expansion order] [-r len_scale ]  " << std::endl;
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

    // Initialize comm:
    GComm::InitComm(&argc, &argv);


    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

#if defined(GEOFLOW_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();
#endif


    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) gbasis[k] = new GLLBasis<GCTYPE,GFTYPE>(np);
    

    // Test some of the coord transformation methods:
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();
    GTPoint<GFTYPE> pdiff;

    GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
    GTPoint<GFTYPE> dX, porigin, P0;

    
    P0.x1= -scalelen ; P0.x2=-scalelen; P0.x3=-scalelen;
    dX   = 1.0e-3;
    gmorton.setIntegralLen(P0,dX);

    // NOTE: only need to find indices for boundary
    //       nodes, and use elem bdy indirection in GGFX/DoOp, 
    //       but for testing, we do:
    GElemList *gelems = &grid.elems();
    GTVector<GTVector<GFTYPE>> *xnodes;
    for ( GSIZET i=0; i<gelems->size(); i++ ) {
      glob_indices.range((*gelems)[i]->igbeg(),(*gelems)[i]->igend()); // restrict to this range
      xnodes = &(*gelems)[i]->xNodes();
      gmorton.key(glob_indices, *xnodes);
std::cout << "main: xNodes[" << i << "]=" << *xnodes << std::endl;
std::cout << "main: glob_indices[" << i << "]=" << glob_indices << std::endl;
    }
    glob_indices.range(0,grid.ndof()-1); // must reset to full range

    // Compute connectivity map:
    bret = ggfx.Init(glob_indices);
    errcode += !bret ? 2 : 0;

    imult = 1.0;
    bret = ggfx.DoOp<GFTYPE>(imult, GGFX_OP_SUM);
    errcode += !bret ? 3 : 0;
    
    for ( GSIZET j=0; j<imult.size(); j++ ) imult[j] = 1.0/imult[j];
    std::cout << "main: imult=" << imult << std::endl;


    f = 1.0;
    massop.opVec_prod(f,g);
    std::cout << "main: mass_prod=" << g << std::endl;
    std::cout << "main: mass_prod_sum=" << g.sum() << std::endl;
  
    // Multiply f by inverse multiplicity:
    g.pointProd(imult);

    #if defined(_G_IS2D)
    std::cout << "main: integral=" << g.sum() << "; area=" << 4.0*PI*pow(radiusi,2.0) << std::endl;
    #elif defined(_G_IS3D)
    std::cout << "main: integral=" << g.sum() << "; volume=" << 4.0*PI*pow(radiusi,3.0)/3.0 << std::endl;
    #endif

    // Accumulate errors:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm);

 
#if defined(GEOFLOW_USE_GPTL)
    GPTLpr_file("timing.txt");
    GPTLfinalize();
#endif


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
