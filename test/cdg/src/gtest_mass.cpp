//==================================================================================
// Module       : gtest_mass.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW test of GMass classes
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
#include "ggfx.hpp"
#include "gmorton_keygen.hpp"


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   errcode, gerrcode;
    GINT   nlat=10, nlong=20;
    GINT   np=1;    // elem 'order'
    GFTYPE radiusi=1, radiuso=2;
    GTVector<GINT> ne[3]; // # elements in each direction in 3d
    GC_COMM comm = GC_COMM_WORLD;

    for ( GSIZET j=0; j<3; j++ ) ne[j] = 10;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:j:k;l:o:p:q:r:v:h")) != -1) {
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
          argv[0] << " [-h] [-i #Elems in r] [-j #Elems in lat]  [-k #Elems in long] [-l refine level] -p expansion order] [-q rad_inner] [-r rad_outer] " << std::endl;
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

    errcode = 0;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);


    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();



    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(np);
std::cout << "main: gbasis [" << k << "]_order=" << gbasis [k]->getOrder() << std::endl;
    }
    
    GGrid grid;

    GPTLstart("gen_base_grid");

    #if defined(_G_IS2D)
    GGridIcos gen_icos(radiusi, ilevel, gbasis, nprocs);
    #elif defined(_G_IS3D)
    GGridIcos gen_icos(radiusi, radiuso, ne, gbasis, nprocs);
    #endif

    GPTLstop("gen_base_grid");

    // Test some of the coord transformation methods:
    GFTYPE xlat, xlatc, xlong, xlongc;
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();
    GTPoint<GFTYPE> pdiff;
    GTVector<GTPoint<GFTYPE>> cart(1), cref(1), sph(1), tcart(1), tsph(1);



    // Generate grid:
    gen_icos.do_grid(grid, myrank);

    GFTYPE gminsep, minsep = grid.minsep(); 
    GFTYPE gmaxsep, maxsep = grid.maxsep(); 
    GComm::Allreduce(&minsep, &gminsep, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm);
    GComm::Allreduce(&maxsep, &gmaxsep, 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
    std::cout << "main: min grid sep=" << gminsep << std::endl;
    std::cout << "main: max grid sep=" << gmaxsep << std::endl;

    GTVector<GFTYPE> f(grid.ndof());
    GTVector<GFTYPE> g(grid.ndof());
    GTVector<GFTYPE> imult(grid.ndof());
    GMass massop(grid);

#if 0
    // Generate interface node indices:
    GTVector<GNODEID> glob_indices(grid.ndof());
    GGFX ggfx;

    GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
    GTPoint<GFTYPE> dX, porigin, P0;

    
    P0.x1= -radiusi; P0.x2=-radiusi; P0.x3=-radiusi;
//  P1.x1=  radiusi; P1.x2= radiusi; P1.x3= radiusi;
    dX   = 1.0e-5*gminsep;
    gmorton.setIntegralLen(P0,dX);

    GPTLstart("gen_morton_indices");
    // NOTE: only need to find indices for boundary
    //       nodes, and use elem bdy indirection in GGFX/DoOp, 
    //       but for testing, we do:
    gelems = &grid.elems();
    GTVector<GTVector<GFTYPE>> *xnodes;
    for ( GSIZET i=0; i<gelems->size(); i++ ) {
      glob_indices.range((*gelems)[i]->igbeg(),(*gelems)[i]->igend()); // restrict to this range
      xnodes = &(*gelems)[i]->xNodes();
      gmorton.key(glob_indices, *xnodes);
#if 0
std::cout << "main: xNodes[" << i << "]=" << *xnodes << std::endl;
std::cout << "main: glob_indices[" << i << "]=" << glob_indices << std::endl;
#endif
    }
    glob_indices.range(0,grid.ndof()-1); // must reset to full range
    GPTLstop("gen_morton_indices");

    GPTLstart("ggfx_init");
    // Compute connectivity map:
    bret = ggfx.Init(glob_indices);
    errcode += !bret ? 2 : 0;
    GPTLstop("ggfx_init");

    GPTLstart("ggfx_doop");
    imult = 1.0;
    bret = ggfx.DoOp<GFTYPE>(imult, GGFX_OP_SUM);
    errcode += !bret ? 3 : 0;
    GPTLstop("ggfx_doop");
    
    for ( GSIZET j=0; j<imult.size(); j++ ) imult[j] = 1.0/imult[j];
#if 0
    std::cout << "main: imult=" << imult << std::endl;
#endif
#endif

    massop.init();

    f = 1.0;
    GPTLstart("massop_prod");
    massop.opVec_prod(f,g);
    GPTLstop("massop_prod");
#if 0
    std::cout << "main: mass_prod=" << g << std::endl;
#endif
    std::cout << "main: mass_prod_sum=" << g.sum() << std::endl;
  
#if 0
    // Multiply f by inverse multiplicity:
    g.pointProd(imult);
#endif

    GFTYPE integral=g.sum();
    GFTYPE gintegral;
    GComm::Allreduce(&integral, &gintegral, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm);

    std::ifstream itst;
    std::ofstream ios;
    itst.open("mass_err.txt");
    ios.open("mass_err.txt",std::ios_base::app);

    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    #if defined(_G_IS2D)
      ios << "# order_xy   level  area_computed  area_analytic   diff " << std::endl;
    #elif defined(_G_IS3D)
      ios << "# order_xyz  level  area_computed  area_analytic   diff " << std::endl;
    #endif
    }
    itst.close();

    GFTYPE aintegral;
    #if defined(_G_IS2D)
    aintegral = 4.0*PI*pow(radiusi,2.0);
    std::cout << "main: integral=" << gintegral << "; area=" << aintegral << std::endl;
    ios << np  << "  "  << "  " << ilevel 
        << "  " << gintegral << "  " << aintegral << "  " << fabs(gintegral-aintegral) << std::endl;
    #elif defined(_G_IS3D)
    aintegral = 4.0*PI*pow(radiusi,3.0)/3.0;
    std::cout << "main: integral=" << gintegral << "; volume=" << aintegral << std::endl;
    ios << np  << " " <<  ilevel 
        << " " << gintegral << " " << aintegral << " " << fabs(gintegral-aintegral) << std::endl;
    #endif
    ios.close();
    

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
