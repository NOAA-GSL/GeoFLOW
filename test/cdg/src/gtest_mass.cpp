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
#include "ggrid_factory.hpp"
#include "gmtk.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"

using namespace geoflow::tbox;
using namespace std;


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   errcode, gerrcode;
    GINT   nlat=10, nlong=20;
    GFTYPE radiusi=1, radiuso=2;
    GTVector<GINT> ne[3]; // # elements in each direction in 3d
    GC_COMM comm = GC_COMM_WORLD;

    for ( GSIZET j=0; j<3; j++ ) ne[j] = 10;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:j:k;l:o:q:r:v:h")) != -1) {
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
      case 'q': // inner radius for 2d/3d
          radiusi = atoi(optarg);
          break;
      case 'r': // outer radius for 3d
          radiuso = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-i #Elems in r] [-j #Elems in lat]  [-k #Elems in long] [-l refine level] [-q rad_inner] [-r rad_outer] " << std::endl;
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

    // Get minimal property tree:
    PropertyTree ptree; 
    GString sgrid;
    GGrid *grid;
    std::vector<GINT> pstd(GDIM);

    ptree.load_file("input.jsn");
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");

    assert(sgrid == "grid_icos" && "Must use ICOS grid for now");

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }

    // Generate grid:
    GPTLstart("gen_grid");
    grid = GGridFactory::build(ptree, gbasis, comm);
    GPTLstop("gen_grid");


    // Test some of the coord transformation methods:
    GFTYPE xlat, xlatc, xlong, xlongc;
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();
    GTPoint<GFTYPE> pdiff;
    GTVector<GTPoint<GFTYPE>> cart(1), cref(1), sph(1), tcart(1), tsph(1);

    GTVector<GFTYPE> f(grid->ndof());
    GTVector<GFTYPE> g(grid->ndof());
    GTVector<GTVector<GFTYPE>*> utmp(1);
    GTVector<GFTYPE> imult(grid->ndof());
    GMass massop(*grid);


    f = 1.0;
    GPTLstart("massop_prod");
    massop.opVec_prod(f,utmp,g);
    GPTLstop("massop_prod");
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
      ios << "#elems" << "  ";
      for ( GSIZET j=0; j<GDIM; j++ ) ios << "p" << j+1 << "  ";
      ios <<  "ndof    area_computed    area_analytic    diff " << std::endl;
    }
    itst.close();

    // Get global no elements and dof:
    GTVector<GSIZET> lsz(2), gsz(2);
    lsz[0] = grid->nelems();
    lsz[1] = grid->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm);


    GFTYPE aintegral;
    #if defined(_G_IS2D)
    aintegral = 4.0*PI*pow(radiusi,2.0);
    #elif defined(_G_IS3D)
    aintegral = 4.0*PI*pow(radiusi,3.0)/3.0;
    #endif
    std::cout << "main: integral=" << gintegral << "; area=" << aintegral << std::endl;

    ios <<   gsz[0] << "  ";
    for ( GSIZET j=0; j<GDIM; j++ ) ios << pstd[j] << "  ";
    ios << gsz[1]     << "  "
        << gintegral  << "  " 
        << aintegral  << "  "
        << fabs(gintegral-aintegral) << std::endl;
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
