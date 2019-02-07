//==================================================================================
// Module       : gtest_burgers.cpp
// Date         : 12/24/18 (DLR)
// Description  : GeoFLOW test of Burgers PDE
// Copyright    : Copyright 2018-19. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggrid.hpp"
#include "ggrid_factory.hpp"
#include "ggfx.hpp"
#include "gmorton_keygen.hpp"
#include "gburgers.hpp"
#include "pdeint/equation_base.hpp"

using namespace geoflow::pdeint;


template<
typename StateType,
typename ValueType = GFTYPE,
typename DerivType = StateType,
typename TimeType  = ValueType,
typename JacoType  = NULLPTR
typename SizeType  = GSIZET
>
struct EquationTypes {
        using State      = StateType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};

void init_grid() {
};

void compute_analytic(GGrid &grid, Time &t, State &uanalyt);

int main(int argc, char **argv)
{

    // Get types used for equations and solver
    using MyTypes = EquationTypes<GTVector<GTVectorGFTYPE>*>>; // Define types used
    using EqnBase = EquationBase<MyTypes>;                     // Equation Base Type
    using EqnImpl = GBurgers<MyTypes>;                         // Equation Implementa
    using ObsBase = ObserverBase<EqnBase>;                     // Observer Base Type
    using ObsImpl = NullObserver<EqnBase>;                     // Observer Implementation Type
    using IntImpl = Integrator<EqnBase>;                       // Integrator Implementation Type


    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   errcode, gerrcode;
    GINT   nlat=10, nlong=20;
    GINT   np=1;    // elem 'order'
    GFTYPE radiusi=1, radiuso=2;
    GTVector<GINT> ne(3); // # elements in each direction in 3d
    GString spde; // pde equation being solved
    GString sgrid;// name of JSON grid object to use
    GC_COMM comm = GC_COMM_WORLD;

    typename MyTypes::State u;
    typename MyTypes::State uerr;
    typename MyTypes::Time  t  = 0;
    typename MyTypes::Time  dt = 0.1;

    // Read Burgers prop tree; may ovewrite with
    // certain command line args:
    PropertyTree ptree;       // main param file structure
    PropertyTree eqptree;     // equation props
    PropertyTree gridptree;   // grid props
    PropertyTree dissptree;   // dissipation props
    ptree.load("gburgers.jsn");

    np = ptree.getValue<GINT>("exp_prder");
    eqptree     = ptree.getPropertyTree("adv_equation_traits");
    sgrid       = ptree.getValue<GString>("grid_type");
    gridptree   = ptree.getPropertyTree(sgrid);
    dissptree   = ptree.getPropertyTree("dissipation_traits");
    
#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:j:k;l:o:p:h")) != -1) {
      switch (iopt) {
      case 'i': // get # elements in r/x
          ne[0] = atoi(optarg);
          break;
      case 'j': // get # elements in lat/y
          ne[1] = atoi(optarg);
          break;
      case 'k': // get # elements in long/z
          ne[2] = atoi(optarg);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
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
    
    GPTLstart("gen_grid");

    // Create grid:
    GGrid *grid = GGridFactory(gridptree, gbasis, comm);

    GPTLstop("gen_grid");

    // Create observer(s), equations, integrator:
    std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(wave_speed, grid));
    std::shared_ptr<EqnBase> eqn_base = eqn_impl;

    // Create the Integrator Implementation
    IntImpl::Traits traits;
    traits.cfl_min = 0;
    traits.cfl_max = 9999;
    traits.dt_min  = 0;
    traits.dt_max  = 9999;

    traits.doheat  = ptree.

/*
    std::shared_ptr<ObsImpl> obs_impl(new ObsImpl());
    std::shared_ptr<ObsBase> obs_base = obs_impl;
*/

    std::shared_ptr<IntImpl> int_impl(new IntImpl(eqn_base,obs_base,traits));



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
