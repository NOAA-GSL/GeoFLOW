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
#include "ggrid_factory.hpp"
#include "ggfx.hpp"
#include "gmorton_keygen.hpp"
#include "gburgers.hpp"
#include "pdeint/equation_base.hpp"

using namespace geoflow::pdeint;


template<
typename StateType = GTVector<GTVectorGFTYPE>*>,
typename ValueType = GFTYPE,
typename DerivType = StateType,
typename TimeType  = ValueType,
typename JacoType  = NULLPTR,
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


void compute_analytic(GGrid &grid, Time &t, const tbox::PropertyTree& ptree,  State &ua);

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
    GINT   np=1;    // elem 'order'
    GSISET maxSteps;
    GFTYPE radiusi=1, radiuso=2;
    GTVector<GINT> ne(3); // # elements in each direction in 3d
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
    PropertyTree stepptree;   // stepper props
    PropertyTree dissptree;   // dissipation props
    PropertyTree tintptree;   // time integration props
    ptree.load("gburgers.jsn");

    // Create other prop trees for various objects:
    np = ptree.getValue<GINT>("exp_prder");
    eqptree     = ptree.getPropertyTree("adv_equation_traits");
    sgrid       = ptree.getValue<GString>("grid_type");
    gridptree   = ptree.getPropertyTree(sgrid);
    stepptree   = ptree.getPropertyTree("stepper_props");
    dissptree   = ptree.getPropertyTree("dissipation_traits");
    tintptree   = ptree.getPropertyTree("time_integration");
    
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

    // Overwrite prop tree traits based on command line args:
    std::vector<GINT> stdne(ne.size());
    for ( auto j=0; j<ne.size; j++ ) stdne.push_back(ne[j]);
    gridptree.setArray("num_elems",stdne);
    ptree.setValue("exp_order",np);

    // Set solver traits from prop tree:
    GBurgers::Traits solver_traits;
    solver_traits.doheat     = eqptree  .getValue("doheat");
    solver_traits.bpureadv   = eqptree  .getValue("bpureadv");
    solver_traits.bconserved = eqptree  .getValue("bconserved");
    solver_traits.itorder    = stepptree.getValue("time_deriv_order");
    solver_traits.inorder    = stepptree.getValue("extrap_order");
    GTVector<GString> ssteppers;
    for ( GSIZET j=0; j<GSTEPPER_MAX; j++ ) ssteppeers.push_back(sGStepperType[j]);
    GSIZET itype; 
    assert(ssteppers.contains(stepptree.getValue("stepping_method"),itype) 
           && "Invalide stepping method in JSON file");
    }
    solver_traits.steptype   = static_case<GStepperType>(itype);
    

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
    std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(*grid, u, solver_traits, utmp));
    std::shared_ptr<EqnBase> eqn_base = eqn_impl;


    // Create the Integrator Implementation
    IntImpl::Traits int_traits;
    int_traits.cfl_min = 0;
    int_traits.cfl_max = 9999;
    int_traits.dt_min  = 0;
    int_traits.dt_max  = 9999;

    std::shared_ptr<ObsImpl> obs_impl(new ObsImpl());
    std::shared_ptr<ObsBase> obs_base = obs_impl;

    std::shared_ptr<IntImpl> int_impl(new IntImpl(eqn_base,obs_base,int_traits));
    dt       = tintptree.getVaue("dt"); 
    maxSteps = tintptree.getValue("cycle_end");

    GPTLstart("time_loop");
    for( GSIZET i=0; i<maxSteps; i++ ){
      eqn_base->step(t,dt,u);
      t += dt;
    }
    GPTLststop("time_loop");

#if 0
    GTVector<GFTYPE> lerrnorm(3), gerrnorm(3);
    compute_analytic(*grid, t, eqptree, ua);
    for ( GSIZET j=0; j<u.size(); i++ ) {
      *utmp[0] = *u[j] - *ua[j];
       lerrnorm[0] = utmp[0]->L1norm (); // inf-norm
       lerrnorm[1] = utmp[0]->Eucnorm(); // Euclidean-norm
       lerrnorm[1] *= lerrnorm[1]; 
       // Add spatial intregration for L2 computation
    }
    // Accumulate errors:
    GComm::Allreduce(lerrnorm.data()  , gerrnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
    GComm::Allreduce(lerrnorm.data()+1, gerrnorm.data()+1, 2, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm)
   
    GSIZET lnelems=grid->nelems();
    GSIZET gnelems;
    GComm::Allreduce(&lnelems, &gnelems, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm)

    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open(burgers_err.txt");
    ios.open("burgers_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p  num_elems   L1_err   Eucl_err  L2_err" << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  "
        << "  " << gerrnorm[0] << "  " << gerrnorm[1] << "  " << gerrnorm[2]
        << std::endl;
    ios.close();

#endif
 
    GPTLpr_file("timing.txt");
    GPTLfinalize();


    GComm::TermComm();
    if ( grid != NULLPTR ) delete grid;

    return(0);

} // end, main
