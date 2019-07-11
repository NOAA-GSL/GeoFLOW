//==================================================================================
// Module       : geoglow.cpp
// Date         : 7/7/19 (DLR)
// Description  : GeoFLOW main driver
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include <typeinfo>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "gburgers.hpp"
#include "ggrid_box.hpp"
#include "ggrid_factory.hpp"
#include "gposixio_observer.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/integrator_factory.hpp"
#include "pdeint/stirrer_factory.hpp"
#include "pdeint/observer_factory.hpp"
#include "pdeint/null_observer.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "gtools.h"

using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;


template< // default template arg types
typename StateType = GTVector<GTVector<GFTYPE>*>,
typename GridType  = GGrid,
typename ValueType = GFTYPE,
typename DerivType = StateType,
typename TimeType  = ValueType,
typename JacoType  = StateType,
typename SizeType  = GSIZET
>
struct EquationTypes {
        using State      = StateType;
        using Grid       = GridType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};

using MyTypes     = EquationTypes<>;           // Define types used
using EqnBase     = EquationBase<MyTypes>;     // Equation Base Type
using EqnBasePtr  = std::shared_ptr<EqnBase>;  // Equation Base Ptr
using EqnImpl     = GBurgers<MyTypes>;         // Equation Implementa
using StirBase    = StirrerBase<MyTypes>;      // Stirring Base Type
using StirBasePtr = std::shared_ptr<StirBase>; // Stirring Base Ptr
using ObsBase     = ObserverBase<EqnBase>;     // Observer Base Type
using BasisBase   = GTVector<GNBasis<GCTYPE,GFTYPE>*>; // Basis pool type


GBOOL  bench_=FALSE;
GINT   nsolve_;  // number *solved* 'state' arrays
GINT   nstate_;  // number 'state' arrays
GINT   ntmp_  ;  // number tmp arrays
GGrid *grid_;
State u_;
State c_;
State ub_;
State uf_;
State utmp_;
GTVector<GFTYPE> nu_(3);
BasisBase gbasis_(GDIM);
PropertyTree ptree;       // main prop tree
GGFX<GFTYPE> ggfx_;       // DSS operator
GC_COMM      comm_ ;      // communicator


void allocate(const PropertyTree &ptree);
void deallocate();
void update_dirichlet(const Time &t, State &u, State &ub);
void update_forcing(const Time &t, State &u, State &uf);
void steptop_callback(const Time &t, State &u, const Time &dt);
void create_observers(PropertyTree &ptree, GSIZET icycle, Time time, 
std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers);
void create_equation(PropertyTree &ptree, EqnBasePtr &pEqn);
void create_stirrer(PropertyTree &ptree, StirBasePtr &pStirrer);
void gresetart(PropertyTree &ptree);
void do_bench(GString sbench, GSIZET ncyc);

//#include "init_pde.h"


int main(int argc, char **argv)
{
    GString serr ="geoflow: ";
    GINT    iopt;
    GSIZET  itindex=0;      // restart flag/index
    GSIZET  icycle=0;       // curr time cycle
    std::vector<GINT> pstd(GDIM);  // order in each direction
    GTMatrix<GINT> p; // needed for restart, but is dummy

    typename MyTypes::Time  t  = 0;
    typename MyTypes::Time  dt = 0.1;

    // Initialize comm & global environment:
    mpixx::environment  env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc,argv); 
    GlobalManager::startup();
    comm_  = world; // need this for solver(s) & grid

    // Read main prop tree; may ovewrite with
    // certain command line args:
    EH_MESSAGE("geoflow::call load prop tree...");
    ptree    = InputManager::getInputPropertyTree();       // main param file structure
    EH_MESSAGE("geoflow: prop tree loaded.");

    // Create other prop trees for various objects:
    itindex     = ptree.getValue <GSIZET>("restart_index");
    pstd        = ptree.getArray   <GINT>("exp_order");
    bench_      = ptree.getValue  <GBOOL>("benchmark");

    // Parse command line. ':' after char
    // option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "i:bh")) != -1) {
      switch (iopt) {
      case 'i': // handled by InputManager
          break;
      case 'b': // do benchmark
          bench_ = TRUE;
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-b]" << std::endl;
          std::cout << "Note: '-b' sets benchmark flag" << std::endl;
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


#if defined(_G_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);
#endif
    // Initialize timer:
    GTimerInit();

    //***************************************************
    // Create basis pool:
    //***************************************************
    EH_MESSAGE("geoflow: create basis pool...");
    create_basis_pool(ptree, gbasis_);

    //***************************************************
    // Create grid:
    //***************************************************
    EH_MESSAGE("geoflow: build grid...");
    GTimerStart("gen_grid");

    grid_ = GGridFactory::build(ptree, gbasis_, comm_);
    GTimerStop("gen_grid");

    //***************************************************
    // Initialize gather/scatter operator:
    //***************************************************
    EH_MESSAGE("geoflow: initialize gather/scatter...");
    GTimerStart("init_ggfx_op");

    init_ggfx(ptree, *grid_, ggfx_);

    GTimerStop("init_ggfx_op");
    EH_MESSAGE("geoflow: gather/scatter initialized.");

    //***************************************************
    // Create state and tmp space:
    //***************************************************
    EH_MESSAGE("geoflow: allocate tmp space...");
    allocate(ptree);

    //***************************************************
    // Create equation set:
    //***************************************************
    EH_MESSAGE("geoflow: create equation...");
    EqnBasePtr pEqn;
    create_equation(ptree, pEqn);

    //***************************************************
    // Create the stirrer (to update forcing)
    //***************************************************
    EH_MESSAGE("geoflow: create stirrer...");
    StirBasePtr pStirrer;
    create_stirrer(ptree, pStirrer);

    //***************************************************
    // Initialize state:
    //***************************************************
    EH_MESSAGE("geoflow: Initializing state...");
    if ( itindex == 0 ) { // start new run
      icycle = 0; t = 0.0; 
      compute_analytic(*grid_, t, ptree, u_);
    }
    else {                // restart run
      gio_restart(ptree, 0, u_, p, icycle, t, comm_);
    }

    //***************************************************
    // Create observers: 
    //***************************************************
    EH_MESSAGE("geoflow: create observers...");
    std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> pObservers(new std::vector<std::shared_ptr<ObserverBase<MyTypes>>>());
    create_observers(ptree, icycle, t, pObservers);

    //***************************************************
    // Create integrator:
    //***************************************************
    EH_MESSAGE("geoflow: create integrator...");
    auto pIntegrator = IntegratorFactory<MyTypes>::build(ptree, eqn_base, pStirrer, pObservers, *grid_);
    pIntegrator->get_traits().cycle = icycle;


    //***************************************************
    // Do time integration (output included
    // via observer(s)):
    //***************************************************
    EH_MESSAGE("geoflow: do time stepping...");
    GTimerReset();
    GTimerStart("time_loop");

    pIntegrator->time_integrate(t, uf_, ub_, u_);

    GTimerStop("time_loop");
    EH_MESSAGE("geoflow: time stepping done.");

    //***************************************************
    // Do benchmarking if required:
    //***************************************************
    do_bench("benchmark.txt", pIntegrator->get_numsteps());
 
#if defined(_G_USE_GPTL)
//  GPTLpr(myrank);
    GPTLpr_file("timings.txt");
    GPTLpr_summary();
#endif
    GTimerFinal();

    //***************************************************
    // Do shutdown, cleaning:
    //***************************************************
    EH_MESSAGE("geoflow: do shutdown...");
    GlobalManager::shutdown();
    GlobalManager::finalize();
    deallocate();

    return(0);

} // end, geoflow


//**********************************************************************************
//**********************************************************************************
// METHOD : update_dirichlet
// DESC   : update/set Dirichlet vectors, ub
// ARGS   : t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void update_dirichlet(const Time &t, State &u, State &ub)
{

  Time  tt = t;
  
/*
  if ( (*igbdy)[GBDY_DIRICHLET].size() > 0 ) {
    compute_analytic(*grid_, tt, ptree,  ua_);
  }

  // ...GBDY_DIRICHLET:
  // Set from State vector, ua_:
  for ( auto k=0; k<u.size(); k++ ) { 
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size(); j++ ) {
      (*ub[k])[j] = (*ua_[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
  }
*/

} // end of method update_dirichlet


//**********************************************************************************
//**********************************************************************************
// METHOD: steptop_callback
// DESC  : Top-of-time-step callback ('backdoor') function. 
//         This function might, e.g. update the linear advection 
//         vel. components as a function of time, since they
//         are not solved for in a PDE, but are, rather, prescribed.
// ARGS  : 
//         t  : time
//         u  : current state
//         dt : time step
//**********************************************************************************
void steptop_callback(const Time &t, State &u, const Time &dt)
{
  
#if 0
  Time tt = t;
  compute_analytic(*grid_, tt, ptree, ua_);
#endif

} // end, method steptop_callback


//**********************************************************************************
//**********************************************************************************
// METHOD: create_equation
// DESC  : Create equation implementation
// ARGS  : ptree   : Main property tree
//         pEqn    : EqnBasePtr pointer that is configured and returned
//**********************************************************************************
void create_equation(PropertyTree &ptree, EqnBasePtr &pEqn)
{
  pEqn = EquationFactory<MyTypes>::build(ptree, *grid_, utmp_);

  // Set PDE callback functions, misc:
  std::function<void(const Time &t, State &u, 
                                      State &ub)>  
      fcallback = update_dirichlet; // set tmp function with proper signature for...
  std::function<void(const Time &t, State &u, const Time &dt)> 
      stcallback = steptop_callback; // set tmp function with proper signature for...
  pEqn->set_bdy_update_callback(fcallback); // bdy update callback
  pEqn->set_steptop_callback(stcallback);   // 'back-door' callback
  pEqn->set_nu(nu_);                        // dissipation (may be altered in initial conditions)

} // end method create_equation


//**********************************************************************************
//**********************************************************************************
// METHOD: create_stirrer
// DESC  : Create forcing functions from main ptree
// ARGS  : ptree   : Main property tree
//         pStirrer: StirBasePtr pointer that is configured and returned
//**********************************************************************************
void create_stirrer(PropertyTree &ptree, StirBasePtr &pStirrer)
{

  pStirrer = StirrerFactory<MyTypes>::build(ptree, *grid_);

  // Set stirrer update callback functions:
  std::function<void(const Time &t, State &u, State &uf)>  
      fcallback = update_forcing; // set tmp function with proper signature for...
  pStirrer->set_update_callback(fcallback); // forcing update callback

} // end method create_stirrer


//**********************************************************************************
//**********************************************************************************
// METHOD: create_observers
// DESC  : Create observer list from main ptree
// ARGS  : grid      : GGrid object
//         icycle    : initial icycle
//         time      : initial time
//         pObservers: gather/scatter op, GGFX
//**********************************************************************************
void create_observers(PropertyTree &ptree, GSIZET icycle, Time time,
std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers)
{
    GINT    ivers;
    GSIZET  rest_ocycle;       // restart output cycle
    GSIZET  deltac;            // cycle interval
    GFTYPE  ofact;             // output freq in terms of restart output
    Time    deltat;            // time interval
    PropertyTree obsptree;     // observer props 
    GString dstr = "none";
    GString ptype;
    GString ctype;

    if ( bench_ ) return;

    std::vector<GString> default_obslist; default_obslist.push_back(dstr);
    std::vector<GString> obslist = ptree.getArray<GString>("observer_list",default_obslist);
    dstr = "constant";
    ptype = ptree.getValue<GString>("exp_order_type",dstr);

    // Tie cadence_type to restart type:
    obsptree = ptree.getPropertyTree("posixio_observer");
    ctype    = obsptree.getValue<GString>("cadence_type");
   
    // If doing a restart, set observer output
    // cycles to value relative to the restart output cycle:
    rest_ocycle = ptree.getValue <GSIZET>("restart_index");
    if ( "constant" == ptype ) ivers = 0;
    if ( "variable" == ptype ) ivers = 1;
    for ( GSIZET j=0; j<obslist.size(); j++ ) {
      if ( "none" != obslist[j] ) {
        obsptree = ptree.getPropertyTree(obslist[j]);
        // Set output version based on exp_order_type:
        if ( "constant" == ptype 
         && "posixio_observer" == obslist[j]  ) obsptree.setValue<GINT>("misc",ivers);

        ofact       = obsptree.getValue<GDOUBLE>("interval_freq_fact",1.0);
        deltat      = obsptree.getValue<GDOUBLE>("time_interval",0.01);
        deltac      = obsptree.getValue<GDOUBLE>("cycle_interval",1);
        // Set current time and output cycle so that observer can initialize itself
        // These should be hidden from the config file:
        if ( "posixio_observer" == obslist[j]  ) ofact = 1.0;
        obsptree.setValue <GSIZET>("start_ocycle",MAX(0.0,rest_ocycle*ofact));
        obsptree.setValue <GFTYPE>("start_time"  ,time);
        obsptree.setValue <GFTYPE>("time_interval", MAX(0.0,deltat/ofact));
        obsptree.setValue <GSIZET>("cycle_interval",MAX(1.0,deltac/ofact));
        obsptree.setValue<GString>("cadence_type",ctype);

        pObservers->push_back(ObserverFactory<MyTypes>::build(obsptree,*grid_));
      }
    }

    for ( GSIZET j=0; j<pObservers->size(); j++ ) (*pObservers)[j]->set_tmp(utmp_);

} // end method create_observers


//**********************************************************************************
//**********************************************************************************
// METHOD: create_basis_pool
// DESC  : Create basis pool from prop tree
// ARGS  : ptree     : main property tree
//         gbasis    : array of allowed basis objects
//**********************************************************************************
void create_basis_pool(PropertyTree &ptree, BasisBase &gbasis);
{
    
  // Eventually, this may become an actual pool, from which
  // solvers will determine basis in each direction. For now...
  for ( auto j=0; j<GDIM; j++ ) {
    gbasis [j] = new GLLBasis<GCTYPE,GFTYPE>(pstd[j]);
  }

} // end method create_basis_pool


//**********************************************************************************
//**********************************************************************************
// METHOD: do_bench
// DESC  : Do benchmark from timers
// ARGS  : fname     : filename
//         ncyc      : number time cycles to average over
//**********************************************************************************
void do_bench(GString fname, GSIZET ncyc)
{
    if ( !bench_ ) return;

#if defined(_G_USE_GPTL)

    GINT   myrank   = GComm::WorldRank(comm_);
    GINT   ntasks   = GComm::WorldSize(comm_);
    GINT   nthreads = 0;
    GFTYPE dxmin, lmin;
    GFTYPE ttotal;
    GFTYPE tggfx;
    GFTYPE texch;
    std::ifstream itst;
    std::ofstream ios;
    GTVector<GSIZET> lsz(2), gsz(2);

    // Get global no elements and dof:
    lsz[0] = grid_->nelems();
    lsz[1] = grid_->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);
    if ( myrank == 0 ) {
      itst.open(fname);
      ios.open(fname,std::ios_base::app);
  
      // Write header, if required:
      if ( itst.peek() == std::ofstream::traits_type::eof() ) {
        ios << "#nelems"  << "  ";
        ios << "ndof"     << "  ";
        ios << "dxmin"    << "  ";
        ios << "lmin"     << "  ";
        ios << "ntasks"   << "  ";
        ios << "nthreads" << "  ";
        ios << "ttotal"   << "  ";
        ios << "tggfx"    << "  ";
        ios << "texch"           ;
        ios << endl;
      }
      itst.close();

      GPTLget_wallclock("time_loop"     , 0,  &ttotal); ttotal /= ncyc;
      GPTLget_wallclock("ggfx_doop"     , 0,  &tggfx ); tggfx  /= ncyc;
      GPTLget_wallclock("ggfx_doop_exch", 0,  &texch ); texch  /= ncyc;

      dxmin = grid_->minnodedist();
      lmin  = grid_->minlength();
  
      ios << gsz[0]          << "   " ;
      ios << gsz[1]          << "   " ;
      ios << dxmin           << "   " ;
      ios << lmin            << "   " ;
      ios << ntasks          << "   " ;
      ios << nthreads        << "   ";
      ios << ttotal          << "   ";
      ios << tggfx           << "   ";
      ios << texch                   ;
      ios << endl;

      ios.close();
    }
#endif

    return;

} // end method do_bench


//**********************************************************************************
//**********************************************************************************
// METHOD: allocate
// DESC  : Allocate state, tmp arrays
// ARGS  : ptree:  main prop tree
//**********************************************************************************
void allocate(const PropertyTree &ptree)
{

  GBOOL        doheat, bpureadv, bforced;
  GINT         nadv, nforced;
  std::vector<GINT>
               ibounded, iforced, diforced;
  std::string  sgrid;
  std::string  eq_name   = ptree.getValue("pde_name");
  PropertyTree eqn_ptree = ptree.getPropertyTree(equation_name);
  PropertyTree stp_ptree = ptree.getPropertyTree("stepper_props");

  assert("3##$%!62ahTze32934Plq1C4" != eq_name
      && "pde_name required");

  if ( "pde_burgers" == eq_name ) {
    sgrid     = ptree.getValue<GString>  ("grid_type");
    doheat    = eqn_ptree.getValue<bool> ("doheat",false);
    bpureadv  = eqn_ptree.getValue<bool> ("bpureadv",false);
    bforced   = eqn_ptree.getValue<bool> ("use_forcing",false);
    for ( auto i=0; i<GDIM; i++ ) diforced.push_back(i);
    iforced   = eqn_ptree.getArray<GINT> ("forcing_comp", diforced);
    nadv      = sgrid == "grid_icos" ? 3 : GDIM;
    nsolve_   = GDIM;
    nstate_   = GDIM;
    if ( doheat || bpureadv ) {
      nsolve_   = 1;
      nstate_   = nadv + nsolve_;
    }
    if ( "grid_icos" != sgrid ) {
      ibounded.resize(nsolve_);
      for ( auto i=0; i<nsolve_; i++ ) ibounded.push_back(i);
    }
    ntmp_     = 24;
  }
  
  nforced = MIN(nsolved_,iforced.size());

  u_   .resize(nstate_);                // state
  ub_  .resize(nstate_); ub_ = NULLPTR  // bdy state array
  uf_  .resize(nstate_); uf_ = NULLPTR; // forcing array
  utmp_.resize(ntmp_);                  // tmp array
  c_   .resize(nadv);                   // adv. velocity

  for ( auto j=0; j<u_      .size(); j++ ) u_             [j] = new GTVector<GFTYPE>(grid_->size());

  for ( auto j=0; j<ibounded.size(); j++ ) ub_  [ibounded[j]] = new GTVector<GFTYPE>(grid_->nbdydof());

  if ( bforced ) {
    for ( auto j=0; j<nforced      ; j++ ) uf_  [iforced[j]] = new GTVector<GFTYPE>(grid_->ndof());
  }

  for ( auto j=0; j<utmp_    .size(); j++ ) utmp_        [j] = new GTVector<GFTYPE>(grid_->size());

  if ( doheat || bpureadv ) { // assign linear adv velocity:
    for ( auto j=0; j<c_.size(); j++ ) c_   [j] = u_[j+1];
  }

} // end method allocate


//**********************************************************************************
//**********************************************************************************
// METHOD: deallocate
// DESC  : De-allocate state, tmp arrays
//**********************************************************************************
void deallocate()
{

  if ( grid_ != NULLPTR )                   delete grid_;
  for ( auto j=0; j<gbasis_.size(); j++ ) delete gbasis_[j];
  for ( auto j=0; j<utmp_  .size(); j++ ) delete utmp_  [j];
  for ( auto j=0; j<u_     .size(); j++ ) delete u_     [j];
  for ( auto j=0; j<ub_    .size(); j++ ) delete ub_    [j];
  for ( auto j=0; j<uf_    .size(); j++ ) delete uf_    [j];

} // end method deallocate
