//==================================================================================
// Module       : geoglow.cpp
// Date         : 7/7/19 (DLR)
// Description  : GeoFLOW main driver
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "geoflow.h"

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
    ptree_    = InputManager::getInputPropertyTree();       // main param file structure
    EH_MESSAGE("geoflow: prop tree loaded.");

    // Create other prop trees for various objects:
    itindex     = ptree_.getValue <GSIZET>("restart_index");
    pstd        = ptree_.getArray   <GINT>("exp_order");
    bench_      = ptree_.getValue  <GBOOL>("benchmark");

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
    create_basis_pool(ptree_, gbasis_);

    //***************************************************
    // Create grid:
    //***************************************************
    EH_MESSAGE("geoflow: build grid...");
    GTimerStart("gen_grid");

    grid_ = GGridFactory::build(ptree_, gbasis_, comm_);
    GTimerStop("gen_grid");

    //***************************************************
    // Initialize gather/scatter operator:
    //***************************************************
    EH_MESSAGE("geoflow: initialize gather/scatter...");
    GTimerStart("init_ggfx_op");

    init_ggfx(ptree_, *grid_, ggfx_);
    grid_->set_ggfx(*ggfx_);

    GTimerStop("init_ggfx_op");
    EH_MESSAGE("geoflow: gather/scatter initialized.");

    //***************************************************
    // Create state and tmp space:
    //***************************************************
    EH_MESSAGE("geoflow: allocate tmp space...");
    allocate(ptree_);

    //***************************************************
    // Create equation set:
    //***************************************************
    EH_MESSAGE("geoflow: create equation...");
    create_equation(ptree_, pEqn_);

    //***************************************************
    // Create the mixer (to update forcing)
    //***************************************************
    EH_MESSAGE("geoflow: create mixer...");
    MixBasePtr pMixer;
    create_mixer(ptree_, pMixer);

    //***************************************************
    // Create observers: 
    //***************************************************
    EH_MESSAGE("geoflow: create observers...");
    std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> pObservers(new std::vector<std::shared_ptr<ObserverBase<MyTypes>>>());
    create_observers(pEqn_, ptree_, icycle, t, pObservers);

    //***************************************************
    // Create integrator:
    //***************************************************
    EH_MESSAGE("geoflow: create integrator...");
    auto pIntegrator = IntegratorFactory<MyTypes>::build(ptree_, pEqn_, pMixer, pObservers, *grid_);
    pIntegrator->get_traits().cycle = icycle;

    //***************************************************
    // Initialize state:
    //***************************************************
    EH_MESSAGE("geoflow: Initializing state...");
    if ( itindex == 0 ) { // start new run
      icycle = 0; t = 0.0; 
      init_state(ptree_, *grid_, pEqn_, t, utmp_, u_, ub_);
      init_bdy  (ptree_, *grid_, pEqn_, t, utmp_, u_, ub_);
      init_force(ptree_, *grid_, pEqn_, t, utmp_, u_, uf_);
    }
    else {                // restart run
      gio_restart(ptree_, 0, u_, p, icycle, t, comm_);
    }

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
  

} // end, method steptop_callback


//**********************************************************************************
//**********************************************************************************
// METHOD: create_equation
// DESC  : Create equation implementation
// ARGS  : ptree   : Main property tree
//         pEqn    : EqnBasePtr pointer that is configured and returned
//**********************************************************************************
void create_equation(const PropertyTree &ptree, EqnBasePtr &pEqn)
{
  pEqn = EquationFactory<MyTypes>::build(ptree, *grid_, utmp_);

  // Set PDE callback functions, misc:
  std::function<void(const Time &t, State &u, State &ub)>  
      fcallback = [](const Time &t, State &u, State &ub)
                  {update_boundary(t, u, ub);}; // set tmp function with proper signature for...
  pEqn->set_bdy_update_callback(fcallback); // bdy update callback

#if 0
  std::function<void(const Time &t, State &u, const Time &dt)> 
      stcallback = [](const Time &t, State &u, const Time &dt)
                   {steptop_callback(t, u, dt);}; // set tmp function with proper signature for...
  pEqn->set_steptop_callback(stcallback);   // 'back-door' callback
#endif

} // end method create_equation


//**********************************************************************************
//**********************************************************************************
// METHOD: create_mixer
// DESC  : Create forcing functions from main ptree
// ARGS  : ptree   : Main property tree
//         pMixer: MixBasePtr pointer that is configured and returned
//**********************************************************************************
void create_mixer(PropertyTree &ptree, MixBasePtr &pMixer)
{

  pMixer = MixerFactory<MyTypes>::build(ptree, *grid_);

#if 0
  // Set mixer update callback functions:
  std::function<void(const Time &t, State &u, State &uf)>  
      fcallback = update_forcing; // set tmp function with proper signature for...
  pMixer->set_update_callback(fcallback); // forcing update callback
#endif

} // end method create_mixer


//**********************************************************************************
//**********************************************************************************
// METHOD: create_observers
// DESC  : Create observer list from main ptree
// ARGS  : grid      : GGrid object
//         icycle    : initial icycle
//         time      : initial time
//         pObservers: gather/scatter op, GGFX
//**********************************************************************************
void create_observers(EqnBasePtr &pEqn, PropertyTree &ptree, GSIZET icycle, Time time,
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

        pObservers->push_back(ObserverFactory<MyTypes>::build(obsptree, pEqn, *grid_));
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
void create_basis_pool(PropertyTree &ptree, BasisBase &gbasis)
{
  std::vector<GINT> pstd(GDIM);  // order in each direction

  pstd = ptree.getArray<GINT>("exp_order");
    
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
  std::string  eqn_name  = ptree.getValue<GString>("pde_name");
  PropertyTree eqn_ptree = ptree.getPropertyTree  (eqn_name);
  PropertyTree stp_ptree = ptree.getPropertyTree  ("stepper_props");
  bforced                = ptree.getValue<GBOOL>  ("use_forcing",false);

  assert("3##$%!62ahTze32934Plq1C4" != eqn_name
      && "pde_name required");

  if ( "pde_burgers" == eqn_name ) {
    sgrid     = ptree.getValue<GString>  ("grid_type");
    doheat    = eqn_ptree.getValue<bool> ("doheat",false);
    bpureadv  = eqn_ptree.getValue<bool> ("bpureadv",false);
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
  
  nforced = MIN(nsolve_,iforced.size());

  u_   .resize(nstate_);                // state
  ub_  .resize(nstate_); ub_ = NULLPTR; // bdy state array
  uf_  .resize(nstate_); uf_ = NULLPTR; // forcing array
  utmp_.resize  (ntmp_);                // tmp array

  for ( auto j=0; j<u_. size(); j++ ) u_                [j] = new GTVector<GFTYPE>(grid_->size());

  for ( auto j=0; j<ibounded.size(); j++ ) ub_[ibounded[j]] = new GTVector<GFTYPE>(grid_->nbdydof());

  if ( bforced ) {
    for ( auto j=0; j<nforced      ; j++ ) uf_ [iforced[j]] = new GTVector<GFTYPE>(grid_->ndof());
  }

  for ( auto j=0; j<utmp_   .size(); j++ ) utmp_        [j] = new GTVector<GFTYPE>(grid_->size());


} // end method allocate


//**********************************************************************************
//**********************************************************************************
// METHOD: deallocate
// DESC  : De-allocate state, tmp arrays
// ARGS  : none.
//**********************************************************************************
void deallocate()
{

  if ( grid_ != NULLPTR )                 delete grid_;
  for ( auto j=0; j<gbasis_.size(); j++ ) delete gbasis_[j];
  for ( auto j=0; j<utmp_  .size(); j++ ) delete utmp_  [j];
  for ( auto j=0; j<u_     .size(); j++ ) delete u_     [j];
  for ( auto j=0; j<ub_    .size(); j++ ) delete ub_    [j];
  for ( auto j=0; j<uf_    .size(); j++ ) delete uf_    [j];

} // end method deallocate


//**********************************************************************************
//**********************************************************************************
// METHOD: init_state
// DESC  : top-level method to set initial conditions.
// ARGS  : ptree: main prop tree
//         grid : grid object
//         peqn : pointer to EqnBase 
//         t    : initial time
//         utmp : vector of tmp vectors
//         u    : full state vector
//         ub   : full boundary state vector
//**********************************************************************************
void init_state(const PropertyTree &ptree, GGrid &grid, EqnBasePtr &peqn, Time &t, State &utmp, State &u, State &ub)
{
  GBOOL bret;

  bret = GInitStateFactory<MyTypes>::init(ptree, grid, peqn, t, utmp, ub, u);

  assert(bret && "state initialization failed");

} // end method init_state


//**********************************************************************************
//**********************************************************************************
// METHOD: init_force
// DESC  : top-level method to set initial forcing.
// ARGS  : ptree: main prop tree
//         grid : grid object
//         peqn : pointer to EqnBase 
//         t   : initial time
//         utmp: vector of tmp vectors 
//         u   : full state vector
//         uf  : full boundary state vector
//**********************************************************************************
void init_force(const PropertyTree &ptree, GGrid &grid, EqnBasePtr &peqn, Time &t, State &utmp, State &u, State &uf)
{
  GBOOL bret;

  bret = GInitForceFactory<MyTypes>::init(ptree, grid, peqn, t, utmp, u, uf);

  assert(bret && "forcing initialization failed");
  
} // end method init_force


//**********************************************************************************
//**********************************************************************************
// METHOD: init_bdy
// DESC  : top-level method to set initial bdy conditions.
// ARGS  : ptree: main prop tree
//         grid : grid object
//         peqn : pointer to EqnBase 
//         t    : initial time
//         utmp : vector of tmp vectors 
//         u    : full state vector
//         ub   : full boundary state vector
//**********************************************************************************
void init_bdy(const PropertyTree &ptree, GGrid &grid, EqnBasePtr &peqn, Time &t, State &utmp, State &u, State &ub)
{
  GBOOL bret;

  bret = GInitBdyFactory<MyTypes>::init(ptree, grid, peqn, t, utmp, u, ub);

  assert(bret && "boundary initialization failed");
  
} // end method init_bdy


//**********************************************************************************
//**********************************************************************************
// METHOD : update_boundary
// DESC   : update/set boundary vectors, ub
// ARGS   : 
//          t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void update_boundary(const Time &t, State &u, State &ub)
{
  GBOOL  bret;
  GFTYPE tt = t;

  bret = GUpdateBdyFactory<MyTypes>::update(ptree_, *grid_, pEqn_, tt, utmp_, u, ub);
  
  assert(bret && "boundary update failed");
  
} // end of method update_boundary


//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : ptree   : main property tree
//         grid    : GGrid object pointer, instantiated here
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(PropertyTree &ptree, GGrid &grid, GGFX<GFTYPE> *&ggfx)
{
  GString                        serr = "init_ggfx: ";
  GFTYPE                         delta;
  GFTYPE                         rad;
  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
  GTPoint<GFTYPE>                dX, porigin, P0;
  GTVector<GNODEID>              glob_indices;
  GTVector<GTVector<GFTYPE>>    *xnodes;
  GString                        sgrid;
  std::vector<GFTYPE>            pstd;
  PropertyTree                   gtree;

  sgrid = ptree.getValue<GString>("grid_type");
  gtree = ptree.getPropertyTree(sgrid);

  if ( sgrid == "grid_box" ) {
    pstd = gtree.getArray<GFTYPE>("xyz0");
    P0   = pstd;
  }
  if ( sgrid == "grid_icos" ) {
    rad   = gtree.getValue<GFTYPE>("radius");
    P0.x1 = -rad ;
    P0.x2 = -rad ;
    P0.x3 = -rad ;
  }
  if ( sgrid == "grid_sphere" ) {
    rad   = gtree.getValue<GFTYPE>("radiuso");
    P0.x1 = -rad ;
    P0.x2 = -rad ;
    P0.x3 = -rad ;
  }

  // First, periodize coords if required to, 
  // before labeling nodes:
  if ( typeid(grid) == typeid(GGridBox) ) { 
    static_cast<GGridBox*>(&grid)->periodize();
  }

  delta  = grid.minnodedist();
  dX     = 0.1*delta;
  xnodes = &grid.xNodes();
  glob_indices.resize(grid.ndof());

  // Integralize *all* internal nodes
  // using Morton indices:
  gmorton.setIntegralLen(P0,dX);

  gmorton.key(glob_indices, *xnodes);

  // Initialize gather/scatter operator:
  GBOOL bret;
  ggfx = new GGFX<GFTYPE>();
  assert(ggfx != NULLPTR && "Cannot instantiate GGFX operator");
  bret = ggfx->init(glob_indices);
  assert(bret && "Initialization of GGFX operator failed");

  // Unperiodize nodes now that connectivity map is
  // generated, so that coordinates mean what they should:
  if ( typeid(grid) == typeid(GGridBox) ) { // periodize coords
    static_cast<GGridBox*>(&grid)->unperiodize();
  }

} // end method init_ggfx

