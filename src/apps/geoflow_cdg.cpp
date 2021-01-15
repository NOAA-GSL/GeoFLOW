//==================================================================================
// Module       : geoflow_cdg.cpp
// Date         : 7/7/19 (DLR)
// Description  : GeoFLOW main driver for CG and DG initial value,
//                boundary value problems
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "geoflow_cdg.h"

int main(int argc, char **argv)
{
    GString serr ="geoflow: ";
    GINT    iopt;
    GSIZET  itindex=0;      // restart flag/index
    GSIZET  icycle=0;       // curr time cycle
    std::vector<GINT> pstd(GDIM);  // order in each direction
    GTMatrix<GINT> p; // needed for restart, but is dummy
    ObsTraitsType binobstraits;
    CommandLine   cline_;

    typename MyTypes::Time  t  = 0;
    typename MyTypes::Time  dt = 0.1;

    // Initialize comm & global environment:
    mpixx::environment  env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc,argv); 
    GlobalManager::startup();
    comm_  = world; // need this for solver(s) & grid
    GEOFLOW_TRACE();

    // Read main prop tree; may ovewrite with
    // certain command line args:
    EH_MESSAGE("geoflow::call load prop tree...");
    ptree_ = InputManager::getInputPropertyTree();
    EH_MESSAGE("geoflow: prop tree loaded.");

    // Create other prop trees for various objects:
    itindex     = ptree_.getValue <GSIZET>("restart_index");
    pstd        = ptree_.getArray   <GINT>("exp_order");
    bench_      = ptree_.getValue  <GBOOL>("benchmark");

    // Process Command Line Flags
    cline_ = InputManager::getInputCommandLine();
    bench_ = bench_ || cline_.exists("b","bench");


    //***************************************************
    // Create basis pool:
    //***************************************************
    EH_MESSAGE("geoflow: create basis pool...");
    create_basis_pool(ptree_, gbasis_);

    //***************************************************
    // Create grid:
    //***************************************************
    EH_MESSAGE("geoflow: build grid...");
    GEOFLOW_TRACE_START("gen_grid");

    ObserverFactory<MyTypes>::get_traits(ptree_, "gio_observer", binobstraits); 
    grid_ = GGridFactory<MyTypes>::build(ptree_, gbasis_, pIO_, binobstraits, comm_);
    GEOFLOW_TRACE_STOP();

    //***************************************************
    // Initialize gather/scatter operator:
    //***************************************************
    EH_MESSAGE("geoflow: initialize gather/scatter...");
    GEOFLOW_TRACE_START("init_ggfx_op");

    init_ggfx(ptree_, *grid_, ggfx_);
    grid_->set_ggfx(*ggfx_);

    GEOFLOW_TRACE_STOP();
    EH_MESSAGE("geoflow: gather/scatter initialized.");

    //***************************************************
    // Set grid terrain:
    //***************************************************
    EH_MESSAGE("geoflow: set grid terrain...");
    GEOFLOW_TRACE_START("do_terrain");

    do_terrain(ptree_, *grid_);

    GEOFLOW_TRACE_STOP();
    EH_MESSAGE("geoflow: terrain added.");

    //***************************************************
    // Create equation set:
    //***************************************************
    EH_MESSAGE("geoflow: create equation...");
    create_equation(ptree_, pEqn_);

    //***************************************************
    // Create state and tmp space:
    //***************************************************
    EH_MESSAGE("geoflow: allocate tmp space...");
    allocate(ptree_);

    //***************************************************
    // Initialize PDE:
    //***************************************************
    EH_MESSAGE("geoflow: initialize PDE...");
    pEqn_->init(u_, utmp_);


    // BTF Early Shutdown
    deallocate();
    EH_MESSAGE("geoflow: do shutdown...");
    GlobalManager::shutdown();
    GlobalManager::finalize();
    GComm::TermComm();
    // BTF Early Shutdown


    //***************************************************
    // Create the mixer (to update forcing)
    //***************************************************
    EH_MESSAGE("geoflow: create mixer...");
    create_mixer(ptree_, pMixer_);

    //***************************************************
    // Create observers: 
    //***************************************************
    EH_MESSAGE("geoflow: create observers...");
    create_observers(pEqn_, ptree_, icycle, t, pObservers_);

    //***************************************************
    // Create integrator:
    //***************************************************
    EH_MESSAGE("geoflow: create integrator...");
    pIntegrator_ = IntegratorFactory<MyTypes>::build(ptree_, pEqn_, pMixer_, pObservers_, *grid_);
    pIntegrator_->get_traits().cycle = icycle;

    //***************************************************
    // Initialize state:
    //***************************************************
    GComm::Synch();
    EH_MESSAGE("geoflow: Initializing state...");
    if ( itindex == 0 ) { // start new run
      icycle = 0; t = 0.0; 
      init_state(ptree_, *grid_, pEqn_, t, utmp_, u_, ub_);
    }
    else {                // restart run
      do_restart(ptree_, *grid_, u_, p, icycle, t);
    }
    init_force(ptree_, *grid_, pEqn_, t, utmp_, u_, uf_);

    //***************************************************
    // Do time integration (output included
    // via observer(s)):
    //***************************************************
    GComm::Synch();
    EH_MESSAGE("geoflow: do time stepping...");
    GEOFLOW_TRACE_START("time_loop");

    pIntegrator_->time_integrate(t, uf_, ub_, u_);

    GEOFLOW_TRACE_STOP();
    EH_MESSAGE("geoflow: time stepping done.");

    //***************************************************
    // Do benchmarking if required:
    //***************************************************
    GEOFLOW_TRACE_START("benchmark_timer");
    do_bench("benchmark.txt", pIntegrator_->get_numsteps());
    GEOFLOW_TRACE_STOP();

    //***************************************************
    // Compare solution if required:
    //***************************************************
    compare(ptree_, *grid_, pEqn_, t, utmp_, ub_, u_);
 
    //***************************************************
    // Do shutdown, cleaning:
    //***************************************************
    deallocate();
    EH_MESSAGE("geoflow: do shutdown...");
    GlobalManager::shutdown();
    GlobalManager::finalize();
    GComm::TermComm();


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
  GEOFLOW_TRACE();
  pEqn = EquationFactory<MyTypes>::build(ptree, *grid_);

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
  GEOFLOW_TRACE();
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
// DESC  : Create IO object and observer list from main ptree. 
// ARGS  : grid      : GGrid object
//         icycle    : initial icycle
//         time      : initial time
//         pObservers: observer list, returned
//**********************************************************************************
void create_observers(EqnBasePtr &pEqn, PropertyTree &ptree, GSIZET icycle, Time time,
std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers)
{
	GEOFLOW_TRACE();
    GINT    ivers;
    GSIZET  rest_ocycle;       // restart output cycle
    GSIZET  deltac, cyc_ref;   // cycle interval
    GFTYPE  ofact ;            // output freq in terms of restart output
    Time    deltat, delt_ref;  // time interval
    GString dstr = "none";
    GString spref;
    GString ctype;
    GString ptype;
    ObserverBase<MyTypes>::Traits
            obstraits;         // observer traits
    PropertyTree 
            obsptree;          // observer props 
    StateInfo 
            stateinfo;         // StateInfo structure
    State   dummy(1);

    if ( bench_ ) return; // don't need IO

    stateinfo = pEqn->stateinfo();

    std::vector<GString> default_obslist; default_obslist.push_back(dstr);
    std::vector<GString> obslist = ptree.getArray<GString>("observer_list",default_obslist);
    dstr = "constant";
    ptype = ptree.getValue<GString>("exp_order_type",dstr);

    // Tie cadence_type to restart type:
    obsptree = ptree.getPropertyTree("gio_observer");
    ctype    = obsptree.getValue<GString>("cadence_type");
   
    // If doing a restart, set observer output
    // cycles to value relative to the restart output cycle:
    rest_ocycle = ptree.getValue <GSIZET>("restart_index");
    if ( "constant" == ptype ) ivers = 0;
    if ( "variable" == ptype ) ivers = 1;

    // Find iobserver, get cadence, and use this to tie 
    // other observer cadences to it:
    for ( GSIZET j=0; j<obslist.size(); j++ ) {
      obsptree = ptree.getPropertyTree(obslist[j]);
      if ( "gio_observer" == obslist[j]  ) {
        delt_ref      = obsptree.getValue<GDOUBLE>("time_interval",0.01);
        cyc_ref       = obsptree.getValue <GSIZET>("cycle_interval",1);
      }
    }

    for ( GSIZET j=0; j<obslist.size(); j++ ) {
      if ( "none" != obslist[j] ) {
        obsptree = ptree.getPropertyTree(obslist[j]);
        // Set output version based on exp_order_type:
        if ( "constant" == ptype 
         && "gio_observer" == obslist[j]  ) obsptree.setValue<GINT>("misc",ivers);

        ofact       = obsptree.getValue<GDOUBLE>("interval_freq_fact",1.0);
        deltat      = obsptree.getValue<GDOUBLE>("time_interval",0.01);
        deltac      = obsptree.getValue <GSIZET>("cycle_interval",1);
        // Set current time and output cycle so that observer can initialize itself
        // These could/should be hidden from the config file:
        if ( "gio_observer" == obslist[j]  ) ofact = 1.0;
        obsptree.setValue <GSIZET>("start_ocycle",MAX(0.0,rest_ocycle*ofact));
        obsptree.setValue <GFTYPE>("start_time"  ,time);

        // Link each observer cadence to I/O cadence:
        if ( "gio_observer" != obslist[j]  ) {
          obsptree.setValue <GFTYPE>("time_interval", MAX(0.0,delt_ref/ofact));
          deltac = (GSIZET)( ((GDOUBLE)cyc_ref)/ofact );
          obsptree.setValue <GSIZET>("cycle_interval",MAX(1  ,deltac));
          obsptree.setValue<GString>("cadence_type",ctype);
        }
        ptree.setPropertyTree(obslist[j],obsptree); // set obs tree with new values

        // Create pIO object factory call:
        if ( pIO_ == NULLPTR ) {
          pIO_ = IOFactory<MyTypes>::build(ptree, *grid_, comm_);
        }
        pObservers->push_back(ObserverFactory<MyTypes>::build(ptree, obslist[j], pEqn, *grid_, pIO_));
        (*pObservers)[j]->set_tmp(utmp_);
        
        if ( "gio_observer" == obslist[j]  ) {
          spref            = obsptree.getValue<std::string>("agg_state_name","state");
          stateinfo.sttype = 0; // grid type filename format
          stateinfo.svars  = obsptree.getArray<std::string>("state_names");
          stateinfo.idir   = obsptree.getValue<std::string>("idir");
          stateinfo.index  = rest_ocycle;
          // If doing a restart, initialize observer with stateinfo data:
          if ( rest_ocycle > 0 ) { 
//          obstraits = (*pObservers)[j]->get_traits();
            pIO_->read_state(spref, stateinfo, dummy, false);
          }
          irestobs_ = j;
        }
        (*pObservers)[j]->init(stateinfo);
      }
    }


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
  GEOFLOW_TRACE();
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
	GEOFLOW_TRACE();
    if ( !bench_ ) return;

#if defined(_G_USE_GPTL)

    GINT   myrank   = GComm::WorldRank(comm_);
    GINT   ntasks   = GComm::WorldSize(comm_);
    GINT   nthreads = 1;
    GFTYPE dxmin, lmin;
    GFTYPE ttotal;
    GFTYPE tggfx;
    GFTYPE texch;
    GFTYPE tgrid;
    GFTYPE tggfxinit;
    std::ifstream itst;
    std::ofstream ios;
    GTVector<GSIZET> lsz(2), gsz(2);

    #pragma omp parallel //num_threads(3)
    {
     nthreads = omp_get_num_threads();
    }

    // Get global no elements and dof & lengths:
    lsz[0] = grid_->nelems();
    lsz[1] = grid_->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);
    dxmin = grid_->minnodedist();
    lmin  = grid_->minlength();
    if ( myrank == 0 ) {
      itst.open(fname);
      ios.open(fname,std::ios_base::app);
  
      // Write header, if required:
      if ( itst.peek() == std::ofstream::traits_type::eof() ) {
        ios << "#nelems"  << "  ";
        ios << "ndof"     << "  ";
        ios << "dxmin"    << "  ";
        ios << "elmin"    << "  ";
        ios << "ntasks"   << "  ";
        ios << "nthreads" << "  ";
        ios << "ttotal"   << "  ";
        ios << "tggfx"    << "  ";
        ios << "texch"    << "  ";
        ios << "tgrid"    << "  ";
        ios << "tggfxinit"       ;
        ios << endl;
      }
      itst.close();

      GPTLget_wallclock("time_loop"     , 0,  &ttotal); ttotal /= ncyc;
      GPTLget_wallclock("ggfx_doop"     , 0,  &tggfx ); tggfx  /= ncyc;
      GPTLget_wallclock("ggfx_doop_exch", 0,  &texch ); texch  /= ncyc;
      GPTLget_wallclock("gen_grid"      , 0,  &tgrid); 
      GPTLget_wallclock("init_ggfx_op"  , 0,  &tggfxinit); 

      ios << gsz[0]          << "   " ;
      ios << gsz[1]          << "   " ;
      ios << dxmin           << "   " ;
      ios << lmin            << "   " ;
      ios << ntasks          << "   " ;
      ios << nthreads        << "   " ;
      ios << ttotal          << "   " ;
      ios << tggfx           << "   " ;
      ios << texch           << "   " ;
      ios << tgrid           << "   " ;
      ios << tggfxinit                ;
      ios << endl                     ;

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
  GEOFLOW_TRACE();
  std::vector<GINT> *iforced;

  nsolve_ =  pEqn_->solve_size();
  nstate_ =  pEqn_->state_size();
  ntmp_   =  pEqn_->tmp_size();
  iforced = &pEqn_->iforced();
  
  u_   .resize(nstate_);                // state
  ub_  .resize(nstate_); ub_ = NULLPTR; // bdy state array
  uf_  .resize(nstate_); uf_ = NULLPTR; // forcing array
  utmp_.resize  (ntmp_);                // tmp array

  for ( auto j=0; j<u_. size(); j++ ) u_                [j] = new GTVector<GFTYPE>(grid_->size());

  for ( auto j=0; j<iforced->size(); j++ ) uf_ [(*iforced)[j]] = new GTVector<GFTYPE>(grid_->ndof());

  for ( auto j=0; j<utmp_   .size(); j++ ) utmp_        [j] = new GTVector<GFTYPE>(grid_->size());

  // If linear adv. prescribed var is set, 
  // point to correct area of u_:
  for ( GINT j=0; j<c_.size(); j++ ) c_[j] = u_[j+1];


} // end method allocate


//**********************************************************************************
//**********************************************************************************
// METHOD: deallocate
// DESC  : De-allocate state, tmp arrays
// ARGS  : none.
//**********************************************************************************
void deallocate()
{
  GEOFLOW_TRACE();
  if ( grid_ != NULLPTR )                 delete grid_;
  if ( ggfx_ != NULLPTR )                 delete ggfx_;
  for ( auto j=0; j<gbasis_.size(); j++ ) delete gbasis_[j];
  for ( auto j=0; j<utmp_  .size(); j++ ) delete utmp_  [j];
  for ( auto j=0; j<u_     .size(); j++ ) delete u_     [j];
  for ( auto j=0; j<ub_    .size(); j++ ) delete ub_    [j];
  for ( auto j=0; j<uf_    .size(); j++ ) delete uf_    [j];

} // end method deallocate


//**********************************************************************************
//**********************************************************************************
// METHOD: init_state
// DESC  : Top-level method to set initial conditions.
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
  GEOFLOW_TRACE();
  GBOOL bret;

  bret = GInitStateFactory<MyTypes>::init(ptree, grid, peqn->stateinfo(), t, utmp, ub, u);

  assert(bret && "state initialization failed");

} // end method init_state


//**********************************************************************************
//**********************************************************************************
// METHOD: init_force
// DESC  : Top-level method to set initial forcing.
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
	GEOFLOW_TRACE();
  GBOOL bret;

  bret = GInitForceFactory<MyTypes>::init(ptree, grid, peqn->stateinfo(), t, utmp, u, uf);

  assert(bret && "forcing initialization failed");
  
} // end method init_force


//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : ptree   : main property tree
//         grid    : GGrid object pointer, instantiated here
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
//void init_ggfx(PropertyTree &ptree, GGrid &grid, GGFX<GFTYPE> *&ggfx)
//{
//	GEOFLOW_TRACE();
//  GString                        serr = "init_ggfx: ";
//  GBOOL                          bret;
//  GINT                           pmax;
//  GFTYPE                         rad;
//  GFTYPE                         tiny = 100.0*std::numeric_limits<GFTYPE>::epsilon();
//  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
//  GTPoint<GFTYPE>                dX, porigin, P0;
//  GTVector<GNODEID>              glob_indices;
//  GTVector<GFTYPE>               delta, ldelta;
//  GTVector<GTVector<GFTYPE>>    *xnodes;
//  State                          cart;
//  State                          xkey;
//  GString                        sgrid;
//  std::vector<GFTYPE>            pstd;
//  PropertyTree                   gtree;
//  State                          tmp;
//
//
//  sgrid = ptree.getValue<GString>("grid_type");
//  gtree = ptree.getPropertyTree(sgrid);
//
//  pmax = 0;
//  for ( auto j=0; j<gbasis_.size(); j++ ) pmax = MAX(pmax, gbasis_[j]->getOrder());
//
//  P0.resize(GDIM);
//  dX.resize(GDIM);
//  delta.resize(GDIM);
//  xnodes = &grid.xNodes();
//  glob_indices.resize(grid_->ndof());
//
//  // Cannot use utmp_ here because it
//  // may not have been allocated yet, so
//  // use local variable:
//  tmp.resize(GDIM);
//  for ( auto j=0; j<GDIM; j++ ) tmp[j] = new GTVector<GFTYPE>(grid_->ndof());
//
//  // If (x, y, z) < epsilon, set to 0:
//  GMTK::zero(*xnodes);
//
//  if ( sgrid == "grid_box" ) {
//    pstd   = gtree.getArray<GFTYPE>("xyz0");
//    P0     = pstd;
//    xkey.resize(GDIM);
//    for ( auto j=0; j<GDIM; j++ ) xkey[j] = &(*xnodes)[j];
//    dX     = 0.05*grid.minnodedist();
////  gmorton.setDoLog(TRUE);
//    gmorton.setType(GMORTON_STACKED);
//  }
//  if ( sgrid == "grid_icos" ) {
//#if 1
//    // Set indices from lat/lon coords:
//    P0.resize(GDIM);
//    dX.resize(GDIM);
//    cart.resize(xnodes->size());
//    xkey.resize(GDIM);
//    P0.x1 = 0.0 ; // lat starting point
//    P0.x2 = 0.0 ; // lon starting point
//    for ( auto j=0; j<xnodes->size(); j++ ) cart[j] = &(*xnodes)[j];
//    for ( auto j=0; j<GDIM; j++ ) xkey[j] = tmp[j];
//    GMTK::cart2latlon(cart, xkey);
//    for ( auto j=0; j<xkey[0]->size(); j++ ) (*xkey[0])[j] = 0.5*PI - (*xkey[0])[j];
//    for ( auto j=0; j<xkey[1]->size(); j++ ) (*xkey[1])[j] += (*xkey[1])[j] < 0.0 ? 2.0*PI : 0.0;
//    rad   = gtree.getValue<GFTYPE>("radius");
//    delta[0] = 0.5*grid.minlength()/(rad*pmax*pmax);
//    delta[1] = grid.minlength()/(rad*pmax*pmax);
//    for ( auto j=0; j<dX.size(); j++ ) dX[j] = 0.025 *delta[j];
////  gmorton.setType(GMORTON_STACKED);
////  gmorton.setType(GMORTON_INTERLEAVE);
//#else
//    P0.resize(GDIM+1);
//    dX.resize(GDIM+1);
//    xkey.resize(GDIM+1);
//    rad   = gtree.getValue<GFTYPE>("radius");
//    P0.x1 = -rad-tiny; P0.x2 = -rad-tiny; P0.x3 = -rad-tiny;
//    for ( auto j=0; j<xkey.size(); j++ ) xkey[j] = tmp[j];
//    delta[0] = grid.minlength()/(pmax*pmax);
//    delta[1] = grid.minlength()/(pmax*pmax);
//    delta[2] = grid.minlength()/(pmax*pmax);
//    for ( auto j=0; j<dX.size(); j++ ) dX[j] = 0.01 *delta[j];
////  dX     = 0.1*grid.minnodedist();
////  gmorton.setDoLog(TRUE);
//    gmorton.setType(GMORTON_INTERLEAVE);
////  gmorton.setType(GMORTON_STACKED);
//#endif
//  }
//  if ( sgrid == "grid_sphere" ) {
//    rad   = gtree.getValue<GFTYPE>("radiusi");
//    P0.x1 = 0.0; // lat starting point
//    P0.x2 = 0.0; // long starting point
//    P0.x3 = rad; // radius starting point
//    cart.resize(xnodes->size());
//    xkey.resize(GDIM);
//    for ( auto j=0; j<xnodes->size(); j++ ) cart[j] = &(*xnodes)[j];
//    for ( auto j=0; j<GDIM; j++ ) xkey[j] = tmp[j];
//    GMTK::rcart2sphere(cart, xkey);
//    for ( auto j=0; j<GDIM; j++ ) ldelta[j] = xkey[j]->amindiff(tiny);
//    GComm::Allreduce(ldelta.data(), delta.data(), GDIM, T2GCDatatype<GFTYPE>(), GC_OP_MIN, comm_);
//    for ( auto j=0; j<GDIM; j++ ) dX[j] = 0.025*delta[j];
////  gmorton.setDoLog(TRUE);
//    gmorton.setType(GMORTON_STACKED);
//  }
//
//  // First, periodize coords if required to,
//  // before labeling nodes:
//  if ( typeid(grid) == typeid(GGridBox) ) {
//    static_cast<GGridBox*>(&grid)->periodize();
//  }
//
//  xnodes = &grid.xNodes();
//  glob_indices.resize(grid.ndof());
//
//  // Integralize *all* internal nodes
//  // using Morton indices:
////gmorton.setDoLog(TRUE);
//  gmorton.setType(GMORTON_STACKED);
////gmorton.setType(GMORTON_INTERLEAVE);
//  gmorton.setIntegralLen(P0,dX);
//  gmorton.key(glob_indices, xkey);
//
//  // Initialize gather/scatter operator:
//  ggfx = new GGFX<GFTYPE>();
//  assert(ggfx != NULLPTR && "Cannot instantiate GGFX operator");
//  bret = ggfx->init(glob_indices);
//  assert(bret && "Initialization of GGFX operator failed");
//
//  // Unperiodize nodes now that connectivity map is
//  // generated, so that coordinates mean what they should:
//  if ( typeid(grid) == typeid(GGridBox) ) { // periodize coords
//    static_cast<GGridBox*>(&grid)->unperiodize();
//  }
//
//  for ( auto j=0; j<GDIM; j++ ) delete tmp[j];
//
//} // end method init_ggfx


//**********************************************************************************
//**********************************************************************************
// METHOD: compare
// DESC  : Top-level method to do a comparison of
//         integrated solution with specified initial 
//         solution and write metrics to a file
// ARGS  : ptree: main prop tree
//         grid : grid object
//         peqn : pointer to EqnBase 
//         t    : current time
//         utmp : vector of tmp vectors
//         u    : full state vector
//**********************************************************************************
void compare(const PropertyTree &ptree, GGrid &grid, EqnBasePtr &peqn, Time &t, State &utmp, State &ub, State &u)
{
	GEOFLOW_TRACE();
  GBOOL             bret, bvardt;
  GINT              myrank, ntasks;
  GFTYPE            dxmin, lmin, tt;
  GTVector<GFTYPE>  lnorm(3), gnorm(3), maxerror(3);
  GTVector<GFTYPE>  nnorm(nsolve_);
  GTVector<GString> savars, scvars, sdvars;
  State             ua(nstate_);
  std::vector<GINT> pstd(GDIM);
  GString           sdir = ".";
  GString           sfile;
  char              stmp[1024];
  PropertyTree      vtree = ptree.getPropertyTree("stepper_props");
  
  bret   = ptree.getValue<GBOOL>("do_comparison");
  if ( !bret ) return;

  ntasks = GComm::WorldSize(comm_);
  myrank = GComm::WorldRank(comm_);
  bvardt = vtree.getValue<GBOOL>("variable_dt",FALSE);
  sfile  = ptree.getValue<GString>("compare_file","compare.txt");
  pstd   = ptree.getArray<GINT>("exp_order");

  // Create analytic solution array:
  for ( GINT j=0; j<ua.size(); j++ ) ua[j] = new GTVector<GFTYPE>(grid.ndof());

  // Set up some output variables:
  for ( GSIZET j=0; j<u.size(); j++ ) {
    sprintf(stmp, "u%lua", j+1);
    savars.push_back(stmp);
    sprintf(stmp, "diff%lu", j+1);
    sdvars.push_back(stmp);
  }
  for ( GSIZET j=0; j<c_.size(); j++ ) {
    sprintf(stmp, "c%lu", j+1);
    scvars.push_back(stmp);
  }


  // Compute analytic solution, do comparisons:
    
  maxerror = 0.0;
  lnorm    = 0.0;  
  nnorm    = 1.0;


  tt = t;
  bret = GInitStateFactory<MyTypes>::init(ptree, grid, peqn->stateinfo(), tt, utmp, ub, ua);
  assert(bret && "state initialization failed");
  for ( GSIZET j=0; j<nsolve_; j++ ) { // local errors
   *utmp [1] = *ua [j]; utmp [1]->rpow(2);
    nnorm[j] = grid.integrate(*utmp   [1],*utmp [0]) ; // L2 norm of analyt soln at t=0
    nnorm[j] = nnorm[j] > std::numeric_limits<GFTYPE>::epsilon() ? nnorm[j] : 1.0;
    cout << "main: nnorm[" << j << "]=" << nnorm[j] << endl;
  }
    
  GTVector<GINT> istate(nsolve_);
  GTVector<GINT> cstate(c_.size());
  for ( GINT j=0; j<nsolve_; j++ ) istate[j] = j;
  for ( GINT j=0; j<c_.size(); j++ ) cstate[j] = j;

  // Compute analytic solution at t:
  tt = t;
  bret = GInitStateFactory<MyTypes>::init(ptree, grid, peqn->stateinfo(), tt, utmp, ub, ua);
  assert(bret && "state initialization failed");

#if 0
  // Set up and output the analytic solution
  // and difference solution as well as
  // advection velocity if one exists:
  std::stringstream  format;
  StateInfo          stateinfo;
  ObsTraitsType      binobstraits = (*pObservers_)[irestobs_]->get_traits();

  stateinfo.sttype   = 1; // state variable type
  stateinfo.svars.resize(binobstraits.state_names.size());
  stateinfo.svars    = binobstraits.state_names;
  stateinfo.idir     = binobstraits.idir;
  stateinfo.odir     = binobstraits.odir;
  stateinfo.index    = 0;
  stateinfo.cycle    = 0;
  stateinfo.time     = t;

  // Get data:
  stateinfo.svars.resize(nsolve_);
  for ( GINT j=0; j<nsolve_; j++ ) { 
    format .str(""); format .clear();
    format << "u" << j+1 << "a";
    stateinfo.svars[j] = format.str();
  }
  pIO_->write_state("ua", stateinfo, ua);

  for ( GINT j=0; j<nsolve_; j++ ) { 
    *utmp[j] = *u[j] - *ua[j];
    format .str(""); format .clear();
    format << "diff" << j+1;
    stateinfo.svars[j] = format.str();
  }
  pIO_->write_state("diff", stateinfo, utmp);
  
  stateinfo.svars.resize(c_.size());
  for ( GINT j=0; j<c_.size(); j++ ) { 
    format .str(""); format .clear();
    format << "c" << j+1;
    stateinfo.svars[j] = format.str();
  }

  pIO_->write_state("c", stateinfo, c_);
#endif

  // Compute error norms:
  for ( GINT j=0; j<nsolve_; j++ ) { //local errors
   *utmp [0] = *u[j] - *ua[j];
   *utmp [1]  = *utmp[0]; utmp [1]->abs();
   *utmp [2]  = *utmp[0]; utmp [2]->rpow(2);
    lnorm[0]  = utmp [0]->infnorm (); // inf-norm
    gnorm[1]  = grid.integrate(*utmp[1],*utmp[0]); // L1-norm numerator
    gnorm[2]  = grid.integrate(*utmp[2],*utmp[0]); // L2-norm numerator
    // Accumulate to find global errors for this field:
    GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm_);
    gnorm[1] =  gnorm[1]/nnorm[j];
    gnorm[2] =  sqrt(gnorm[2]/nnorm[j]);
    // now find max errors of each type for each field:
    for ( GINT i=0; i<3; i++ ) maxerror[i] = MAX(maxerror[i],fabs(gnorm[i]));
  }

  // Compute some global quantities for output:
  dxmin = grid.minnodedist();
  lmin  = grid.minlength();
  if ( myrank == 0 ) {
    cout << "main: maxerror = " << maxerror << endl;
  }
   
  GTVector<GSIZET> lsz(2), gsz(2);
  lsz[0] = grid.nelems();
  lsz[1] = grid.ndof();
  GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);

  // Print convergence data to file:
  std::ifstream itst;
  std::ofstream ios;

  if ( myrank == 0 ) {
    itst.open(sfile);
    ios.open(sfile,std::ios_base::app);
  
    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
      ios << "#ntasks" << "  ";
      ios << "ncyc"    << "  ";
      ios << "var_dt"  << "  ";
      for ( GSIZET j=0; j<GDIM; j++ ) ios << "p" << j+1 << "  ";
      ios << "num_elems    dx_min   EL_min     inf_err     L1_err      L2_err" << std::endl;
    }
    itst.close();

    ios << ntasks << "  " ;
    ios << pIntegrator_->get_numsteps()  << "  ";
    ios << bvardt << "  ";
    for ( GINT j=0; j<GDIM; j++ ) ios << pstd[j] << "  ";
    ios << gsz[0] << "  " << dxmin       << "  " << lmin
                  << "  " << maxerror[0] << "  " << maxerror[1] 
                  << "  " << maxerror[2]
        << std::endl;
    ios.close();
  }

  for ( GINT j=0; j<ua.size(); j++ ) delete ua[j];

} // end method compare


//**********************************************************************************
//**********************************************************************************
// METHOD : do_terrain
// DESC   : Set the user-specified terrain in grid
// ARGS   : 
//          ptree: main prop tree
//          grid : grid object
// RETURNS: none.
//**********************************************************************************
void do_terrain(const PropertyTree &ptree, GGrid &grid)
{
	GEOFLOW_TRACE();
  GBOOL bret, bterr;
  GINT  iret, nc;
  State xb, tmp, utmp;
  
  nc = grid.xNodes().size();
  xb.resize(nc);

  // Cannot use utmp_ here because it
  // may not have been allocated yet, so
  // use local variable:
  utmp.resize(2*nc);
  tmp.resize(nc);
  for ( auto j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(grid_->ndof());

  // Set terrain & tmp arrays from tmp array pool:
  for ( auto j=0; j<nc        ; j++ ) xb [j] = utmp[j];
  for ( auto j=0; j<tmp.size(); j++ ) tmp[j] = utmp[j+nc];

  bret = GSpecTerrainFactory<MyTypes>::spec(ptree, grid, tmp, xb, bterr);
  assert(bret);

  if ( bterr ) grid.add_terrain(xb, tmp);

  for ( auto j=0; j<utmp.size(); j++ ) delete utmp[j];

} // end of method do_terrain


//**********************************************************************************
//**********************************************************************************
// METHOD : do_restart
// DESC   : Set state for restart
// ARGS   : 
//          ptree: main prop tree
//          grid : grid object
//          u    : read-in state
//          cycle: time cycle from restart
//          t    : time from restart
// RETURNS: none.
//**********************************************************************************
void do_restart(const PropertyTree &ptree, GGrid &, State &u, 
                GTMatrix<GINT>&p,  GSIZET &cycle, Time &t)
{
	GEOFLOW_TRACE();
  assert(pIO_ != NULLPTR && "IO operator not set!");

  GBOOL              bret;
  GSIZET             itindex;
  GFTYPE             tt = t;
  std::stringstream  format;
  StateInfo          stateinfo;
  ObsTraitsType      binobstraits = (*pObservers_)[irestobs_]->get_traits();

  itindex            = ptree.getValue<GSIZET>("restart_index", 0);
  stateinfo.icomptype.resize(pEqn_->stateinfo().icomptype.size());
  stateinfo.icomptype= pEqn_->stateinfo().icomptype;
  stateinfo.sttype   = 0; // state variable type
  stateinfo.svars.resize(binobstraits.state_names.size());
  stateinfo.svars    = binobstraits.state_names;
  stateinfo.idir     = binobstraits.idir;
  stateinfo.odir     = binobstraits.odir;
  stateinfo.index    = itindex;

  // Get data:
  pIO_->read_state(binobstraits.agg_state_name, stateinfo, u);
  p.resize(stateinfo.porder.size(1),stateinfo.porder.size(2));
  p = stateinfo.porder;
 
  // Assign time and cycle from traits, for return:
  cycle = stateinfo.cycle;
  t     = stateinfo.time;

  
} // end of method do_restart






void init_ggfx(PropertyTree& ptree, GGrid& grid, GGFX<GFTYPE>*& ggfx)
{
  GEOFLOW_TRACE();

  std::vector<std::array<GFTYPE,GDIM>> xyz(grid_->ndof());
  GEOFLOW_TRACE_START("ReOrder XYZ Arrays");
  const auto ndof = grid_->ndof();
  for(std::size_t i = 0; i < ndof; i++){
	  for(std::size_t d = 0; d < GDIM; d++){
		  xyz[i][d] = grid.xNodes()[d][i];
	  }
  }
  GEOFLOW_TRACE_STOP();

  // Create GGFX
  ASSERT(ggfx == nullptr);
  ggfx = new GGFX<GFTYPE>();
  ASSERT(ggfx != nullptr);
  pio::pout << "Calling ggfx->init(xyz)" << std::endl;
  ggfx->init(0.25*grid.minnodedist(), xyz);

} // end method init_ggfx_2









