//==================================================================================
// Module       : gtest_burgers.cpp
// Date         : 12/24/18 (DLR)
// Description  : GeoFLOW test of (linear-only) advection
// Copyright    : Copyright 19. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gtvector.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "gadvect.hpp"
#include "ghelmholtz.hpp"
#include "gmass.hpp"
#include "gexrk_stepper.hpp"
#include "ggrid_factory.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "gio.h"
#include "gtools.h"

#if defined(_G_USE_GPTL)
    #include "gptl.h"
#endif

using namespace geoflow::tbox;
using namespace std;

typedef GTVector<GTVector<GFTYPE>*> State;
typedef GTVector<GTVector<GFTYPE>*> Derivative;
typedef GFTYPE Time;


GGrid      *grid_;
GAdvect    *gadvect_;
GHelmholtz *ghelm_;
GMass      *gimass_;
GGFX       *ggfx_;
State u_;
State c_;
State ua_;
State ub_;
State utmp_;
State urktmp_;
State urhstmp_;
State uoptmp_;
State uold_;
GTVector<GFTYPE> nu_(3);
PropertyTree ptree;       // main prop tree


void compute_analytic(GGrid &grid, GFTYPE &t, const PropertyTree& ptree, GTVector<GTVector<GFTYPE>*> &u);
void update_dirichlet(const Time &t, State &u, State &ub);
void apply_bc(const Time &t, State &u, const State &ub);
void dudt(const Time &t, const State &u,
          const Time &dt, State &dudt);


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GBOOL  bgridwritten=FALSE;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   n, np=1;    // elem 'order'
    GINT   nstate=GDIM;  // number 'state' arrays
    GINT   nsolve=GDIM;  // number *solved* 'state' arrays
    GSIZET maxSteps;
    GFTYPE radiusi=1, radiuso=2, dt, t, tt;
    std::vector<GINT> ne(3); // # elements in each direction in 3d
    GString sgrid;// name of JSON grid object to use
    GTVector<GString> svars, savars, sdu;
    char stmp[1024];
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm & global environment:
//  GComm::InitComm(&argc,&argv);
    mpixx::environment env(argc, argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc, argv); 
    GlobalManager::startup();
    comm = world; // need this for solver(s) & grid


    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree eqptree;     // equation props
    PropertyTree gridptree;   // grid props
    PropertyTree stepptree;   // stepper props
    PropertyTree dissptree;   // dissipation props
    PropertyTree tintptree;   // time integration props

cout << "main: call load_file..." << endl;
    ptree    = InputManager::getInputPropertyTree();       // main param file structure
cout << "main: load_file done." << endl;
    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    np          = ptree.getValue<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);
    stepptree   = ptree.getPropertyTree("stepper_props");
    dissptree   = ptree.getPropertyTree("dissipation_traits");
    tintptree   = ptree.getPropertyTree("time_integration");

    ne          = gridptree.getArray<GINT>("num_elems");  // may be modified by command line


    // Parse command line. ':' after char
    // option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "i:j:k:l:m:p:h")) != -1) {
      switch (iopt) {
      case 'i': // handled by InputManager
          break;
      case 'j': // get # elements in r/x
          ne[0] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'k': // get # elements in lat/y
          ne[1] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'm': // get # elements in long/z
          ne[2] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          gridptree.setValue<GINT>("ilevel",ilevel);
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
          ptree.setValue<GINT>("exp_order",np);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-j #Elems in x/r] [-k #Elems in y/lat]  [-m #Elems in z/long] [-l refine level] -p expansion order] " << std::endl;
          std::cout << "Note: Icos grid only requires refine level to specify number of elements. " << std::endl;
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


    // Set solver traits from prop tree:
    GINT   itorder    = stepptree.getValue <GINT>("time_deriv_order");
    GFTYPE nu_scalar  = dissptree.getValue<GFTYPE>("nu");
    nu_.resize(1); 
    nu_ = nu_scalar; 

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
    grid_ = GGridFactory::build(gridptree, gbasis, comm);
    GPTLstop("gen_grid");


    GPTLstart("do_gather_op");
    // Initialize gather/scatter operator:
    ggfx_ = new GGFX();
    init_ggfx(ptree, *grid_, *ggfx_);
    GPTLstop("do_gather_op");


    // Create state and tmp space:
    nstate = nsolve = 1; // 1-state + GDIM  v-components
    utmp_.resize(24);
    uold_.resize(nsolve);
    u_.resize(nsolve);
    ua_.resize(nsolve);
    ub_.resize(nsolve);
    c_ .resize(GDIM);
    urktmp_ .resize(nstate*(itorder+2)+1); // RK stepping work space
    urhstmp_.resize(1); // work space for RHS
    GSIZET nop = utmp_.size()-uold_.size()-urktmp_.size()-urhstmp_.size();
    assert(nop > 0 && "Invalid operation tmp array specification");
    uoptmp_ .resize(nop); // RHS operator work space

    for ( GSIZET j=0; j<uold_.size(); j++ ) uold_[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<utmp_.size(); j++ ) utmp_[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<u_   .size(); j++ ) u_   [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<ua_  .size(); j++ ) ua_  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<ub_  .size(); j++ ) ub_  [j] = new GTVector<GFTYPE>(grid_->nbdydof());
    for ( GSIZET j=0; j<c_   .size(); j++ ) c_   [j] = new GTVector<GFTYPE>(grid_->size());
    n = 0;
    for ( GSIZET j=0; j<urktmp_ .size(); j++, n++ ) urktmp_ [j] = utmp_[n];
    for ( GSIZET j=0; j<urhstmp_.size(); j++, n++ ) urhstmp_[j] = utmp_[n];
    for ( GSIZET j=0; j<uoptmp_ .size(); j++, n++ ) uoptmp_ [j] = utmp_[n];


    // Create RK stepper, rhs objects:
    GExRKStepper<GFTYPE> gexrk(itorder);

    gadvect_ = new GAdvect(*grid_);
    ghelm_   = new GHelmholtz(*grid_);
    gimass_  = new GMass(*grid_, TRUE);

    // Set Dirichlet bdy state update function:
    std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, 
                                        GTVector<GTVector<GFTYPE>*> &ub)>  
        updatebc = update_dirichlet; // set tmp function with proper signature for...

    // Set bdy condition application function:
    std::function<void(const Time &t,         
                       State  &uin,
                       const State &ub)> applybc = apply_bc;

    // Set RHS function:
   std::function<void(const Time &t,   
                      const State  &uin,
                      const Time &dt,
                      State &dudt)> rhs = dudt;

    gexrk.set_update_bdy_callback(updatebc);
    gexrk.set_apply_bdy_callback(applybc);
    gexrk.setRHSfunction(rhs);
    gexrk.set_ggfx(ggfx_);

    dt         = tintptree.getValue<GFTYPE>("dt"); 
    maxSteps = tintptree.getValue<GSIZET>("cycle_end");

    // Initialize state:
    t = 0.0;
    compute_analytic(*grid_, t, ptree, u_);
    for ( auto j=0; j<u_.size(); j++ ) {
      sprintf(stmp, "u%d", j+1);
      svars.push_back(stmp);
      sprintf(stmp, "u%da", j+1);
      savars.push_back(stmp);
      sprintf(stmp, "du%d", j+1);
      sdu.push_back(stmp);
    }
    gio(*grid_, u_, 1, 0, 0.0, svars, comm, bgridwritten);

cout << "main: u(t=0)=" << *u_[0] << endl;
   // Set advection velocity (may be 
   // a functino of time):
   *c_[0] = 1.0;
   *c_[1] = 0.0;
   if ( GDIM == 3 ) *c_[2] = 0.0;

    GPTLstart("time_loop");
    for( GSIZET i=0; i<maxSteps; i++ ){
      gexrk.step(t, u_, ub_, dt, urktmp_);
      gio(*grid_, u_, 1, i, t, svars, comm, bgridwritten);
      t += dt;
    }
    GPTLstop("time_loop");
    
    tt = 0.0;
    compute_analytic(*grid_, tt, ptree, ua_); // analyt soln at t=0

    // Write binary output:
    gio(*grid_, u_, u_.size(), maxSteps, t, svars, comm, bgridwritten);

#if 1
    // Compute analytic solution, do comparisons:
    GTVector<GFTYPE> lnorm(3), gnorm(3), maxerror(3);
    GTVector<GFTYPE> nnorm(nsolve);
    
    maxerror = 0.0;
    lnorm    = 0.0;  
    nnorm    = 1.0;


    for ( GSIZET j=0; j<nsolve; j++ ) { //local errors
      ua_[j]->pow(2);
      nnorm[j] = grid_->integrate(*ua_  [j],*utmp_[0]) ; // L2 norm of analyt soln at t=0
cout << "main: ua(t=0)[" << j << "]_max=" << ua_[j]->max() << endl;
cout << "main: ua(t=0)[" << j << "]_imax=" << ua_[j]->imax() << endl;
cout << "main: u (t=" << t << ")[" << j << "]_max=" << u_[j]->max() << endl;
cout << "main: u (t=" << t << ")[" << j << "]_imax=" << u_[j]->imax() << endl;
cout << "main: nnorm[" << j << "]=" << nnorm[j] << endl;
    }
    
    compute_analytic(*grid_, t, ptree, ua_); // analyt soln at t
    gio(*grid_, ua_, ua_.size(),  maxSteps, t, savars, comm, bgridwritten);
    for ( GSIZET j=0; j<nsolve; j++ ) { //local errors
cout << "main: u [t=" << t << "]=" << *u_ [ j] << endl;
cout << "main: ua[t=" << t << "]=" << *ua_ [j] << endl;
cout << "main: imax(u)=" << u_[j]->imax() << " imax(ua)=" << ua_[j]->imax() << endl;
cout << "main: max (u)=" << u_[j]->max() << " max (ua)=" << ua_[j]->max() << endl;
      *utmp_[0] = *u_[j] - *ua_[j];
cout << "main: diff=u-ua[" << j << "]=" << *utmp_[0] << endl;
      *utmp_[1] = *utmp_[0]; utmp_[1]->abs();
      *utmp_[2] = *utmp_[0]; utmp_[2]->pow(2);
      lnorm   [0]  = utmp_[0]->infnorm (); // inf-norm
      gnorm   [1]  = grid_->integrate(*utmp_[1],*utmp_[0])/sqrt(nnorm[j]) ; // L1-norm
      gnorm   [2]  = grid_->integrate(*utmp_[2],*utmp_[0]) ; // L2-norm
      // Accumulate to find global errors for this field:
      GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
      gnorm[2] =  sqrt(gnorm[2]/nnorm[j]);
      // now find max errors of each type for each field:
      for ( GSIZET i=0; i<3; i++ ) maxerror[i] = MAX(maxerror[i],fabs(gnorm[i]));
    }

    cout << "main: maxerror = " << maxerror << endl;
   
    GSIZET lnelems=grid_->nelems();
    GSIZET gnelems;
    GComm::Allreduce(&lnelems, &gnelems, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm);

    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("burgers_err.txt");
    ios.open("burgers_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p     num_elems      inf_err     L1_err      L2_err" << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  "
        << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
        << std::endl;
    ios.close();

#endif
 
    GPTLpr_file("timing.txt");
    GPTLfinalize();

//  GComm::TermComm();
    GlobalManager::shutdown();
    GlobalManager::finalize();
    if ( grid_    != NULLPTR ) delete grid_;
    if ( gadvect_ != NULLPTR ) delete gadvect_;
    if ( ghelm_   != NULLPTR ) delete ghelm_;
    if ( ggfx_    != NULLPTR ) delete ggfx_;

    return(0);

} // end, main



//**********************************************************************************
//**********************************************************************************
// METHOD: compute_dirichlet_adv
// DESC  : Compute solution to pure advection equation with 
//         GBDY_DIRICHLET bcs, a Gaussian 'lump'. Must use box grid.
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_dirichlet_adv(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GBOOL            bContin;
  GINT             n;
  GFTYPE           argxp; 
  GFTYPE           nxy, sig0, u0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  // Check bdy conditioins:
  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_DIRICHLET") >= 2*GDIM 
      && "Dirichlet boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  r0.x1 = heatptree.getValue<GFTYPE>("x0"); 
  r0.x2 = heatptree.getValue<GFTYPE>("y0"); 
  r0.x3 = heatptree.getValue<GFTYPE>("z0"); 
  sig0  = heatptree.getValue<GFTYPE>("sigma"); 
  u0    = heatptree.getValue<GFTYPE>("u0"); 

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  *c_[0]  = 1.0;
  *c_[1]  = 0.0;
  if ( GDIM > 2 ) *c_[2]  = 0.0;

  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<GDIM; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 2.0*t*nu_[0]); // scalar viscosity only
    si   [k] = 0.5/(sig[k]*sig[k]);
    ufact[k] = u0*pow(sig0,2)/pow(sig[k],2);
  }

  // Ok, return to assumption of isotropic nu: 
  for ( GSIZET j=0; j<nxy; j++ ) {
    for ( GSIZET i=0; i<GDIM; i++ ) xx[i] = (*xnodes)[i][j] - r0[i] - (*c_[i])[j]*t;
    argxp = 0.0;
    for ( GSIZET i=0; i<GDIM; i++ ) argxp += -pow(xx[i],2.0)*si[i];
   (*ua[0])[j] = ufact[0]*exp(argxp);
  }
  
} // end, compute_dirichlet_adv


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_periodic_adv
// DESC  : Compute solution to pure advection equation with 
//         GBDY_PERIODIC bcs, a Gaussian 'lump'. Must use box grid.
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_periodic_adv(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GBOOL            bContin;
  GINT             n;
  GSIZET           i, j, k;
  GFTYPE           argxp, argxm, psum, rat, sum, eps;
  GFTYPE           nxy, pint, sig0, u0;
  GTVector<GFTYPE> f(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  eps = std::numeric_limits<GFTYPE>::epsilon();

  // Get periodicity length, gL:
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
  std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; r0 = dxyz; gL = r0;

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_PERIODIC") >= 2*GDIM
      && "Periodic boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  r0.x1 = heatptree.getValue<GFTYPE>("x0"); 
  r0.x2 = heatptree.getValue<GFTYPE>("y0"); 
  r0.x3 = heatptree.getValue<GFTYPE>("z0"); 
  sig0  = heatptree.getValue<GFTYPE>("sigma"); 
  u0    = heatptree.getValue<GFTYPE>("u0"); 

  // Set adv velocity components. Note:
  // First state elem is the scalar solution, and
  // the remainder are the velocity components:

  for ( j=0; j<GDIM; j++ ) {
    sig[j] = sqrt(sig0*sig0 + 2.0*t*nu_[0]);
    si [j] = 0.5/(sig[j]*sig[j]);
  }

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  *c_[0]  = 1.0;
  *c_[1]  = 0.0;
  if ( GDIM > 2 ) *c_[2]  = 0.0;

  for ( j=0; j<nxy; j++ ) {

    for ( k=0, argxp=0.0; k<GDIM; k++ ) {
      f [k]  = modf((*c_[k])[j]*t/gL[k],&pint);
//    f [k]  = (*c_[k])[j]*t/gL[k];
      xx[k]  = (*xnodes)[k][j] - r0[k] - f[k]*gL[k];
    }

    sum  = 0.0;

    n    = 0;
    rat = 1.0;
    while ( rat > eps ) {
      argxp = argxm = 0.0;
      for ( k=0; k<GDIM; k++ ) {
        argxp   += pow((xx[k]+n*gL[k]),2.0)*si[k];
        argxm   += pow((xx[k]-n*gL[k]),2.0)*si[k];
      }
      psum    =  n==0 ? exp(-argxp) : exp(-argxp) + exp(-argxm);
      sum    +=  psum;
      rat     = psum / sum ;
      n++;
    }
cout << "compute_periodic_adv: .................................n=" << n << endl;
    (*ua[0])[j] = u0*pow(sig0,2)/pow(sig[0],2)*sum;
  }

  
} // end, compute_periodic_adv


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_analytic
// DESC  : Compute analytic solutions based on property tree
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_analytic(GGrid &grid, GFTYPE &t, const PropertyTree& ptree, GTVector<GTVector<GFTYPE>*>  &ua)
{
  GString      sblock = ptree.getValue<GString>("init_block"); // name of initialization block
  GString      sbcs   = ptree.getValue<GString>("bdy_conditions"); // name of initialization block
  PropertyTree blockptree = ptree.getPropertyTree(sblock); // sub-block of main ptree describing initialization type

  
  if ( sbcs == "PERIODIC" ) {
  compute_periodic_adv(grid, t, ptree, ua);
  } else { // is DIRICHLET
  compute_dirichlet_adv(grid, t, ptree, ua);
  }

} // end, compute_analytic


//**********************************************************************************
//**********************************************************************************
// METHOD : update_dirichlet
// DESC   : update/set Dirichlet vectors, ub
// ARGS   : t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void update_dirichlet(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &ub)
{

  GFTYPE tt = t;
  GTVector<GTVector<GSIZET>> *igbdy = &grid_->igbdy();
  
  // If bc is time dependent, update it here.
  // Note: grid_ ptree, and ua_ are global:
  if ( (*igbdy)[GBDY_DIRICHLET].size() > 0 ) {
    compute_analytic(*grid_, tt, ptree,  ua_);
  }

  // ...GBDY_DIRICHLET:
  // Set from State vector, ua_:
  for ( GSIZET k=0; k<ua_.size(); k++ ) { 
    for ( GSIZET j=0; j<(*igbdy)[GBDY_DIRICHLET].size(); j++ ) {
      (*ub[k])[j] = (*ua_[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
  }

} // end of method update_dirichlet


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_bc
// DESC   : apply bcs
// ARGS   : t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void apply_bc(const Time &t, State &u, const State &ub)
{
  GTVector<GTVector<GSIZET>>  *igbdy = &grid_->igbdy();

  // Use indirection to set the global field node values
  // with domain boundary data. ub must be updated outside 
  // of this method.

  // NOTE: This is useful to set Dirichlet-type bcs only. 
  // Neumann bcs type have to be set with the
  // differential operators themselves, though the node
  // points in the operators may still be set from ub

  GBdyType itype;
  for ( GSIZET m=0; m<igbdy->size(); m++ ) { // for each type of bdy in gtypes.h
    itype = static_cast<GBdyType>(m);
    if ( itype == GBDY_NEUMANN
     ||  itype == GBDY_PERIODIC
     ||  itype == GBDY_OUTFLOW
     ||  itype == GBDY_SPONGE
     ||  itype == GBDY_NONE   ) continue;
    for ( GSIZET k=0; k<u.size(); k++ ) { // for each state component
      for ( GSIZET j=0; j<(*igbdy)[m].size(); j++ ) { // set Dirichlet-like value
        (*u[k])[(*igbdy)[m][j]] = (*ub[k])[j];
      }
    }
  }

} // end of method apply_bc


//**********************************************************************************
//**********************************************************************************
// METHOD: dudt
// DESC  : RHS function
// ARGS  : 
//**********************************************************************************
void dudt(const Time &t, const State &u,
          const Time &dt, State &dudt)
{
    gadvect_->apply   (*u[0], c_ , uoptmp_, *dudt[0]); // apply advection
//  ghelm_->opVec_prod(*u[0], uoptmp_, *urhstmp_[0]);  // apply diffusion
//  GMTK::saxpby<GFTYPE>(*urhstmp_[0], -1.0, *dudt[0], -1.0);
    *dudt[0] *= -1.0;
    gimass_->opVec_prod(*urhstmp_[0], uoptmp_, *dudt[0]); // apply M^-1

} // end method dudt

