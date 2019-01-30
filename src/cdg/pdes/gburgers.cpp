//==================================================================================
// Module       : gburgers.cpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE. 
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gburgers.hpp"



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method 
// DESC   : Instantiate with grid + state + tmp 
// ARGS   : grid  : grid object
//          u     : state (i.e., vector of GVectors)
//          tmp   : Array of tmp vector pointers, pointing to vectors
//                  of same size as State
// RETURNS: none
//**********************************************************************************
GBurgers::GBurgers(GGrid &grid, State &u, GTVector<GTVectorGFTYPE>*> &tmp) :
bconserved_ (FALSE),
gadvect_    (NULLPTR),
gmass_      (NULLPTR),
gpdv_       (NULLPTR),
//gflux_      (NULLPTR),
grid_       (&grid)
{
  static_assert(std::is_same<State,GTVector<GTVectorGFTYPE>>>::value,
               "State is of incorrect type"); 
  valid_types_.resize(3);
  valid_types_[0] = GSTEPPER_EXRK23;
  valid_types_[1] = GSTEPPER_ABBDF;
  valid_types_[2] = GSTEPPER_EXTBDF;
  init();
  
} // end of constructor method 



//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GBurgers::~GBurgers()
{
  if ( gmass_   != NULLPTR ) delete gmass_;
//if ( gflux_   != NULLPTR ) delete gflux_;
  if ( ghelm_   != NULLPTR ) delete ghelm_;
  if ( gadvect_ != NULLPTR ) delete gadvect_;
  if ( gpdv_    != NULLPTR ) delete gpdv_;
  ghelm_ = new GHelmholtz(*grid_);
  if ( bconsrved_ ) {
    gpdv_  = new GpdV(*grid_,*gmass_);
//  gflux_ = new GFlux(*grid_);
  }
  else {
    gadvect_ = new GAdvect(*grid_);
  }

} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : dt_impl
// DESC   : Compute time step, assuming a Courant number of 1:
//            dt = min_grid(dx/u)
// ARGS   : t : time
//          u : state
//          dt: timestep, returned
// RETURNS: none.
//**********************************************************************************
void GBurgers::dt_impl(const Time &t, State &u, Time &dt)
{
   GSIZET ibeg, iend;
   GFTYPE dtmin, umax;
   GTVector<GFTYPE> *drmin;

   drmin = &grid_->minnodedist();

   // This is an estimate. The minimum length on each element,
   // computed in GGrid object is divided by the maximum of
   // the state variable on each element:
   dt = 1.0;
   dtmin = std::numeric_limits<GFTYPE>::max();
   for ( auto k=0; k<u.size(); k++ ) { // each u
     for ( auto e=0; e<gelems_.size(); e++ ) { // for each element
       ibeg = gelems_[e]->igbeg(); iebd = gelems_[e]->igend();
       u[k].range(ibeg, iend);
       umax = u[k].max();
       dtmin = MIN(dtmin, (*drmin_)[e] / umax);
     }
     u[k].range_reset();
   }

   GComm::Allreduce(&dtmin, &dt, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm_);

} // end of method dt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_explicit
// DESC   : Compute RHS for explicit schemes
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers::dudt_impl(const Time &t, State &u, Time &dt, Derivative &dudt)
{
  assert(!bconserved_ &&
         "Jacobian computation not implemented"); 

  // If non-conservative, compute RHS from:
  //     du/dt = -u.Grad u + nu nabla u 
  // for each u
  
  for ( auto k=0; k<u.size(); k++ ) {
    gadvect_->apply(u[k], u, utmp_, dudt[k]);
    ghelm_->opVec_prod(u[k],utmp_[0]);
  }
  
} // end of method dudt_impl


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_callback
// DESC   : Set the callback object and method pointers for setting bdy conditions
// ARGS   : ptr2obj : pointer to object that hosts method callback
//          callback: method name
// RETURNS: none.
//**********************************************************************************
void GBurgers::set_bdy_callback(std::function<void(GGrid &)> &callback)
{
  bdycallback_  = &callback;
} // end of method set_bdy_callback

#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initialize equation object
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GBurgers::init()
{
  GString serr = "GBurgers::init: ";
  
  padvect_ = new GAdvect(grid_,

  gmass_ = new GMass(*grid_);
  ghelm_ = new GHelmholtz(*grid_);
  if ( bconsrved_ ) {
    gpdv_  = new GpdV(*grid_,*gmass_);
//  gflux_ = new GFlux(*grid_);
    assert( (gmass_   != NULLPTR
          || ghelm_   != NULLPTR
          || gpdv_    != NULLPTR) && "operators undefined");
  }
  else {
    gadvect_ = new GAdvect(*grid_);
    assert( (gmass_   != NULLPTR
          || ghelm_   != NULLPTR
          || gaevect_ != NULLPTR) && "operators undefined");
  }

} // end of method init


