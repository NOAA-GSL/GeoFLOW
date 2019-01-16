//==================================================================================
// Module       : gburgers_equation.ipp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE. 
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gburgers_equation.hpp"



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
GBurgers_equation::GBurgers_equation(GGrid &grid, State &u, GTVector<GTVectorGFTYPE>*> &tmp) :
bconserved_ (FALSE),
grid_      (&grid)
{
  static_assert(std::is_same<State,GTVector<GTVectorGFTYPE>>>::value,
               "State is of incorrect type"); 

  init();
  
} // end of constructor method 



//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GBurgers_equation::~GBurgers_equation()
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
void GBurgers_equation::dt_impl(const Time &t, State &u, Time &dt)
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
       ibeg = gelems_[e[->igbeg(); iebd = gelems_[e]->igend();
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
// METHOD : dudt_impl
// DESC   : Compute RHS for explicit schemes
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers_equation::dudt_impl(const Time &t, State &u, Time &dt, Derivative &dudt)
{
  assert(!bconserved_ &&
         "Jacobian computation not implemented"); 

  // If non-conservative, compute RHS from:
  //     du/dt = -u.Grad u + nu nabla u 
  // for each u
  
  for ( auto k=0; k<u.size(); k++ ) {
    gadvect_->apply(u[k], u, dudt[k]);
    ghelm_->opVec_prod(u[k],utmp_[0]);
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : dfdu_impl
// DESC   : Compute Jacobian dF/du for F == RHS. Used to linearize
//          equations for purely implicit schemes.
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers_equation::dfdu_impl(const Time &t, State &u, Time &dt, Jacobian &dfdu)
{
  assert(FALSE &&
         "Jacobian computation not implemented"); 
} // end of method dfdu_impl


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_callback
// DESC   : Set the callback object and method pointers for setting bdy conditions
// ARGS   : ptr2obj : pointer to object that hosts method callback
//          callback: method name
// RETURNS: none.
//**********************************************************************************
void GBurgers_equation::set_bdy_callback(std::function<void(GGrid &)> &callback)
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
void GBurgers_equation::init()
{
  GString serr = "GBurgers_equation::init: ";
  
  padvect_ = new GAdvect(grid_,

  gmass_ = new GMass(*grid_);
  ghelm_ = new GHelmholtz(*grid_);
  if ( bconsrved_ ) {
    gpdv_  = new GpdV(*grid_,*gmass_);
//  gflux_ = new GFlux(*grid_);
  }
  else {
    gadvect_ = new GAdvect(*grid_);
  }

} // end of method init


