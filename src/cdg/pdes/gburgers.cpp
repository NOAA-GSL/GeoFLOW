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
// ARGS   : grid     : grid object
//          u        : state (i.e., vector of GVectors)
//          isteptype: stepper type
//          iorder   : vector of integers indicating the 
//                       index 0: time order for du/dt derivaitve for multistep
//                                methods; RK order (= num stages + 1)
//                       index 1: order of approximation for nonlin term
//          tmp      : Array of tmp vector pointers, pointing to vectors
//                     of same size as State
// RETURNS: none
//**********************************************************************************
GBurgers::GBurgers(GGrid &grid, State &u, GStepperType isteptype, GTVector<GINT> iorder, GTVector<GTVectorGFTYPE>*> &tmp) :
bconserved_     (FALSE),
isteptype_  (isteptype),
nsteps_             (0),
itorder_            (0),
inorder_            (0),
gadvect_      (NULLPTR),
gmass_        (NULLPTR),
gpdv_         (NULLPTR),
//gflux_      (NULLPTR),
grid_           (&grid)
{
  static_assert(std::is_same<State,GTVector<GTVectorGFTYPE>>>::value,
               "State is of incorrect type"); 
  valid_types_.resize(4);
  valid_types_[0] = GSTEPPER_RK2;
  valid_types_[1] = GSTEPPER_RK4;
  valid_types_[2] = GSTEPPER_BDFAB;
  valid_types_[3] = GSTEPPER_BDFEXT;
 
  assert(valid_types_.contains(isteptype) && "Invalid stepper type"); 

  // Find multistep/multistage coefficients:
  GMultilevel_coeffs_base *tcoeff_obj; // time deriv coeffs
  GMultilevel_coeffs_base *acoeff_obj; // adv op. coeffs
  switch ( isteptype_ ) {
    case GSTEPPER_RK2:
      itorder_   = iorder[0];
      break;
    case GSTEPPER_RK4:
      itorder_   = iorder[0];
      break;
    case GSTEPPER_BDFAB:
      itorder_   = iorder[0];
      inorder_   = iorder[1];
      tcoeff_obj = new G_BDF(itorder_);
      acoeff_obj = new G_AB (inorder_);
      tcoeffs_.resize(tcoeff_obj.getCoeffs().size());
      acoeffs_.resize(acoeff_obj.getCoeffs().size());
      tcoeffs_ = tcoeff_obj;
      acoeffs_ = acoeff_obj;
      break;
    case GSTEPPER_BDFEXT:
      itorder_   = iorder[0];
      inorder_   = iorder[1];
      tcoeff_obj = new G_BDF(itorder_);
      acoeff_obj = new G_EXT(inorder_);
      tcoeffs_.resize(tcoeff_obj.getCoeffs().size());
      acoeffs_.resize(acoeff_obj.getCoeffs().size());
      tcoeffs_ = tcoeff_obj;
      acoeffs_ = acoeff_obj;
      break;
  }
  delete tcoeff_obj;
  delete acoeff_obj'
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
// METHOD : dudt_impl
// DESC   : Compute RHS for explicit schemes
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers::dudt_impl(const Time &t, State &u, Time &dt, Derivative &dudt)
{
  assert(!bconserved_ &&
         "conservation not yet supported"); 

  // If non-conservative, compute RHS from:
  //     du/dt = -u.Grad u + nu nabla u 
  // for each u
  
  for ( auto k=0; k<u.size(); k++ ) {
    gadvect_->apply(u[k], u, utmp_, dudt[k]);
    ghelm_->opVec_prod(u[k],utmp_[0]);
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_extbdf
// DESC   : EXT/BDF time stepping: apply EXT to advection terms, and BDF to time
//          derivative
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers::step_extbdf(const Time &t, State &uin, Time &dt, Derivative &uout)
{

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


