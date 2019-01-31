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
  valid_types_[0] = GSTEPPER_EXRK2;
  valid_types_[1] = GSTEPPER_EXRK4;
  valid_types_[2] = GSTEPPER_BDFAB;
  valid_types_[3] = GSTEPPER_BDFEXT;
 
  assert(valid_types_.contains(isteptype_) && "Invalid stepper type"); 

  init(u, iorder);
  
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

  // Make sure that, in init(), Helmholtz op is using only
  // weak Laplacian, or there will be problems:
  for ( auto k=0; k<u.size(); k++ ) {
    gadvect_->apply(u[k], u, utmp_, dudt[k]);
    ghelm_->opVec_prod(u[k],utmp_[0]);
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl
// DESC   : Step implementation method
// ARGS   : t   : time
//          u   : state
//          dt  : time step
//          uout: updated state
// RETURNS: none.
//**********************************************************************************
void GBurgers::step_extbdf(const Time &t, State &uin, Time &dt, Derivative &uout)
{

  switch ( isteptype_ ) {
    case GSTEPPER_EXRK2:
    case GSTEPPER_EXRK4:
      step_exrk(t, uin, dt, uout);
      break;
    case GSTEPPER_BDFAB:
    case GSTEPPER_BDFEXT:
      step_imex(t, uin, dt, uout);
      break;
  }
  
} // end of method step_extbdf


//**********************************************************************************
//**********************************************************************************
// METHOD : step_multistep
// DESC   : Carries out multistep update. The time derivative and 
//          advection terms are handlex using a multistep expansion. The
//          advection term is 'extrapolated' to the new time level
//          using known state data so as to obviate the need for a fully 
//          implicit treatment. The dissipation term is handled implicitly.
// ARGS   : t   : time
//          u   : state
//          dt  : time step
//          uout: updated state
// RETURNS: none.
//**********************************************************************************
void GBurgers::step_multistep(const Time &t, State &uin, Time &dt, Derivative &uout)
{

  // If non-conservative, compute RHS from:
  //     du/dt = -u.Grad u + nu nabla u 
  // for each u
  
  for ( auto k=0; k<u.size(); k++ ) {
    gadvect_->apply(u[k], u, utmp_, dudt[k]);
    ghelm_->opVec_prod(u[k],utmp_[0]);
  }
  
} // end of method step_multistep


//**********************************************************************************
//**********************************************************************************
// METHOD : step_exrk
// DESC   : Take a step using Explicit RK method
// ARGS   : t   : time
//          u   : state
//          dt  : time step
//          uout: updated state
// RETURNS: none.
//**********************************************************************************
void GBurgers::step_exrk(const Time &t, State &uin, Time &dt, Derivative &uout)
{

  


} // end of method step_exrk



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
// ARGS   : u     : State variable
//          iorder: time stepping trunction order vector. If using an explicit 
//                  scheme, only the first vector member is used. If using
//                  semi-implicit schemes, the first slot is for the time 
//                  derivative order, and the second for the nonlinear 
//                  term 'extrapolation' order. 
// RETURNS: none.
//**********************************************************************************
void GBurgers::init(State &u, GTVector<GINT> &iorder)
{
  GString serr = "GBurgers::init: ";

  // Find multistep/multistage time stepping coefficients:
  GMultilevel_coeffs_base *tcoeff_obj; // time deriv coeffs
  GMultilevel_coeffs_base *acoeff_obj; // adv op. coeffs
  switch ( isteptype_ ) {
    case GSTEPPER_EXRK2:
    case GSTEPPER_EXRK4:
      gextk_.setOrder(itorder_);
      gextk_.setRHSfunction(this->dudt);
      break;
    case GSTEPPER_BDFAB:
      itorder_   = iorder[0];
      inorder_   = iorder[1];
      tcoeff_obj = new G_BDF(itorder_);
      acoeff_obj = new G_AB (inorder_);
      tcoeffs_.resize(tcoeff_obj.getCoeffs().size());
      acoeffs_.resize(acoeff_obj.getCoeffs().size());
      tcoeffs_ = tcoeff_obj; acoeffs_ = acoeff_obj;
      break;
    case GSTEPPER_BDFEXT:
      itorder_   = iorder[0];
      inorder_   = iorder[1];
      tcoeff_obj = new G_BDF(itorder_);
      acoeff_obj = new G_EXT(inorder_);
      tcoeffs_.resize(tcoeff_obj.getCoeffs().size());
      acoeffs_.resize(acoeff_obj.getCoeffs().size());
      tcoeffs_ = tcoeff_obj; acoeffs_ = acoeff_obj;
      break;
  }
  delete tcoeff_obj;
  delete acoeff_obj'
  
  // Instantiate spatial discretization operators:
  padvect_ = new GAdvect(grid_);
  gmass_   = new GMass(*grid_);
  ghelm_   = new GHelmholtz(*grid_);
  
  // If doing semi-implicit time stepping; handle viscous term 
  // (linear) inplicitly, which implies using full Helmholtz operator:
  if ( isteptype_ == GSTEPPER_BDFAB || isteptype_ == GSTEPPER_BDFEXT ) {
    ghelm_->set_mass(*gmass_);
  }

  if ( bconserved_ ) {
    assert(FALSE && "Conservation not yet supported");
    gpdv_  = new GpdV(*grid_,*gmass_);
//  gflux_ = new GFlux(*grid_);
    assert( (gmass_   != NULLPTR
          || ghelm_   != NULLPTR
          || gpdv_    != NULLPTR) && "1 or more operators undefined");
  }
  else {
    gadvect_ = new GAdvect(*grid_);
    assert( (gmass_   != NULLPTR
          || ghelm_   != NULLPTR
          || gaevect_ != NULLPTR) && "1 or more operators undefined");
  }

  // If doing a multi-step method, instantiate (deep) space for 
  // required time levels for state:
  u_keep_ .resize(itorder_);
  for ( auto i=0; i<itorder_-1; i++ ) { // for each time level
    for ( auto j=0; j<u.size(); j++ ) u_keep_.resize(u[j].size());
  }

} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD : cycle_keep
// DESC   : Cycle the mult-level states making sure the most
//          recent is at index 0, the next most recent, at index 1, etc...
// ARGS   : u     : State variable providing most recent state
// RETURNS: none.
//**********************************************************************************
void GBurgers::cycle_keep(State &u)
{

  // Make sure following index map contains the 
  // correct time level information:
  //   u_keep[0] <--> time level n (most recent)
  //   u_keep[1] <--> time level n-1
  //   u_keep[2] <--> time level n-2 ...
  u_keep_ .resize(itorder_);
  for ( auto i=itorder_-1; i>=1; i-- ) u_keep[i] = u_keep[i+1];
  u_keep[0] = u;

} // end of method cycle_keep


