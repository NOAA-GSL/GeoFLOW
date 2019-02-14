//==================================================================================
// Module       : gburgers.cpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE:
//                     du/dt + u . Del u = nu Del^2 u
//                This solver can be built in 2D or 3D, and can be configured to
//                remove the nonlinear terms so as to solve only the heat equation;
//
//                The State variable must always be of specific type
//                   GTVector<GTVector<GFTYPE>*>, but the elements rep-
//                resent different things depending on whether
//                the equation is doing nonlinear advection, heat only, or 
//                pure linear advection. If solving with nonlinear advection or
//                the heat equation, the State consists of elements [*u1, *u2, ....]
//                which is a vector solution. If solving the pure advection equation,
//                the State consists of [*u, *c1, *c2, *c3], where u is the solution
//                desired (there is only one, and it's a scalar), and ci 
//                are the constant-in-time Eulerian velocity components.
// 
// 
//                The dissipation coefficient may also be provided as a spatial field.
//
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gab.hpp"
#include "gext.hpp"
#include "gbdf.hpp"
#include "gburgers.hpp"



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method  (1)
// DESC   : Instantiate with grid + state + tmp. Use for fully nonlinear
//          Burgers equation, heat equation.
// ARGS   : ggfx      : gather/scatter operator
//          grid      : grid object
//          u         : state (i.e., vector of GVectors)
//          traits    :
//            steptype  : stepper type
//            itorder   : time order to du/dt derivative for multistep
//                        methods; RK num stages (~ order)
//            inorder   : order of approximation for nonlin term
//            doheat    : do heat equation only? If this is TRUE, then neither 
//                        of the following 2 flags have any meaning.
//            pureadv   : do pure advection? Has meaning only if doheat==FALSE
//            bconserved: do conservative form? Has meaning only if pureadv==FALSE.
//          tmp       : Array of tmp vector pointers, pointing to vectors
//                      of same size as State. Must be MAX(2*DIM+2,iorder+1)
//                      vectors
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GBurgers<TypePack>::GBurgers(GGFX &ggfx, GGrid &grid, State &u, GBurgers::Traits &traits, GTVector<GTVector<GFTYPE>*> &tmp) :
doheat_         (traits.doheat),
bpureadv_     (traits.bpureadv),
bconserved_ (traits.bconserved),
isteptype_    (traits.steptype),
nsteps_                     (0),
itorder_       (traits.itorder),
inorder_       (traits.inorder),
nu_                   (NULLPTR),
gadvect_              (NULLPTR),
gmass_                (NULLPTR),
gpdv_                 (NULLPTR),
//gflux_                (NULLPTR),
gbc_                  (NULLPTR),
grid_                   (&grid),
ggfx_                   (&ggfx),
comm_         (&ggrx.getComm())
{
  static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>>>::value,
               "State is of incorrect type"); 
  valid_types_.resize(4);
  valid_types_[0] = GSTEPPER_EXRK2;
  valid_types_[1] = GSTEPPER_EXRK4;
  valid_types_[2] = GSTEPPER_BDFAB;
  valid_types_[3] = GSTEPPER_BDFEXT;
 
  assert(valid_types_.contains(isteptype_) && "Invalid stepper type"); 

  gbc_ = new GBC(*grid_);
  utmp_.resize(tmp.size()); utmp_ = tmp;
  init(u, traits.iorder);
  
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GBurgers<TypePack>::~GBurgers()
{
  if ( gmass_   != NULLPTR ) delete gmass_;
  if ( gimass_  != NULLPTR ) delete gimass_;
//if ( gflux_   != NULLPTR ) delete gflux_;
  if ( ghelm_   != NULLPTR ) delete ghelm_;
  if ( gadvect_ != NULLPTR ) delete gadvect_;
  if ( gpdv_    != NULLPTR ) delete gpdv_;
  if ( gbc_     != NULLPTR ) delete gbc_;

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
template<typename TypePack>
void GBurgers<TypePack>::dt_impl(const Time &t, State &u, Time &dt)
{
   GSIZET ibeg, iend;
   GFTYPE dtmin, umax;
   GTVector<GFTYPE> *drmin  = &grid_->minnodedist();
   GElemList        *gelems = &grid_->elems();

   // This is an estimate. The minimum length on each element,
   // computed in GGrid object is divided by the maximum of
   // the state variable on each element:
   dt = 1.0;
   dtmin = std::numeric_limits<GFTYPE>::max();

   if ( bpureadv_ ) { // pure (linear) advection
     for ( auto k=1; k<u.size(); k++ ) { // each u
       for ( auto e=0; e<gelems->size(); e++ ) { // for each element
         ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
         u[k].range(ibeg, iend);
         umax = u[k].max();
         dtmin = MIN(dtmin, (*drmin)[e] / umax);
       }
       u[k].range_reset();
     }
   }
   else {             // nonlinear advection
     for ( auto k=0; k<u.size(); k++ ) { // each c
       for ( auto e=0; e<gelems->size(); e++ ) { // for each element
         ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
         u[k].range(ibeg, iend);
         umax = u[k].max();
         dtmin = MIN(dtmin, (*drmin)[e] / umax);
       }
       u[k].range_reset();
     }
   }

   GComm::Allreduce(&dtmin, &dt, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm_);

} // end of method dt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : dudt_impl
// DESC   : Compute RHS for explicit schemes
// ARGS   : t  : time
//          u  : state. If doing full nonlinear problem, the
//               u[0] = u[1] = u[2] are nonlinear fields.
//               If doing pure advection, only the first
//               element is the field being advected; the
//               remainder of the State elements are the 
//               (linear) velocity components that are not
//               updated.
//          dt : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::dudt_impl(const Time &t, State &u, Time &dt, Derivative &dudt)
{
  assert(!bconserved_ &&
         "conservation not yet supported"); 

  // NOTE:
  // Make sure that, in init(), Helmholtz op is using only
  // (mass isn't being used) weak Laplacian, or there will 
  // be problems. This is required for explicit schemes, for
  // which this method is called.

  // Do heat equation RHS:
  if ( doheat_ ) {
    for ( auto k=0; k<u.size(); k++ ) {
      ghelm_->opVec_prod(u[k],utmp_[0]); // apply diffusion
      gimass_->opVec_prod(utmp_[0],dudt[k]); // apply M^-1
    }
    return;
  }

  // Do linear/nonlinear advection + dissipation RHS:

  // If non-conservative, compute RHS from:
  //     du/dt = -u.Grad u + nu nabla^2 u 
  // for each u


  if ( bpureadv_ ) { // pure linear advection
    // Remember: only the first element of u is the variable;
    //           the remainder should be the adv. velocity components: 
    for ( auto j=0; j<u.size()-1; j++ ) c_[j] = u[j+1];
    gadvect_->apply(u[0], c_, utmp_, dudt[0]); // apply advection
    ghelm_->opVec_prod(u[0],utmp_[0]); // apply diffusion
    utmp_[0] -= dudt[0];
    gimass_->opVec_prod(utmp_[0],dudt[0]); // apply M^-1
  }
  else {             // nonlinear advection
    for ( auto k=0; k<u.size(); k++ ) {
      gadvect_->apply(u[k], u, utmp_, dudt[k]);
      ghelm_->opVec_prod(u[k],utmp_[0]); // apply diffusion
      utmp_[0] += dudt[k];
      gimass_->opVec_prod(utmp_[0],dudt[k]); // apply M^-1
    }
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl
// DESC   : Step implementation method  entry point
// ARGS   : t   : time
//          u   : state
//          dt  : time step
//          uout: updated state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::step_impl(const Time &t, State &uin, State &ub, Time &dt, Derivative &uout)
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
  
} // end of method step_impl


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
template<typename TypePack>
void GBurgers<TypePack>::step_multistep(const Time &t, State &uin, State &ub, Time &dt, State &uout)
{
  assert(FALSE && "Multistep methods not yet available");

  
} // end of method step_multistep


//**********************************************************************************
//**********************************************************************************
// METHOD : step_exrk
// DESC   : Take a step using Explicit RK method
// ARGS   : t   : time
//          uin : input state;; must not be modified
//          dt  : time step
//          uout: updated state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::step_exrk(const Time &t, State &uin, State &ub, Time &dt, State &uout)
{

  // If non-conservative, compute RHS from:
  //     du/dt = M^-1 ( -u.Grad u + nu nabla u ):
  // for each u

  // Set tmp arrays from member utmp_ data:
  GTVector<GTVector<GFTYPE>*> utmp1(utmp_.size()-itorder_-1);
  GTVector<GTVector<GFTYPE>*> utmp2(itorder_+1);
  for ( auto j=0; j<itorder_+1; j++ ) utmp2[j] = utmp_[j];
  for ( auto j=0; j<utmp_.size()-itorder_-1; j++ ) utmp1[j] = utmp_[itorder_+1+j];


  // Cycle over stages:
  for ( auto j=0; j<uin.size(); j++ ) *uout[j] = *uin[j];
  for ( auto k=0; k<itorder_; k++ ) {
    apply_bc_impl(t, uout, ub);
    gexrk_.step(t, uout[k], dt, utmp1, utmp2);
    for ( auto j=0; j<uin.size(); j++ ) *uout[k] = *utmp2[k];
  }

} // end of method step_exrk



//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_callback
// DESC   : Set the callback object for updating Dirichlet 
//          boundary state
// ARGS   : 
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::set_bdy_callback(std::function<void(Time &t, State &u,
                                          State &ub)> *callback)
{
  assert(gbc_ != NULLPTR && "Boundary operator not instantiated");

  bdy_update_callback_ = callback;
  gbc_->set_dirichlet_callback(bdy_update_callback_);
  
} // end of method set_bdy_callback


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
template<typename TypePack>
void GBurgers<TypePack>::init(State &u, GBurgers::Traits &traits)
{
  GString serr = "GBurgers<TypePack>::init: ";

  // Find multistep/multistage time stepping coefficients:
  GMultilevel_coeffs_base<GFTYPE> *tcoeff_obj; // time deriv coeffs
  GMultilevel_coeffs_base<GFTYPE> *acoeff_obj; // adv op. coeffs
  switch ( isteptype_ ) {
    case GSTEPPER_EXRK2:
    case GSTEPPER_EXRK4:
      gexrk_.setOrder(itorder_);
      gexrk_.setRHSfunction(this->dudt);
      break;
    case GSTEPPER_BDFAB:
      dthist_.resize(MAX(itorder_,inorder_));
      tcoeff_obj = new G_BDF<GFTYPE>(itorder_, dthist_);
      acoeff_obj = new G_AB<GFTYPE> (inorder_, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj; acoeffs_ = acoeff_obj;
      break;
    case GSTEPPER_BDFEXT:
      dthist_.resize(MAX(itorder_,inorder_));
      tcoeff_obj = new G_BDF<GFTYPE>(itorder_, dthist_);
      acoeff_obj = new G_EXT<GFTYPE>(inorder_, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj; acoeffs_ = acoeff_obj;
      break;
  }
  delete tcoeff_obj;
  delete acoeff_obj;
  
  // Instantiate spatial discretization operators:
  gmass_   = new GMass(*grid_);
  ghelm_   = new GHelmholtz(*grid_);
  
  if ( isteptype_ ==  GSTEPPER_EXRK2 
    || isteptype_ == GSTEPPER_EXRK4 ) {
    gimass_ = new GMass(*grid_, TRUE); // create inverse of mass
  }

  // If doing semi-implicit time stepping; handle viscous term 
  // (linear) inplicitly, which implies using full Helmholtz operator:
  if ( isteptype_ == GSTEPPER_BDFAB || isteptype_ == GSTEPPER_BDFEXT ) {
    ghelm_->set_mass(*gmass_);
  }

  if ( bconserved_ && !doheat_ ) {
    assert(FALSE && "Conservation not yet supported");
    gpdv_  = new GpdV(*grid_,*gmass_);
//  gflux_ = new GFlux(*grid_);
    assert( (gmass_   != NULLPTR
          && ghelm_   != NULLPTR
          && gpdv_    != NULLPTR) && "1 or more operators undefined");
  }
  if ( !bconserved_ && !doheat_ ) {
    gadvect_ = new GAdvect(*grid_);
    assert( (gmass_   != NULLPTR
          && ghelm_   != NULLPTR
          && gadvect_ != NULLPTR) && "1 or more operators undefined");
  }

  // If doing linear advection, set up a helper vector for
  // linear velocity:
  if ( bpureadv_ ) { 
    c_.resize(u.size()-1);
    for ( auto j=0; j<c_.size(); j++ ) c_[j] = u[j+1];
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
template<typename TypePack>
void GBurgers<TypePack>::cycle_keep(State &u)
{

  // Make sure following index map contains the 
  // correct time level information:
  //   u_keep[0] <--> time level n (most recent)
  //   u_keep[1] <--> time level n-1
  //   u_keep[2] <--> time level n-2 ...
  u_keep_ .resize(itorder_);
  for ( auto i=itorder_-1; i>=1; i-- ) u_keep_[i] = u_keep_[i+1];
  u_keep_[0] = u;

} // end of method cycle_keep


//**********************************************************************************
//**********************************************************************************
// METHOD : set_nu
// DESC   : Set viscosity quantity. This may be a field, or effectively a
//          scalar (where only element 0 contains valid data). If not set, 
//          Helmholtz op creates a default of nu = 1.
//          nu : viscosity vector (may be of length 1). Will be managed
//               by caller; only a pointer is used by internal operators.
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::set_nu(GTVector<GFTYPE> &nu)
{
  assert(ghelm_ != NULLPTR && "Init must be called first");
  nu_ = &nu; // Not sure this class actually needs this. May be removed later
  ghelm_->set_Lap_scalar(*nu_);

} // end of method set_nu


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_bc_impl
// DESC   : Apply global domain boundary conditions, ub
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::apply_bc_impl(const Time &t, State &u, State &ub)
{
  GTVector<GINT>     *igbdy     = &grid_->igbdy();

  // Use indirection to set the global field node values
  // with domain boundary data. ub must be updated outside 
  // of this method.

  // NOTE: This is useful to set Dirichlet-type bcs only. 
  // Neumann bcs type have to be set with the
  // differential operators themselves
 
  for ( GSIZET k=0; k<u.size(); k++ ) {
    for ( GSIZET j=0; j<igbdy->size(); j++ ) {
      u[k][(*igbdy)[j]] = ub[k][j];
    } 
  } 
  
} // end of method apply_bc_impl

