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
// Derived From : EquationBase.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gab.hpp"
#include "gext.hpp"
#include "gbdf.hpp"
#include "gexrk_stepper.hpp"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method  (1)
// DESC   : Instantiate with grid + state + tmp. Use for fully nonlinear
//          Burgers equation, heat equation.
// ARGS   : ggfx      : gather/scatter operator
//          grid      : grid object
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
GBurgers<TypePack>::GBurgers(GGFX<GFTYPE> &ggfx, Grid &grid, GBurgers<TypePack>::Traits &traits, GTVector<GTVector<GFTYPE>*> &tmp) :
doheat_          (traits.doheat),
bpureadv_      (traits.bpureadv),
bconserved_  (traits.bconserved),
bforced_        (traits.bforced),
bupdatebc_               (FALSE),
bsteptop_                (FALSE),
bvariabledt_ (traits.variabledt),
nsteps_                      (0),
itorder_        (traits.itorder),
inorder_        (traits.inorder),
courant_        (traits.courant),
nu_                    (NULLPTR),
gadvect_               (NULLPTR),
gmass_                  (NULLPTR),
gpdv_                  (NULLPTR),
//gflux_                 (NULLPTR),
gbc_                   (NULLPTR),
grid_                    (&grid),
ggfx_                    (&ggfx),
update_bdy_callback_   (NULLPTR),
steptop_callback_      (NULLPTR)
{
  static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
                "State is of incorrect type"); 
  assert(tmp.size() >= req_tmp_size() && "Insufficient tmp space provided");
  assert(!(doheat_ && bpureadv_) && "Invalid PDE configuration");

  // Check if specified stepper type is valid:
  GBOOL bfound;
  valid_types_.resize(4);
  valid_types_[0] = "GSTEPPER_EXRK";
  valid_types_[1] = "GSTEPPER_BDFAB";
  valid_types_[2] = "GSTEPPER_BDFEXT";
  bfound = valid_types_.contains(traits.steptype, isteptype_);
  assert( bfound && "Invalid stepping method specified");

  comm_ = ggfx_->getComm();

  gbc_ = new GBC(*grid_);
  utmp_.resize(tmp.size()); utmp_ = tmp;
  init(traits);
  
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
  if ( gexrk_   != NULLPTR ) delete gexrk_;

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
   GString serr = "GBurgers<TypePack>::dt_impl: ";
   GSIZET ibeg, iend;
   GFTYPE dtmin, dt1, umax;
   GFTYPE        drmin  = grid_->minnodedist();
   GElemList    *gelems = &grid_->elems();

   // This is an estimate. The minimum length on each element,
   // computed in GGrid object is divided by the maximum of
   // the state variable on each element:
   dtmin = std::numeric_limits<GFTYPE>::max();

   if ( bpureadv_ ) { // pure (linear) advection
     for ( auto k=1; k<u.size(); k++ ) { // each advecting u
       for ( auto e=0; e<gelems->size(); e++ ) { // for each element
         ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
         u[k]->range(ibeg, iend);
         umax = u[k]->amax();
         dtmin = MIN(dtmin, drmin / umax);
       }
       u[k]->range_reset();
     }
   }
   else {             // nonlinear advection
     for ( auto k=0; k<u.size(); k++ ) { // each advecting u
       for ( auto e=0; e<gelems->size(); e++ ) { // for each element
         ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
         u[k]->range(ibeg, iend);
         umax = u[k]->amax();
         dtmin = MIN(dtmin, drmin / umax);
       }
       u[k]->range_reset();
     }
   }

   GComm::Allreduce(&dtmin, &dt1, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm_);

   // Limit any timestep-to-timestep increae to 10%:
   dt = MIN(dt1*courant_, 1.1*dt);

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
//          uf : forcing tendency state vector
//          ub : bdy state vector
//          dt : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::dudt_impl(const Time &t, const State &u, const State &uf,  const State &ub, const Time &dt, Derivative &dudt)
{
  assert(!bconserved_ &&
         "conservation not yet supported"); 

  GString serr = "GBurgers<TypePack>::dudt_impl: ";

  // NOTE:
  // Make sure that, in init(), Helmholtz op is using only
  // weak Laplacian (q * mass isn't being used), or there will 
  // be problems. This is required for explicit schemes, for
  // which this method is called.
  // Do heat equation RHS:
  //     du/dt = nu nabla^2 u -->
  //     du/dt = -nu L u (in descrete form) 
  if ( doheat_ ) {
    for ( auto k=0; k<u.size(); k++ ) {
      ghelm_->opVec_prod(*u[k], uoptmp_, *urhstmp_[0]); // apply diffusion
     *urhstmp_[0] *= -1.0; // weak Lap op is neg on RHS
      gimass_->opVec_prod(*urhstmp_[0], uoptmp_, *dudt[k]); // apply M^-1
      if ( bforced_ && uf[k] != NULLPTR ) *dudt[k] += *uf[k];
    }
    return;
  }

  // Do linear/nonlinear advection + dissipation RHS:

  // If non-conservative, compute RHS from:
  //     du/dt = -c(u).Grad u + nu nabla^2 u 
  // for each u, where c may be indep of u (pure advection):

  if ( bpureadv_ ) { // pure linear advection
    // Remember: adv. velocity, c_ should already be set in
    //           main entry point, step_impl() method:
    GTimerStart("advection_time");
    gadvect_->apply   (*u[0], c_ , uoptmp_, *dudt[0]); // apply advection
    GTimerStop("advection_time");
    ghelm_->opVec_prod(*u[0], uoptmp_, *urhstmp_[0]);  // apply diffusion
    GMTK::saxpby<GFTYPE>(*urhstmp_[0], -1.0, *dudt[0], -1.0);
    gimass_->opVec_prod(*urhstmp_[0], uoptmp_, *dudt[0]); // apply M^-1
    if ( bforced_ ) *dudt[0] += *uf[0];
  }
  else {             // nonlinear advection
    for ( auto k=0; k<GDIM; k++ ) {
      gadvect_->apply(*u[k], u, uoptmp_, *dudt[k]);     // apply advection
      ghelm_->opVec_prod(*u[k], uoptmp_, *urhstmp_[0]); // apply diffusion
      GMTK::saxpby<GFTYPE>(*urhstmp_[0], -1.0, *dudt[k], -1.0);
      gimass_->opVec_prod(*urhstmp_[0], uoptmp_, *dudt[k]); // apply M^-1
      if ( bforced_ ) *dudt[k] += *uf[k];
    }
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (1)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//                If doing pure advection, this state contains the *unevolved*
//                advection velocities, ci, which must be separated from the 
//                single evolved state variable:
//                     uin = [u_evolve, c1, c2, c3]
//          uf  : force-tendency vector
//          ub  : bdy vector
//          dt  : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::step_impl(const Time &t, State &uin, State &uf, State &ub, const Time &dt)
{

  // If there's a top-of-the-timestep callback, 
  // call it here:
  if ( bsteptop_ ) {
    steptop_callback_(t, uin, dt);
  }

  // Set evolved state vars from input state.
  // These are not deep copies:
  if ( bpureadv_ || doheat_ ) {
    uevolve_[0] = uin[0]; 
    for ( auto j=0; j<c_.size(); j++ ) c_ [j] = uin[j+1];
  }
  else {
    for ( auto j=0; j<uin.size(); j++ ) uevolve_ [j] = uin[j];
  }

  switch ( isteptype_ ) {
    case GSTEPPER_EXRK:
      for ( auto j=0; j<uold_.size(); j++ ) *uold_[j] = *uevolve_[j];
      step_exrk(t, uold_, uf, ub, dt, uevolve_);
      break;
    case GSTEPPER_BDFAB:
    case GSTEPPER_BDFEXT:
      step_multistep(t, uevolve_, uf, ub, dt);
      break;
  }

} // end of method step_impl (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (2)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//                If doing pure advection, this state contains the *unevolved*
//                advection velocities, which must be separated from the 
//                single evolved state variable.
//          ub  : bdy vector
//          dt  : time step
//          uout: output state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::step_impl(const Time &t, const State &uin, State &uf,  State &ub, const Time &dt, State &uout)
{
  assert(FALSE && "step_impl(2) not available");

} // end of method step_impl (2)


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
//          uf  : force tendency vector
//          ub  : bdy vector
//          dt  : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::step_multistep(const Time &t, State &uin, State &uf, State &ub, const Time &dt)
{
  assert(FALSE && "Multistep methods not yet available");
  
} // end of method step_multistep


//**********************************************************************************
//**********************************************************************************
// METHOD : step_exrk
// DESC   : Take a step using Explicit RK method
// ARGS   : t   : time
//          uin : input state; must not be modified
//          uf  : force tendency vector
//          ub  : bdy vector
//          dt  : time step
//          uout: output/updated state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::step_exrk(const Time &t, State &uin, State &uf, State &ub, const Time &dt, State &uout)
{
  assert(gexrk_ != NULLPTR && "GExRK operator not instantiated");

  // If non-conservative, compute RHS from:
  //     du/dt = M^-1 ( -u.Grad u + nu nabla u ):
  // for each u

  // GExRK stepper steps entire state over one dt:
  gexrk_->step(t, uin, uf, ub, dt, urktmp_, uout);
  GMTK::constrain2sphere(*grid_, uout);

} // end of method step_exrk


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initialize equation object
// ARGS   : traits: GBurgers::Traits variable
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GBurgers<TypePack>::init(GBurgers::Traits &traits)
{
  GString serr = "GBurgers<TypePack>::init: ";

  GBOOL  bmultilevel = FALSE;
  GSIZET n, nsolve = bpureadv_ || doheat_ ? 1 : GDIM;
  GSIZET nc = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GINT   nop;

  // Find multistep/multistage time stepping coefficients:
  GMultilevel_coeffs_base<GFTYPE> *tcoeff_obj=NULLPTR; // time deriv coeffs
  GMultilevel_coeffs_base<GFTYPE> *acoeff_obj=NULLPTR; // adv op. coeffs

  // If doing pure advection, set advection 
  // components from input state vector, so
  // allocate size:
  uevolve_.resize(nsolve); // state var to evolve
  if ( bpureadv_ || doheat_ ) {
    c_.resize(nc);    // adevective vel components
  }

  std::function<void(const Time &t,                    // RHS callback function
                     const State  &uin,
                     const State  &uf,
                     const State  &ub,
                     const Time &dt,
                     State &dudt)> rhs
                  = [this](const Time &t,           
                     const State  &uin, 
                     const State  &uf, 
                     const State  &ub, 
                     const Time &dt,
                     State &dudt){dudt_impl(t, uin, uf, ub, dt, dudt);}; 

  std::function<void(const Time &t,                    // Bdy apply callback function
                     State  &uin, 
                     State &ub)> applybc 
                  = [this](const Time &t,              
                     State  &uin, 
                     State &ub){apply_bc_impl(t, uin, ub);}; 

  switch ( isteptype_ ) {
    case GSTEPPER_EXRK:
      gexrk_ = new GExRKStepper<GFTYPE>(itorder_);
      gexrk_->setRHSfunction(rhs);
      gexrk_->set_apply_bdy_callback(applybc);
      gexrk_->set_ggfx(ggfx_);
      // Set 'helper' tmp arrays from main one, utmp_, so that
      // we're sure there's no overlap:
      uold_   .resize(nsolve); // solution at time level n
      urktmp_ .resize(nsolve*(itorder_+1)+1); // RK stepping work space
      urhstmp_.resize(1); // work space for RHS
      nop = utmp_.size()-uold_.size()-urktmp_.size()-urhstmp_.size();
      assert(nop > 0 && "Invalid operation tmp array specification");
      uoptmp_ .resize(nop); // RHS operator work space
      // Make sure there is no overlap between tmp arrays:
      n = 0;
      for ( GSIZET j=0; j<nsolve         ; j++, n++ ) uold_   [j] = utmp_[n];
      for ( GSIZET j=0; j<urktmp_ .size(); j++, n++ ) urktmp_ [j] = utmp_[n];
      for ( GSIZET j=0; j<urhstmp_.size(); j++, n++ ) urhstmp_[j] = utmp_[n];
      for ( GSIZET j=0; j<uoptmp_ .size(); j++, n++ ) uoptmp_ [j] = utmp_[n];
      break;
    case GSTEPPER_BDFAB:
      dthist_.resize(MAX(itorder_,inorder_));
      tcoeff_obj = new G_BDF<GFTYPE>(itorder_, dthist_);
      acoeff_obj = new G_AB<GFTYPE> (inorder_, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj->getCoeffs(); 
      acoeffs_ = acoeff_obj->getCoeffs();
      uold_   .resize(nsolve); // solution at time level n
      for ( GSIZET j=0; j<nsolve; j++ ) uold_[j] = utmp_[j];
      bmultilevel = TRUE;
      break;
    case GSTEPPER_BDFEXT:
      dthist_.resize(MAX(itorder_,inorder_));
      tcoeff_obj = new G_BDF<GFTYPE>(itorder_, dthist_);
      acoeff_obj = new G_EXT<GFTYPE>(inorder_, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj->getCoeffs(); 
      acoeffs_ = acoeff_obj->getCoeffs();
      urhstmp_.resize(utmp_.size()-urktmp_.size());
      for ( GSIZET j=0; j<utmp_.size(); j++ ) urhstmp_[j] = utmp_[j];
      bmultilevel = TRUE;
      break;
    default:
      assert(FALSE && "Invalid stepper type");
  }
  if ( tcoeff_obj != NULLPTR ) delete tcoeff_obj;
  if ( acoeff_obj != NULLPTR ) delete acoeff_obj;
  
  // Instantiate spatial discretization operators:
  gmass_   = new GMass(*grid_);
  ghelm_   = new GHelmholtz(*grid_);
  
  if ( isteptype_ ==  GSTEPPER_EXRK ) {
    gimass_ = new GMass(*grid_, TRUE); // create inverse of mass
  }

  // If doing semi-implicit time stepping; handle viscous term 
  // (linear) inplicitly, which implies using full Helmholtz operator:
  if ( isteptype_ == GSTEPPER_BDFAB || isteptype_ == GSTEPPER_BDFEXT ) {
    assert(FALSE && "Implicit time stepping not yet supported");
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

  // If doing a multi-step method, instantiate (deep) space for 
  // required time levels for state:
  if ( bmultilevel ) {
    ukeep_ .resize(itorder_);
    for ( auto i=0; i<itorder_-1; i++ ) { // for each time level
      ukeep_[i].resize(nsolve);
      for ( auto j=0; j<nsolve; j++ ) ukeep_[i][j] = new GTVector<GFTYPE>(grid_->ndof());
    }
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
  //   ukeep[0] <--> time level n (most recent)
  //   ukeep[1] <--> time level n-1
  //   ukeep[2] <--> time level n-2 ...
  ukeep_ .resize(itorder_);
  for ( auto i=itorder_-1; i>=1; i-- ) ukeep_[i] = ukeep_[i+1];
  ukeep_[0] = u;

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
void GBurgers<TypePack>::apply_bc_impl(const Time &t, State &u, const State &ub)
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
      if ( ( ub[k] == NULLPTR ) continue;
      for ( GSIZET j=0; j<(*igbdy)[m].size(); j++ ) { // set Dirichlet-like value
        (*u[k])[(*igbdy)[m][j]] = (*ub[k])[j];
      } 
    } 
  } 
  
} // end of method apply_bc_impl

//**********************************************************************************
//**********************************************************************************
// METHOD : req_tmp_size
// DESC   : Find required tmp size on GLL grid
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
GINT GBurgers<TypePack>::req_tmp_size()
{
  GINT isize = 0;
  GINT nstate = 0;
 
  isize  = 2*GDIM + 3;
  nstate = GDIM;
  if ( doheat_ || bpureadv_ ) nstate = 1;

  if ( isteptype_ == GSTEPPER_EXRK ) {
    isize += nstate * itorder_; 
  }

  return isize;
  
} // end of method req_tmp_size

