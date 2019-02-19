//==================================================================================
// Module       : gexrk_stepper.ipp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order/ # stages
// ARGS   : nstage: number stages not necessarily == truncation order
//**********************************************************************************
template<typename T>
GExRKStepper<T>::GExRKStepper(GSIZET nstage)
:
nstage_               (nstage),
rhs_callback_         (NULLPTR),
bdy_update_callback_  (NULLPTR),
bdy_apply_callback_   (NULLPTR)
{
  butcher_ .setOrder(nstage_);
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : 
//**********************************************************************************
template<typename T>
GExRKStepper<T>::~GExRKStepper()
{
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : step
// DESCRIPTION: Computes one RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//
//              Given Butcher tableau, alpha_i, beta_ij, and c_i, , num stages, M,
//              the update is:
//                 u(t^n+1) = u(t^n) + h Sum_i=1^M c_i k_i,
//              where
//                 k_m = RHS( t^n + alpha_m * dt, u^n + dt Sum_j=1^M-1 beta_mj k_j ),
//              and 
//                 RHS = RHS(t, u).
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, n^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void GExRKStepper<T>::step(const Time &t, const State &uin, State &ub,  
                           const Time &dt, State &tmp, State &uout)
{
  assert(rhs_callback_ != NULLPTR  && "RHS callback not set");

  GSIZET       i, j, n;
  GFTYPE       tt;
  GTVector<T> *isum  ;
  GTVector<T> *alpha = &butcher_.alpha();
  GTMatrix<T> *beta  = &butcher_.beta ();
  GTVector<T> *c     = &butcher_.c    ();
  
  GTVector<State> K(nstage_); // K for each stage (& all state members)
  State u(uin.size());   // tmp pointers of full state size

  // Set temp space:
  isum = tmp[uin.size()];
  for ( j=0,n=0; j<nstage_; j++ ) {
    for ( i=0; i<uin.size(); i++ )  {
      K[j][i] = tmp[uin.size()+1+n];
      n++;
    }
  }
  for ( j=0; j<uin.size(); j++ ) {
    u   [j] = tmp[j];
   *u   [j] = *uin[j]; // deep copy
   *uout[j] = *uin[j];
  }
  
  tt = t+(*alpha)[0]*dt;
  if ( bdy_update_callback_ != NULLPTR ) (*bdy_update_callback_)(tt, u, ub); 
  if ( bdy_apply_callback_  != NULLPTR ) (*bdy_apply_callback_ )(tt, u, ub); 
  (*rhs_callback_)( tt, u, dt, K[0]); // k_1 at stage 1

  for ( i=1; i<nstage_-1; i++ ) { // cycle thru remaining stages minus 1
    // Compute k_m:
    // k_m = RHS( t^n + alpha_m * dt, u^n + dt Sum_j=1^M-1 beta_mj k_j ),
    tt = t+(*alpha)[i]*dt;
    for ( n=0; n<uin.size(); n++ ) { // for each state member, u
      for ( j=0,*isum=0.0; j<i; j++ ) *isum += (*K[j][n]) * ( (*beta)(i,j)*dt );
     *u[n]  = (*uin[n]) + (*isum);
    }
    if ( bdy_update_callback_ != NULLPTR ) (*bdy_update_callback_)(tt, u, ub); 
    if ( bdy_apply_callback_  != NULLPTR ) (*bdy_apply_callback_ )(tt, u, ub); 
    (*rhs_callback_)( tt, u, dt, K[i]); // k_i at stage i
    *uout[n] += (*K[i][n])*( (*c)[i]*dt ); // += dt * c_i * k_i
   }

   // Do contrib from final stage, M:
   for ( n=0; n<uin.size(); n++ ) { // for each state member, u
     for ( j=0,*isum=0.0; j<nstage_-1; j++ ) *isum += (*K[j][n]) * ( (*beta)(nstage_-1,j)*dt );
     *u[n] = (*uin[n]) + (*isum);
   }
   tt = t+(*alpha)[nstage_-1]*dt;
   if ( bdy_update_callback_ != NULLPTR ) (*bdy_update_callback_)(tt, u, ub); 
   if ( bdy_apply_callback_  != NULLPTR ) (*bdy_apply_callback_ )(tt, u, ub); 
   (*rhs_callback_)( tt, u, dt, K[0]); // k_M at stage M

   for ( n=0; n<uin.size(); n++ ) { // for each state member, u
    *uout[n] += (*K[0][n])*( (*c)[i]*dt ); // += dt * c_M * k_M
   }

} // end of method step


