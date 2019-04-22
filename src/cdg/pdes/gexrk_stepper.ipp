//==================================================================================
// Module       : gexrk_stepper.ipp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order/ # stages
// ARGS   : nstage: number stages not necessarily == truncation order
//**********************************************************************************
template<typename T>
GExRKStepper<T>::GExRKStepper(GSIZET nstage)
:
bRHS_                 (FALSE),
bapplybc_             (FALSE),
bupdatebc_            (FALSE),
nstage_               (nstage),
ggfx_                 (NULLPTR)
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
// METHOD     : step (1)
// DESCRIPTION: Computes one RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//
//              Given Butcher tableau, nodes, RK matrix, and weights, 
//              alpha_i, beta_ij, and c_i, , num stages, M,
//              the update is:
//                 u(t^n+1) = u(t^n) + h Sum_i=1^M c_i k_i,
//              where
//                 k_m = RHS( t^n + alpha_m * dt, u^n + h Sum_j=1^m-1 beta_mj k_j ),
//              and 
//                 RHS = RHS(t, u); and h = dt/M.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              ub   : bdy tendency
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void GExRKStepper<T>::step(const Time &t, const State &uin, State &uf, State &ub,  
                           const Time &dt, State &tmp, State &uout)
{
  assert(bRHS_  && "(1): RHS callback not set");

  GSIZET       i, j, m, n, nstate=uin.size();
  GFTYPE       h, tt;
  GTVector<T> *isum  ;
  GTVector<T> *alpha = &butcher_.alpha();
  GTMatrix<T> *beta  = &butcher_.beta ();
  GTVector<T> *c     = &butcher_.c    ();
  
  State u(nstate);   // tmp pointers of full state size

  resize(nstate);    // check if we need to resize K_
  
  h = dt ; 

  // Set temp space, initialize iteration:
  for ( n=0; n<nstate; n++ ) {
    u   [n] =  tmp[n];
   *uout[n] = *uin[n]; // deep copy 
  }
  isum = tmp[nstate];
  for ( j=0,n=0; j<MAX(nstage_-1,1); j++ ) {
    for ( i=0; i<nstate; i++ )  {
      K_[j][i] = tmp[nstate+1+n]; // set K storage from tmp space
      n++;
    }
  }

 
  for ( m=0; m<nstage_-1; m++ ) { // cycle thru stages minus 1
    // Compute k_m:
    //   k_m = RHS( t^n + alpha_m * h, u^n + h Sum_j=1^m-1 beta_mj k_j ),
    tt = t + (*alpha)[m]*h;

    for ( n=0; n<nstate; n++ ) { // for each state member, u
      for ( j=0,*isum=0.0; j<m; j++ )
        GMTK::saxpby(*isum,  1.0, *K_[j][n], (*beta)(m,j)*h);
     *u[n]  = (*uin[n]) + (*isum);
    }

    if ( bupdatebc_ ) bdy_update_callback_(tt, u, ub); 
    if ( bapplybc_  ) bdy_apply_callback_ (tt, u, ub); 
    rhs_callback_( tt, u, uf, ub, h, K_[m] ); // k_m at stage m

    // x^n+1 = x^n + h Sum_i=1^m c_i K_i, so
    // accumulate the sum in uout here: 
    for ( n=0; n<nstate; n++ ) { // for each state member, u
      if ( ggfx_ != NULLPTR ) {
        ggfx_->doOp(*K_[m][n], GGFX_OP_SMOOTH);
      }
      *isum     = (*K_[m][n])*( (*c)[m]*h );
      *uout[n] += *isum; // += h * c_m * k_m
    }
  } // end, m-loop

   // Do contrib from final stage, M = nstage_:
   tt = t + (*alpha)[nstage_-1]*h;
   for ( n=0; n<nstate; n++ ) { // for each state member, u
     for ( j=0,*isum=0.0; j<nstage_-1; j++ ) 
       GMTK::saxpby(*isum,  1.0, *K_[j][n], (*beta)(nstage_-1,j)*h);
     *u[n] = (*uin[n]) + (*isum);
   }

#if 0
   if ( ggfx_ != NULLPTR ) {
     for ( n=0; n<nstate; n++ ) { // for each state member, u
       ggfx_->doOp(*u[n], GGFX_OP_SMOOTH);
     }
   }
#endif

   if ( bupdatebc_ ) bdy_update_callback_(tt, u, ub); 
   if ( bapplybc_  ) bdy_apply_callback_ (tt, u, ub); 
   rhs_callback_( tt, u, uf, ub, h, K_[0]); // k_M at stage M

   if ( ggfx_ != NULLPTR ) {
     for ( n=0; n<nstate; n++ ) { // for each state member, uout
       ggfx_->doOp(*K_[0][n], GGFX_OP_SMOOTH);
     }
   }

   // Compute final output state, and set its bcs and
   // H1-smooth it:
   for ( n=0; n<nstate; n++ ) { // for each state member, u
    *isum     = (*K_[0][n])*( (*c)[nstage_-1]*h );
    *uout[n] += *isum; // += h * c_M * k_M
   }
   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout, ub); 

   if ( ggfx_ != NULLPTR ) {
     for ( n=0; n<nstate; n++ ) { // for each state member, uout
       ggfx_->doOp(*uout[n], GGFX_OP_SMOOTH);
     }
   }

   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout, ub); 
  
} // end of method step (1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : step (2)
// DESCRIPTION: Computes one RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//              The input state is overwritten.
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
//               
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void GExRKStepper<T>::step(const Time &t, State &uin, State &uf, State &ub,  
                           const Time &dt, State &tmp)
{
  assert(bRHS_  && "(2) RHS callback not set");

  GSIZET       i, m, j, n, nstate=uin.size();
  GFTYPE       h, tt;
  GTVector<T> *isum  ;
  GTVector<T> *alpha = &butcher_.alpha();
  GTMatrix<T> *beta  = &butcher_.beta ();
  GTVector<T> *c     = &butcher_.c    ();
  
  State u(nstate);        // tmp pointers of full state size
  State uout(nstate);     // tmp pointers of full output state size

  resize(nstate);         // check if we need to resize K_
  
  h = dt;

  // Set temp space: 
  //  size(tmp) = [nstate, nstate, 1, nstate*nstate]:
  for ( j=0; j<nstate; j++ ) {
    u   [j] = tmp[j];
    uout[j] = tmp[nstate+j];
   *uout[j] = *uin[j]; // deep copy
  }
  isum = tmp[2*nstate];
  for ( j=0,n=0; j<MAX(nstage_-1,1); j++ ) {
    for ( i=0; i<nstate; i++ )  {
      K_[j][i] = tmp[2*nstate+1+n]; // set K storage from tmp space
      n++;
    }
  }
  
  for ( m=0; m<nstage_-1; m++ ) { // cycle thru stages minus 1
    // Compute k_m:
    //   k_m = RHS( t^n + alpha_m * h, u^n + h Sum_j=1^m-1 beta_mj k_j ),
    tt = t+(*alpha)[m]*h;
    for ( n=0; n<nstate; n++ ) { // for each state member, u
      for ( j=0,*isum=0.0; j<m; j++ ) *isum += (*K_[j][n]) * ( (*beta)(m,j)*h );
     *u[n]  = (*uin[n]) + (*isum);
    }

    if ( bupdatebc_ ) bdy_update_callback_(tt, u, ub); 
    if ( bapplybc_  ) bdy_apply_callback_ (tt, u, ub); 
    if ( ggfx_ != NULLPTR ) {
      for ( n=0; n<nstate; n++ ) { // for each state member, u
        ggfx_->doOp(*u[n], GGFX_OP_SMOOTH);
      }
    }
    rhs_callback_( tt, u, uf, ub, h, K_[m] ); // k_m at stage m

    // x^n+1 = x^n + h Sum_i=1^m c_i K_i, so
    // accumulate the sum in uout here: 
    for ( n=0; n<nstate; n++ ) { // for each state member, u
      *uout[n] += (*K_[m][n])*( (*c)[m]*h ); // += h * c_m * k_m
    }
  }

   // Do contrib from final stage, M = nstage_:
   tt = t+(*alpha)[nstage_-1]*h;
   for ( n=0; n<nstate; n++ ) { // for each state member, u
     for ( j=0,*isum=0.0; j<nstage_-1; j++ ) *isum += (*K_[j][n]) * ( (*beta)(nstage_-1,j)*h );
     *u[n] = (*uin[n]) + (*isum);
      if ( ggfx_ != NULLPTR ) ggfx_->doOp(*u[n], GGFX_OP_SMOOTH);
   }
   if ( bupdatebc_ ) bdy_update_callback_(tt, u, ub); 
   if ( bapplybc_  ) bdy_apply_callback_ (tt, u, ub); 
   rhs_callback_( tt, u, uf, ub, h, K_[0]); // k_M at stage M

   for ( n=0; n<nstate; n++ ) { // for each state member, u
    *uout[n] += (*K_[0][n])*( (*c)[nstage_-1]*h ); // += h * c_M * k_M
   }

   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout, ub); 
   if ( ggfx_ != NULLPTR ) {
     for ( n=0; n<nstate; n++ ) { // for each state member, uouyt
       ggfx_->doOp(*uout[n], GGFX_OP_SMOOTH);
     }
   }
   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout, ub); 

  // deep copy tmp space to uin for return:
  for ( j=0; j<nstate; j++ ) {
   *uin[j] = *uout[j]; 
  }

} // end of method step (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : resize
// DESCRIPTION: Check if need to resize RK member data
// ARGUMENTS  : nstate : no. state vectors in state
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void GExRKStepper<T>::resize(GINT nstate)
{
  if ( K_.size() == 0 ) K_.resize(nstage_);

  for ( GSIZET j=0; j<K_.size(); j++ ) {
    if ( K_[j].size() < nstate ) K_[j].resize(nstate);
  }

} // end of method resize


