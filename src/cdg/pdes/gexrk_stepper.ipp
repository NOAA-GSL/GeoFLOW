//==================================================================================
// Module       : gexrk_stepper.ipp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order/ # stages
// ARGS   : 
//          traits: this::Traits structure 
//          grid  : Grid object
//**********************************************************************************
template<typename Grid,typename T>
GExRKStepper<Grid,T>::GExRKStepper(Traits &traits, Grid &grid)
:
bRHS_                 (FALSE),
bapplybc_             (FALSE),
bSSP_           (traits.bSSP),
norder_       (traits.norder),
nstage_       (traits.nstage),
grid_                 (&grid),
ggfx_               (NULLPTR)
{
  if ( !bSSP_ ) {
    butcher_ .setOrder(norder_); // nstage_ = norder_
  }
  else {
    assert( (norder_==2 && nstage_==2)
        ||  (norder_==3 && nstage_==4)
        ||  (norder_==3 && nstage_==3) ); 
  }

} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order/ # stages
// ARGS   : 
//          traits: this::Traits structure 
//**********************************************************************************
template<typename Grid,typename T>
GExRKStepper<Grid,T>::GExRKStepper(Traits &traits)
:
bRHS_                 (FALSE),
bapplybc_             (FALSE),
bSSP_           (traits.bSSP),
norder_       (traits.norder),
nstage_       (traits.nstage),
grid_               (NULLPTR),
ggfx_               (NULLPTR)
{
  if ( !bSSP_ ) {
    butcher_ .setOrder(norder_); // nstage_ is ignored
  }
  else {
    assert( (norder_==2 && nstage_==2)
        ||  (norder_==3 && nstage_==4)
        ||  (norder_==3 && nstage_==3) ); 
  }

} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : 
//**********************************************************************************
template<typename Grid,typename T>
GExRKStepper<Grid,T>::~GExRKStepper()
{
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : step (1)
// DESCRIPTION: Computes one RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step(const Time &t, const State &uin, State &uf,
                           const Time &dt, State &tmp, State &uout)
{

  assert(bRHS_  && "(1) RHS callback not set");
  if ( bSSP_ ) {
    step_ssp(t, uin, uf, dt, tmp, uout);
  }
  else {
    step_b(t, uin, uf, dt, tmp, uout);
  }

} // end, method step(1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : step (2)
// DESCRIPTION: Computes one RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//              The input state is overwritten.
//
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step(const Time &t, State &uin, State &uf, 
                           const Time &dt, State &tmp)
{
  assert(bRHS_  && "(2) RHS callback not set");
  if ( bSSP_ ) {
    step_ssp(t, uin, uf, dt, tmp);
  }
  else {
    step_b(t, uin, uf, dt, tmp);
  }

} // end, method step(2)

#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD     : step_b (1)
// DESCRIPTION: Computes one (Butcher) RK step at specified timestep. Note: callback 
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
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_b(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &tmp, State &uout)
{

  GSIZET       i, j, m, n, nstate=uin.size();
  GFTYPE       h, tt;
  GTVector<T> *isum  ;
  GTVector<T> *alpha = &butcher_.alpha(); // time stage coeff
  GTMatrix<T> *beta  = &butcher_.beta (); // partial weights for k_i
  GTVector<T> *c     = &butcher_.c    (); // weights for final sum over k_i
  
  State u(nstate);   // tmp pointers of full state size

  resize(nstate);    // check if we need to resize K_
  
  h = dt ; 

  // Set temp space, initialize iteration:
  for ( n=0; n<nstate; n++ ) {
    u   [n] =  tmp[n];
   *uout[n] = *uin[n]; // deep copy 
  }
  isum = tmp[nstate];
  for ( j=0,n=0; j<MAX(nstage_,1); j++ ) {
    for ( i=0; i<nstate; i++ )  {
      K_[j][i] = tmp[nstate+1+n]; // set K storage from tmp space
      n++;
    }
  }

  for ( m=0; m<nstage_; m++ ) { // cycle thru stages minus 1
    // Compute k_m:
    //   k_m = RHS( t^n + alpha_m * h, u^n + h Sum_j=1^m-1 beta_mj k_j ),
    tt = t + (*alpha)[m]*h;

    // Accumulate sum for argument for K_m:
    for ( j=0,*isum=0.0; j<m; j++ )
        GMTK::saxpby<T>(*isum,  1.0, *K_[j][0], (*beta)(m,j)*h);
    *u[0]  = (*uin[0]); *u[0] += (*isum); // no copy const. called
    
    rhs_callback_( tt, u, uf, h, K_[m] ); // k_m at stage m

    // x^n+1 = x^n + h Sum_i=1^m c_i K_i, so
    // accumulate the sum in uout here: 
    *isum     = (*K_[m][0]); *isum *= ( (*c)[m]*h );
    *uout[0] += *isum; // += h * c_m * k_m
  } // end, m-loop over stages

  
} // end of method step_b (1)
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_b (1)
// DESCRIPTION: Computes one (Butcher) RK step at specified timestep. Note: callback 
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
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_b(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &tmp, State &uout)
{

  GSIZET       i, j, m, n, nstate=uin.size();
  GFTYPE       h, tt;
  GTVector<T> *isum  ;
  GTVector<T> *alpha = &butcher_.alpha(); // time stage coeff
  GTMatrix<T> *beta  = &butcher_.beta (); // partial weights for k_i
  GTVector<T> *c     = &butcher_.c    (); // weights for final sum over k_i
  
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

  tt = t ;
  if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
 
  for ( m=0; m<nstage_-1; m++ ) { // cycle thru stages minus 1
    // Compute k_m:
    //   k_m = RHS( t^n + alpha_m * h, u^n + h Sum_j=1^m-1 beta_mj k_j ),
    tt = t + (*alpha)[m]*h;

    // Accumulate sum for argument for K_m:
    for ( n=0; n<nstate; n++ ) { // for each state member, u
      for ( j=0,*isum=0.0; j<m; j++ )
        GMTK::saxpby<T>(*isum,  1.0, *K_[j][n], (*beta)(m,j)*h);
     *u[n]  = (*uin[n]); *u[n] += (*isum); // no copy const. called
    }
    if ( grid_!=NULLPTR ) GMTK::constrain2sphere<Grid,T>(*grid_, u);
    if ( bapplybc_ ) bdy_apply_callback_ (tt, u); 
    rhs_callback_( tt, u, uf, h, K_[m] ); // k_m at stage m

    // x^n+1 = x^n + h Sum_i=1^m c_i K_i, so
    // accumulate the sum in uout here: 
    for ( n=0; n<nstate; n++ ) { // for each state member, u
      if ( ggfx_ != NULLPTR ) {
        ggfx_->doOp(*K_[m][n], GGFX<GFTYPE>::Smooth());
      }
      *isum     = (*K_[m][n]); *isum *= ( (*c)[m]*h );
      *uout[n] += *isum; // += h * c_m * k_m
    }
    if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, uout);
    if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
  } // end, m-loop over stages


   // Do contrib from final stage, M = nstage_:
   tt = t + (*alpha)[nstage_-1]*h;
   for ( n=0; n<nstate; n++ ) { // for each state member, u
     for ( j=0,*isum=0.0; j<nstage_-1; j++ ) 
       GMTK::saxpby<T>(*isum,  1.0, *K_[j][n], (*beta)(nstage_-1,j)*h);
     *u[n]  = (*uin[n]); *u[n] += (*isum); // no copy const. called
   }
   if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, u);
   if ( bapplybc_  ) bdy_apply_callback_ (tt, u); 
   rhs_callback_( tt, u, uf, h, K_[0]); // k_M at stage M

   for ( n=0; ggfx_!=NULLPTR && n<nstate; n++ ) { // for each state member, uout
     ggfx_->doOp(*K_[0][n], GGFX<GFTYPE>::Smooth());
   }

   // Compute final output state, and set its bcs and
   // H1-smooth it:
   for ( n=0; n<nstate; n++ ) { // for each state member, u
    *isum     = (*K_[0][n]); *isum *= ( (*c)[nstage_-1]*h );
    *uout[n] += *isum; // += h * c_M * k_M
   }
   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
   if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, uout);

   for ( n=0; ggfx_!=NULLPTR && n<nstate; n++ ) { // for each state member, uout
     ggfx_->doOp(*uout[n], GGFX<GFTYPE>::Smooth());
   }

   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
  
} // end of method step_b (1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_b (2)
// DESCRIPTION: Computes one (Butcher) RK step at specified timestep. Note: callback 
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
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_b(const Time &t, State &uin, State &uf, 
                           const Time &dt, State &tmp)
{

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

    if ( bapplybc_  ) bdy_apply_callback_ (tt, u); 
    if ( ggfx_ != NULLPTR ) {
      for ( n=0; n<nstate; n++ ) { // for each state member, u
        ggfx_->doOp(*u[n], GGFX<GFTYPE>::Smooth());
      }
    }
    rhs_callback_( tt, u, uf, h, K_[m] ); // k_m at stage m

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
      if ( ggfx_ != NULLPTR ) ggfx_->doOp(*u[n], GGFX<GFTYPE>::Smooth());
   }
   if ( bapplybc_  ) bdy_apply_callback_ (tt, u); 
   rhs_callback_( tt, u, uf, h, K_[0]); // k_M at stage M

   for ( n=0; n<nstate; n++ ) { // for each state member, u
    *uout[n] += (*K_[0][n])*( (*c)[nstage_-1]*h ); // += h * c_M * k_M
   }

   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
   if ( ggfx_ != NULLPTR ) {
     for ( n=0; n<nstate; n++ ) { // for each state member, uouyt
       ggfx_->doOp(*uout[n], GGFX<GFTYPE>::Smooth());
     }
   }
   if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 

  // deep copy tmp space to uin for return:
  for ( j=0; j<nstate; j++ ) {
   *uin[j] = *uout[j]; 
  }

} // end of method step_b (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : resize
// DESCRIPTION: Check if need to resize RK member data
// ARGUMENTS  : nstate : no. state vectors in state
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::resize(GINT nstate)
{
  K_.resize(nstage_);

  for ( auto j=0; j<K_.size(); j++ ) {
    K_[j].resize(nstate);
  }

} // end of method resize



//**********************************************************************************
//**********************************************************************************
// METHOD     : step_ssp (1)
// DESCRIPTION: Computes one SSP RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//              The input state is not overwritten. This is a driver for specific
//              methods that handle different norder_ and nstage_.
//
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_ssp(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &tmp, State &uout)
{
  if      ( norder_ == 2 && nstage_ == 2 ) {
    step_ssp22(t, uin, uf, dt, tmp, uout);
  }
  else if ( norder_ == 3 && nstage_ == 3 ) {
    step_ssp33(t, uin, uf, dt, tmp, uout);
  }
  else if ( norder_ == 3 && nstage_ == 4 ) {
    step_ssp34(t, uin, uf, dt, tmp, uout);
  }
  else {
    assert(FALSE);
  }

} // end, method step_ssp (1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_s34
// DESCRIPTION: Computes one SSP RK34 (3rd order, 4 stages) step at 
//              specified timestep. Note: callback to RHS-computation 
//              function must be set prior to entry.
//
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_ssp34(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &tmp, State &uout)
{

  GSIZET       i, j, m, n, nstate=uin.size();
  GFTYPE       dtt, tt;
  
//State u(nstate);   // tmp pointers of full state size

  resize(nstate);    // check if we need to resize K_
  
  for ( j=0,n=0; j<MAX(nstage_-1,1); j++ ) {
    for ( i=0; i<nstate; i++ )  {
      K_[j][i] = tmp[n]; // set K storage from tmp space
      n++;
    }
  }


  // Stage 1:
  tt  = t;
  dtt = dt;
  step_euler(tt, uin, uf, dtt, K_[0]);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*K_[0][i],  0.5, *uin[i], 0.5);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, K_[0]); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, K_[0]);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*K_[0][i], GGFX<GFTYPE>::Smooth());
    }
  }
 
  // Stage 2:
  tt  = t + 0.5*dt;
  dtt = dt;
  step_euler(tt, K_[0], uf, dtt, K_[1]);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*K_[1][i],  0.5, *K_[0][i], 0.5);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, K_[1]); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, K_[1]);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*K_[1][i], GGFX<GFTYPE>::Smooth());
    }
  }
 
  // Stage 3:
  tt  = t + dt;
  dtt = dt;
  step_euler(tt, K_[1], uf, dtt, K_[2]);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*K_[2][i],  1.0/6.0, *K_[1][i], 1.0/6.0);
    GMTK::saxpby<T>(*K_[2][i],  1.0, *uin[i], 2.0/3.0);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, K_[2]); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, K_[2]);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*K_[2][i], GGFX<GFTYPE>::Smooth());
    }
  }

  // Stage 4:
  tt  = t + 0.5*dt;
  dtt = dt;
  step_euler(tt, K_[2], uf, dtt, uout);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*uout[i],  0.5, *K_[2][i], 0.5);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, uout);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*uout[i], GGFX<GFTYPE>::Smooth());
    }
  }
  
} // end of method step_ssp34


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_ssp33 (1)
// DESCRIPTION: Computes one SSP RK3e (3rd order, 3 stages) step at 
//              specified timestep. Note: callback to RHS-computation 
//              function must be set prior to entry.
//
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_ssp33(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &tmp, State &uout)
{

  GSIZET       i, j, m, n, nstate=uin.size();
  GFTYPE       dtt, tt;
  

  resize(nstate);    // check if we need to resize K_
  
  for ( j=0,n=0; j<MAX(nstage_-1,1); j++ ) {
    for ( i=0; i<nstate; i++ )  {
      K_[j][i] = tmp[n]; // set K storage from tmp space
      n++;
    }
  }


  // Stage 1:
  //   u(1)  = u^n + dt L(u^n);
  tt  = t;
  dtt = dt;
  step_euler(tt, uin, uf, dtt, K_[0]);   
  if ( bapplybc_  ) bdy_apply_callback_ (tt, K_[0]); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, K_[0]);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*K_[0][i], GGFX<GFTYPE>::Smooth());
    }
  }
 
  // Stage 2:
  tt  = t;
  dtt = dt;
  step_euler(tt, K_[0], uf, dtt, K_[1]);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*K_[1][i],  0.25, *uin[i], 0.75);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, K_[1]); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, K_[1]);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*K_[1][i], GGFX<GFTYPE>::Smooth());
    }
  }
 
  // Stage 3:
  tt  = t;
  dtt = dt;
  step_euler(tt, K_[1], uf, dtt, uout);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*uout[i],  2.0/3.0, *uin[i], 1.0/3.0);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, uout);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*uout[i], GGFX<GFTYPE>::Smooth());
    }
  }

} // end of method step_ssp33


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_ssp22 (1)
// DESCRIPTION: Computes one SSP RK2 (2nd order, 2 stages) step at 
//              specified timestep. Note: callback to RHS-computation 
//              function must be set prior to entry.
//
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_ssp22(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &tmp, State &uout)
{

  GSIZET       i, j, m, n, nstate=uin.size();
  GFTYPE       dtt, tt;
  

  resize(nstate);    // check if we need to resize K_
  
  for ( j=0,n=0; j<nstage_-1; j++ ) {
    for ( i=0; i<nstate; i++ )  {
      K_[j][i] = tmp[n]; // set K storage from tmp space
      n++;
    }
  }


  // Stage 1:
  //   u(1)  = u^n + dt L(u^n);
  tt  = t;
  dtt = dt;
  step_euler(tt, uin, uf, dtt, K_[0]);   
  if ( bapplybc_  ) bdy_apply_callback_ (tt, K_[0]); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_,K_[0]);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*K_[0][i], GGFX<GFTYPE>::Smooth());
    }
  }
 
  // Stage 2:
  tt  = t;
  dtt = dt;
  step_euler(tt, K_[0], uf, dtt, uout);   
  for ( i=0; i<nstate; i++ )  {
    GMTK::saxpby<T>(*uout[i],  0.5, *uin[i], 0.5);
  }
  if ( bapplybc_  ) bdy_apply_callback_ (tt, uout); 
  if ( grid_ != NULLPTR  ) GMTK::constrain2sphere<Grid,T>(*grid_, uout);
  for ( i=0; i<nstate; i++ )  {
    if ( ggfx_ != NULLPTR ) {
      ggfx_->doOp(*uout[i], GGFX<GFTYPE>::Smooth());
    }
  }
 
} // end of method step_ssp22


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_ssp (2)
// DESCRIPTION: Computes one SSP RK step at specified timestep. Note: callback 
//              to RHS-computation function must be set prior to entry.
//              The input state is overwritten. This is a driver for specific
//              methods that handle different norder_ and nstage_.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              dt   : time step
//              tmp  : tmp space. Must have at least NState*(M+1)+1 vectors,
//                     where NState is the number of state vectors.
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_ssp(const Time &t, State &uin, State &uf, 
                           const Time &dt, State &tmp)
{

  assert(FALSE && "Not implemented yet");

} // end of method step_ssp (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : step_euler 
// DESCRIPTION: Computes one Euler step of state:
//             
//               F = u^n + dt RHS(u^n,t^n)
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              uf   : forcing tendency
//              dt   : time step
//              uout : updated state, at t^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename Grid,typename T>
void GExRKStepper<Grid,T>::step_euler(const Time &t, const State &uin, State &uf, 
                           const Time &dt, State &uout)
{
  assert(bRHS_  && "RHS callback not set");

  rhs_callback_( t, uin, uf, dt, uout ); 
  for ( auto i=0; i<uout.size(); i++ ) {
    GMTK::saxpby<T>(*uout[i],  dt, *uin[i], 1.0);
  }

} // end, step_euler
