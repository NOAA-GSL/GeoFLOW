//==================================================================================
// Module       : gexrk_stepper.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object representing an Explicit RK stepper of a specified order
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gexrk_stepper.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order/ # stages
// ARGS   : iorder: truncation order
//**********************************************************************************
template<typename T>
GExRKStepper<T>::GExRKStepper(GSIZET iorder)
:
iorder_               (iorder),
rhs_callback_         (NULLPTR)
{
  butcher_.setOrder(iorder_);
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
//                 k_m = RHS( t^n + alpha_m * dt, u^n + h Sum_j=1^M-1 beta_mj k_j ),
//              and 
//                 RHS = RHS(t, u).
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              uin  : initial (entry) state, u^n
//              dt   : time step
//              tmp  : tmp space. Must have at least M+1 = iorder_+1 vectors.
//              uout : updated state, n^n+1
//               
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void GExRKStepper<T>::step(const T &t, GTVector<GTVector<T>*> &uin,
                        T &dt, GTVector<GTVector<T>*> &tmp, 
                        GTVector<GTVector<T>*> &uout)
{
  assert(rhs_callback_ != NULLPTR  && "RHS callback not set");

  GSIZET       i, j, n;
  GTVector<T> *u     = tmp[iorder_-1] ; 
  GTVector<T> *isum  = tmp[iorder_];    
  GTVector<T> *alpha = &butcher_.alpha();
  GTMatrix<T> *beta  = &butcher_.beta ();
  GTVector<T> *c     = &butcher_.c    ();
  
  
  for ( n=0; n<uin.size(); n++ ) { // for each state member
     *uout[n] = *uin[n];
    *isum = 0.0;
    (*rhs_callback_)( t+(*alpha)[i]*dt, *u, dt, *tmp[0]); // k_1 at stage 1
    // Compute k_m:
    // k_m = RHS( t^n + alpha_m * dt, u^n + h Sum_j=1^M-1 beta_mj k_j ),
    for ( i=1; i<iorder_-1; i++ ) { // cycle thru remaining stages minus 1
      for ( j=0,*isum=0.0; j<i; j++ ) *isum += (*tmp[j]) * ( (*beta)(i,j)*dt ); 
     *u    = uin[n] + (*isum);   
        rhs_callback_( t+(*alpha)[i]*dt, *u, dt, *tmp[i]); // k_i at stage i
       *uout[n] += (*tmp[i])*( (*c)[i]*dt ); // += dt * c_i * k_i
    }
    // Do contrib from final stage, M:
    for ( j=0,*isum=0.0; j<i; j++ ) *isum += (*tmp[j]) * ( (*beta)(iorder_-1,j)*dt ); 
    u = tmp[0]; *tmp[0] = uin[n] + (*isum);   
    rhs_callback_( t+(*alpha)[i]*dt, *u, dt, *tmp[1]); // k_i at stage M
   *uout[n] += (*tmp[1])*( (*c)[i]*dt ); // += dt * c_M * k_M
  }

} // end of method step


