//==================================================================================
// Module       : gab.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                and Adams-Bashforth scheme with variable timestep.
//                PDE. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gab.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order + dt history vector
// ARGS   : iorder: truncation order
//**********************************************************************************
template<typename T>
G_AB<T>::G_AB(GINT iorder, GTVector<T> &dthist)
: GMultilevel_coeffs_base<T>(iorder, dthist)
{

  assert(iorder_ >= 1 && 
         iorder_ <= 3 && "Invalid AB order");

  maxorder_ = 3;
  coeffs_.resize(iorder_);
  coeffs_ = 0.0;
  computeCoeffs();
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : 
//**********************************************************************************
template<typename T>
G_AB<T>::~G_AB()
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor
// DESC   : Default copy constructor
// ARGS   : a: G_AB object
//**********************************************************************************
template<typename T>
G_AB<T>::G_AB(const G_AB &a)
{
  iorder_   = a.iorder_;
  maxorder_ = a.maxorder_;
  coeffs_.resize(a.coeffs_.size()); coeffs_  = a.coeffs_;
  dthist_   = a.dthist_;
} // end of copy constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : computeCoeffs
// DESCRIPTION: Computes G_AB coefficients with variable timestep history.
//              NOTE: dthist_ pointer to timestep history buffer must be 
//                    set properly prior to entry, or 'exit' will be called.
//                    'exit' will be called if the timestep hist. buffer == NULL or
//                    if the number of buffer elements is too small as required
//                    by iorder_ variable.
// ARGUMENTS  : none.
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void G_AB<T>::computeCoeffs()
{

  assert(dthist_ != NULLPTR && dthist_->size() < iorder_  && "Invalid dt-history vector");

  T r1, r2, r3;
  if      ( iorder_ == 1 ) {
    coeffs_[0] = 1.0;
  }
  else if ( iorder_ == 2 ) {
    r1         = (*dthist_)[0] / (*dthist_)[1];
    coeffs_[0] = 1.0 + 0.5*r1;
    coeffs_[1] = 0.5*r1;
  }
  else if ( iorder_ == 3 ) {
    r1         = (*dthist_)[0] / (*dthist_)[1];
    r2         = (*dthist_)[1] / (*dthist_)[2];
    r3         = (*dthist_)[0] / (*dthist_)[2];

    coeffs_[1] = -r1* ( r3/3.0 + 0.5*r2 + 0.5);
    coeffs_[2] = -(0.5*r3 + coeffs_[1]*r2) / ( 1.0 + r2);
    coeffs_[0] = 1.0 - coeffs_[1] - coeffs_[2];
  }

} // end of method computeCoeffs


