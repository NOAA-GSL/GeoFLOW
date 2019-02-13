//==================================================================================
// Module       : gext.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                an extrapolation scheme with variable timestep.
//                PDE. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gext.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order + dt history vector
// ARGS   : iorder: truncation order
//**********************************************************************************
template<typename T>
G_EXT<T>::G_EXT(GSIZET iorder, GTVector<T> &dthist)
:
iorder_               (iorder),
maxorder_             (3),
dthist_               (&dthist)
{
  assert(iorder_ >= 1 && iorder_ <= 3 && "Invalid AB order");
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
G_EXT<T>::~G_EXT()
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor
// DESC   : Default copy constructor
// ARGS   : a: G_EXT object
//**********************************************************************************
template<typename T>
G_EXT<T>::G_EXT(const G_EXT &a)
{
  iorder_   = a.iorder_;
  maxorder_ = a.maxorder_;
  coeffs_.resize(a.coeffs_.size()); coeffs_  = a.coeffs_;
  dthist_   = a.dthist_;
} // end of copy constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : computeCoeffs
// DESCRIPTION: Computes G_EXT coefficients with variable timestep history.
//              NOTE: dthist_ pointer to timestep history buffer must be 
//                    set properly prior to entry, or 'exit' will be called.
//                    'exit' will be called if the timestep hist. buffer == NULL or
//                    if the number of buffer elements is too small as required
//                    by iorder_ variable.
// ARGUMENTS  : none.
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void G_EXT<T>::computeCoeffs()
{

  assert(dthist_ != NULLPTR && dthist_->size() < iorder_  && "Invalid dt-history vector");

  T r1, r2, r3, r4, r5;
  if      ( iorder_ == 1 ) {
    coeffs_[0] = 1.0;
  }
  else if ( iorder_ == 2 ) {
    r1         = (*dthist_)[0] / (*dthist_)[1];
    coeffs_[0] = 1.0 + r1;
    coeffs_[1] = -r1;
  }
  else if ( iorder_ == 3 ) {
    r1         = (*dthist_)[0] / (*dthist_)[1];
    r2         = (*dthist_)[1] / (*dthist_)[2];
    r3         = (*dthist_)[0] / (*dthist_)[2];
    r4         = (r2  + r3) / (1.0 + r2);
    r5         = (1.0 + r1) * (1.0 + r2);

    coeffs_[0] = 1.0 - r1*r1 - r3*r4 + r1*r5;;
    coeffs_[1] = r1*r1 - r1*r5;
    coeffs_[2] = r3* r4;
  }

} // end of method computeCoeffs

