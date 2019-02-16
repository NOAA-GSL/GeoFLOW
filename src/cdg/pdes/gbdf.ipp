//==================================================================================
// Module       : gbdf.ipp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multilevel coefficients for 
//                a Backwards Differencing Formula (BDF) scheme with 
//                variable timestep.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with truncation order + dt history vector
// ARGS   : iorder: truncation order
//**********************************************************************************
template<typename T>
G_BDF<T>::G_BDF(GINT iorder, GTVector<T> &dthist)
: GMultilevel_coeffs_base<T>(iorder, dthist)
{
  assert(iorder_ >= 1 && iorder_ <= 3 && "Invalid AB order");

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
G_BDF<T>::~G_BDF()
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor
// DESC   : Default copy constructor
// ARGS   : a: G_BDF object
//**********************************************************************************
template<typename T>
G_BDF<T>::G_BDF(const G_BDF &a)
{
  iorder_   = a.iorder_;
  maxorder_ = a.maxorder_;
  coeffs_.resize(a.coeffs_.size()); coeffs_  = a.coeffs_;
  dthist_   = a.dthist_;
} // end of copy constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : computeCoeffs
// DESCRIPTION: Computes G_BDF coefficients with variable timestep history.
//              NOTE: dthist_ pointer to timestep history buffer must be 
//                    set properly prior to entry, or 'exit' will be called.
//                    'exit' will be called if the timestep hist. buffer == NULL or
//                    if the number of buffer elements is too small as required
//                    by iorder_ variable.
// ARGUMENTS  : none.
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void G_BDF<T>::computeCoeffs()
{

  assert(dthist_ != NULLPTR && dthist_->size() >= iorder_  && "Invalid dt-history vector");

  c_bdf_.resize(4,4);
  c_bdf_ = 0.0;

  
  c_bdf_(0,0) =  1.0;
  c_bdf_(1,0) =  2.0;
  c_bdf_(1,1) = -0.5;
  c_bdf_(2,0) =  3.0;
  c_bdf_(2,1) = -1.5;
  c_bdf_(2,2) =  1.0/3.0;
  c_bdf_(3,0) =  4.0;
  c_bdf_(3,1) = -3.0;
  c_bdf_(3,2) =  4.0/3.0;
  c_bdf_(3,3) = -0.25;

  for ( auto j=0; j<=iorder_; j++ ) coeffs_[j] = c_bdf_(iorder_-1,j);


} // end of method computeCoeffs


