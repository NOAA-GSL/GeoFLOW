//==================================================================================
// Module       : gbutcherrk.hpp
// Date         : 1/28/19 (DLR)
// Description  : Object computing multistage Butcher tableau for
//                explicit RK methods. 
//                  Butcher tableau takes the form:
//                         0 | 
//                    alpha2 | beta21
//                    alpha3 | beta31 beta32
//                    alpha4 | beta41 beta42 beta43
//                       ... | ...
//                    
//                    alpham | betam1 betam2 ... betamm-1
//                   --------------------------------------                
//                           | c1 c2 ...      c_m-1 c_m
//                   This class returns alpha, and c coeffs as vectors, and
//                   the beta as a matrix with lower diagonal filled with the beta
 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gbutcherrk.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with global order
// ARGS   : iorder: truncation order
//**********************************************************************************
template<typename T>
GButcherRK<T>::GButcherRK()
:
iorder_               (0)
{
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Instantiate with global order
// ARGS   : iorder: truncation order
//**********************************************************************************
template<typename T>
GButcherRK<T>::GButcherRK(GSIZET iorder)
:
iorder_               (iorder)
{
  computeCoeffs();
} // end of constructor (2) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : computeCoeffs
// DESCRIPTION: Computes Butcher tableau
// ARGUMENTS  : none.
// RETURNS    : none.
//**********************************************************************************
template<typename T>
void GButcherRK<T>::computeCoeffs()
{

  assert( (iorder_ == 2 || iorder_ == 4) && "Invalid RK order");

  alpha_.resize(iorder_)
  c_    .resize(iorder_)
  beta_ .resize(iorder_,iorder_)

  alpha_ = 0.0; beta_ = 0.0; c_ = 0.0;
  if ( iorder_ == 2 ) {
     alpha_[0] = 0;   c_[0] = 0.0;
     alpha_[1] = 0.5; c_[1] = 1.0;
     beta_(1,0) = 0.5;
    return;
  }
  
  if ( iorder_ == 4 ) {
     alpha_[0] = 0;   c_[0] = 1.0/6.0;
     alpha_[1] = 0.5; c_[1] = 1.0/3.0;
     alpha_[2] = 0.5; c_[2] = 1.0/3.0;
     alpha_[3] = 1.0; c_[2] = 1.0/6.0;
     beta_(1,0) = 0.5;
                       beta_(2,1) = 0.5;
                                          beta_(3,2) = 1.0;
    return;
  }
  
  


} // end of method computeCoeffs


