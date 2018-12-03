//==================================================================================
// Module       : gtpoint
// Date         : 7/1/18 (DLR)
// Description  : Encapsulates the access methods and data associated with
//                defining template 'point' object.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <limits>
#include "gtpoint.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD     : Constructor method, default (1)
// DESCRIPTION: 
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
GTPoint<T>::GTPoint() :
gdim_     (3),
eps_      (100*std::numeric_limits<T>::min()),
x1        (0),
x2        (0),
x3        (0),
x4        (0)
{
  lr_.resize(gdim_);
  nexp_.resize(gdim_);
  px_.resize(gdim_); 
  if ( gdim_ > 0 ) px_[0] = &x1;
  if ( gdim_ > 1 ) px_[1] = &x2;
  if ( gdim_ > 2 ) px_[2] = &x3;
  if ( gdim_ > 3 ) px_[3] = &x4;
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : Constructor method (2)
// DESCRIPTION: 
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
GTPoint<T>::GTPoint(GINT dim) :
gdim_     (dim),
eps_      (100*std::numeric_limits<T>::min()),
x1        (0),
x2        (0),
x3        (0),
x4        (0)
{
  lr_.resize(gdim_);
  nexp_.resize(gdim_);
  px_.resize(gdim_); 
  if ( gdim_ > 0 ) px_[0] = &x1;
  if ( gdim_ > 1 ) px_[1] = &x2;
  if ( gdim_ > 2 ) px_[2] = &x3;
  if ( gdim_ > 3 ) px_[3] = &x4;
} // end of constructor method (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : Constructor method (3)
// DESCRIPTION: 
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
GTPoint<T>::GTPoint(GINT dim, T e) :
gdim_     (dim),
eps_      (e),
x1        (0),
x2        (0),
x3        (0),
x4        (0)
{
  lr_.resize(gdim_);
  nexp_.resize(gdim_);
  px_.resize(gdim_); 
  if ( gdim_ > 0 ) px_[0] = &x1;
  if ( gdim_ > 1 ) px_[1] = &x2;
  if ( gdim_ > 2 ) px_[2] = &x3;
  if ( gdim_ > 3 ) px_[3] = &x4;
} // end of constructor method (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : Copy constructor
// DESCRIPTION: 
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
GTPoint<T>::GTPoint(const GTPoint<T> &e)
{
  gdim_ = e.gdim_;
  eps_  = e.eps_;
  x1 = e.x1;
  x2 = e.x2;
  x3 = e.x3;
  x4 = e.x4;
} // end, copy constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : Destructor
// DESCRIPTION: 
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
GTPoint<T>::~GTPoint()
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD     : setBracket
// DESCRIPTION: set bracket for float-type points for determining bounds (adds
//              'fuzziness')
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
void GTPoint<T>::setBracket(T e)
{
  eps_ = e;
} // end of method setBracket


//**********************************************************************************
//**********************************************************************************
// METHOD     : getBracket
// DESCRIPTION: gets bracket bounds
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
T GTPoint<T>::getBracket()
{
  return eps_; 
} // end of method getBracket

//**********************************************************************************
//**********************************************************************************
// METHOD     : resize
// DESCRIPTION: resizes dimension
// ARGUMENTS  :
// RETURNS    :
//**********************************************************************************
template<typename T>
void GTPoint<T>::resize(GINT dim)
{
  gdim_ = dim;
  lr_.resize(gdim_);
  nexp_.resize(gdim_);
  px_.resize(gdim_); 
  if ( gdim_ > 0 ) px_[0] = &x1;
  if ( gdim_ > 1 ) px_[1] = &x2;
  if ( gdim_ > 2 ) px_[2] = &x3;
  if ( gdim_ > 3 ) px_[3] = &x4;
} // end of method resize

#if 0
template class GTPoint<GDOUBLE>;
template class GTPoint <GFLOAT>;
template class GTPoint   <GINT>;
template class GTPoint  <GLONG>;
template class GTPoint <GLLONG>;
#endif
