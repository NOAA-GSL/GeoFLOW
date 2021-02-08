//==================================================================================
// Module       : gtpoint
// Date         : 7/1/18 (DLR)
// Description  : Encapsulates the access methods and data associated with
//                defining template 'point' object.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

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
x1        (0),
x2        (0),
x3        (0),
x4        (0),
eps_      (0)
{
  GEOFLOW_TRACE();
  if ( std::is_floating_point<T>::value ) {
    eps_ = 100*std::numeric_limits<T>::min();
  }
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
x1        (0),
x2        (0),
x3        (0),
x4        (0),
eps_      (0)
{
  GEOFLOW_TRACE();
  if ( std::is_floating_point<T>::value ) {
    eps_ = 100*std::numeric_limits<T>::min();
  }
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
  GEOFLOW_TRACE();
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
  GEOFLOW_TRACE();
  gdim_ = e.gdim_;
  eps_  = e.eps_;
  x1 = e.x1;
  x2 = e.x2;
  x3 = e.x3;
  x4 = e.x4;
  px_.resize(gdim_);
  if ( gdim_ > 0 ) px_[0] = &x1;
  if ( gdim_ > 1 ) px_[1] = &x2;
  if ( gdim_ > 2 ) px_[2] = &x3;
  if ( gdim_ > 3 ) px_[3] = &x4;
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
  GEOFLOW_TRACE();
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
  GEOFLOW_TRACE();
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
  GEOFLOW_TRACE();
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
  GEOFLOW_TRACE();
  gdim_ = dim;
  px_.resize(gdim_); 
  if ( gdim_ > 0 ) px_[0] = &x1;
  if ( gdim_ > 1 ) px_[1] = &x2;
  if ( gdim_ > 2 ) px_[2] = &x3;
  if ( gdim_ > 3 ) px_[3] = &x4;
} // end of method resize

