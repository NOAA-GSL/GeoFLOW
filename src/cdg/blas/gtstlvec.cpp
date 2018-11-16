//==================================================================================
// Module      : gtvector
// Date        : 1/1/2018 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types. Derives from std::vector.
// Copyright   : Copyright 2018. Colorado State University. All rights reserved
// Derived from: std::vector
//==================================================================================
#include <cstddef>
#include <cstring>
#include "gtstlvec.hpp"
#include "cff_blas.h"

#if !defined(_G_VEC_CACHE_SIZE)
  # define _G_VEC_CACHE_SIZE 16
#endif

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Basic
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<class T> GTSTLVec<T>::GTSTLVec():
std::vector<T>(1),
data_   (NULLPTR),
n_      (1),
icsz_   (_G_VEC_CACHE_SIZE)
{

  data_ = std::vector<T>::data();
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with array size
// ARGS   : GSIZET n: array size
// RETURNS: none
//**********************************************************************************
template<class T> GTSTLVec<T>::GTSTLVec(GSIZET n):
std::vector<T>(n),
data_   (NULLPTR),
n_      (n),
icsz_   (_G_VEC_CACHE_SIZE)
{

  data_ = std::vector<T>::data();
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}

#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with array size, offsets
// ARGS   : 
// RETURNS: none
//**********************************************************************************
template<class T> GTSTLVec<T>::GTSTLVec(GSIZET n, GSIZET ib, GSIZET ie):
std::vector<T>(n),
data_   (NULLPTR),
n_      (n),
icsz_   (_G_VEC_CACHE_SIZE)
{
  data_ = data();
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: GTSTLVec argument
// DESC   : Instantitate with GTSTLVec argument
// ARGS   : GTSTLVec<T> argument
// RETURNS: 
//**********************************************************************************
template<class T> GTSTLVec<T>::GTSTLVec(GTSTLVec<T> &obj):
std::vector<T>(obj.size()),
data_   (NULLPTR),
n_      (obj.size()),
icsz_   (_G_VEC_CACHE_SIZE)
{
  data_ = std::vector<T>::data();
  std::memcpy(data_, obj.data(), sizeof(T)*n_);
  
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
}


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: raw pointer and size args
// DESC   : Instantitate with data block and size
// ARGS   : T *indata: pointer to external data block
//          GSIZET n : size of external data block
// RETURNS: 
//**********************************************************************************
template<class T> GTSTLVec<T>::GTSTLVec(T *indata, GSIZET n):
std::vector<T>(n),
data_   (NULLPTR),
n_      (n),
icsz_   (_G_VEC_CACHE_SIZE)
{
  data_ = std::vector<T>::data();
  std::memcpy(data_, indata, sizeof(T)*n_);
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
  updatedev();
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method
// DESC   : Override degault copy constructor
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<class T> GTSTLVec<T>::GTSTLVec(const GTSTLVec<T> &obj):
std::vector<T>(obj.size()),
data_   (NULLPTR),
n_      (obj.size()),
icsz_   (_G_VEC_CACHE_SIZE)
{
  data_ = std::vector<T>::data();
  std::memcpy(data_, obj.data(), sizeof(T)*n_);
  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   :
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T> GTSTLVec<T>::~GTSTLVec()
{
  #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
}


//**********************************************************************************
//**********************************************************************************
// METHOD : resize 
// DESC   : Resize vector/datablock
// ARGS   : GSIZET n: new vector size
// RETURNS: none.
//**********************************************************************************
template<class T> void GTSTLVec<T>::resize(GSIZET nnew)
{

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc exit data delete( data_[0:n_-1] )
  #endif

  std::vector<T>::resize(nnew);
  data_ = std::vector<T>::data();
  n_ = nnew;

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data create( data_[0:n_-1] )
  #endif
 
} // end, method resize


//**********************************************************************************
//**********************************************************************************
// METHOD : assignment operator= GTSTLVec
// DESC   : Equate to another GTSTLVec
// ARGS   : GTSTLVec<T> & right-hand side arg 
// RETURNS: GTSTLVec & 
//**********************************************************************************
template<class T> GTSTLVec<T> &GTSTLVec<T>::operator=(const GTSTLVec<T> &obj)
{
  if ( this == &obj)
  {
    return *this;
  }
  

  #if !defined(_G_USE_OPENACC)
  if ( data_ != NULLPTR &&  n_ != obj.size() ) 
  {
    delete [] data_;
    n_ = obj.size();
    data_ = new T[n_];
  }
  std::memcpy(data_, obj.data(), sizeof(T)*n_);
  #else
    GSIZET j;
    const T *dobj = obj.data();
    for ( j=0; j<n_; j++ )
    {
      data_[j] = dobj[j];
    }
  #endif

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
  
  return *this;
} // end, operator=(GTSTLVec &)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator = T
// DESC   : Equate to constant
// ARGS   : T arg
// RETURNS: GTSTLVec & 
//**********************************************************************************
template<class T> GTSTLVec<T> &GTSTLVec<T>::operator=(T a)
{
  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] = a;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

  return *this;
} // end, operator=


//**********************************************************************************
//**********************************************************************************
// METHOD : pointProd
// DESC   : point-by-point multiplication, returned in specified GTSTLVec
// ARGS   : ret: GTSTLVec & ret
//          obj: const GTSTLVec<>  factor
// RETURNS: GTSTLVec & 
//**********************************************************************************
template<class T> void GTSTLVec<T>::pointProd(GTSTLVec<T> &ret, const GTSTLVec<T> &obj) 
{

  #if !defined(_G_USE_OPENACC)
  if ( obj.size() < n_ || ret.size() < n_ )
  {
    std::cout << "GTSTLVec<T>::pointProd: " << "incompatible size" << std::endl;
    exit(1);
  }
  #endif

  GSIZET j; 
  const T *dobj=obj.data();
  T *dret=ret.data();
  for ( j=0; j<n_; j++ ) 
  {
    dret[j] = data_[j] * dobj[j];
  }

} // end, pointProd


//**********************************************************************************
//**********************************************************************************
// METHOD : constProd
// DESC   : multiply this by constant, return
// ARGS   : ret: GTSTLVec & ret vector
//          b  : T-type constant 
// RETURNS: GTSTLVec & 
//**********************************************************************************
template<class T> void GTSTLVec<T>::constProd(GTSTLVec<T> &ret, const T b) 
{
  #if !defined(_G_USE_OPENACC)
  if ( ret.size() < n_ )
  {
    std::cout << "GTSTLVec<T>::constProd: " << "incompatible size" << std::endl;
    exit(1);
  }
  #endif

  GSIZET j;
  T *dret=ret.data();

  for ( j=0; j<n_; j++ ) 
  {
    dret[j] = data_[j] * b;
  }

} // end, pointProd



//**********************************************************************************
//**********************************************************************************
// METHOD : add
// DESC   : point-by-point addition, returned in specified GTSTLVec:
//            ret = a*this + b*obj
// ARGS   : ret: GTSTLVec<T> return vector
//          obj: const GTSTLVec<T> summand
// RETURNS: GTSTLVec & 
//**********************************************************************************
template<class T> void GTSTLVec<T>::add(GTSTLVec<T> &ret, const GTSTLVec<T> &obj, const T a, const T b) 
{

  #if !defined(_G_USE_OPENACC)
  if ( obj.size() < n_  || ret.size() < n_ )
  {
    std::cout << "GTSTLVec<T>::add: " << "incompatible size" << std::endl;
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
  GSIZET j;
  const T *dobj=obj.data();
  T *dret=ret.data();
  for ( j=0; j<n_; j++ ) 
  {
    dret[j] = a*data_[j] + b*dobj[j];
  }
  #else
  szaxpby(dret, data_, &a, dobj, &b, &n_, icsz_);
  #endif

} // end, add 


//**********************************************************************************
//**********************************************************************************
// METHOD : operator += (1)
// DESC   :
// ARGS   : GTSTLVec &
// RETURNS: void
//**********************************************************************************
template<class T> void  GTSTLVec<T>::operator+=(GTSTLVec<T> &obj)
{

  T *dobj = obj.data();
  #if !defined(_G_USE_GBLAS)

  for ( GSIZET j=0; j<n_; j++ ) 
  {
    data_[j] += dobj[j];
  }
  #else
    T a=1.0;
    T b=1.0;
    dxaxpby(data_, &a, dobj, &b, &n_, &icsz_);
  #endif

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator+= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator += (2)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T> void GTSTLVec<T>::operator+=(T b)
{

  GSIZET j;
  for ( j=0; j<n_; j++ ) 
  {
    data_[j] += b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator+= (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -=  (1)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T> void GTSTLVec<T>::operator-=(GTSTLVec<T> &b)
{

  GSIZET j;
  T *p = b.data();
  for ( j=0; j<n_; j++ ) 
  {
    data_[j] -= p[j];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator-= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -= (2)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T> void GTSTLVec<T>::operator-=(T b)
{

  GSIZET j;
  for ( j=0; j<n_; j++ ) 
  {
    data_[j] -= b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operatori= (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (1)
// DESC   : point-by-point product
// ARGS   : GTSTLVec &
// RETURNS: void
//**********************************************************************************
template<class T> void GTSTLVec<T>::operator*=(GTSTLVec<T> &b)
{

  GSIZET j;
  T *p = b.data();
  for ( j=0; j<n_; j++ ) 
  {
    data_[j] *= p[j];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator*= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (2)
// DESC   : product of vector and constant
// ARGS   : GTSTLVec &
// RETURNS: void
//**********************************************************************************
template<class T> void GTSTLVec<T>::operator*=(T b)
{

  GSIZET j;
  for ( j=0; j<n_; j++ ) 
  {
    data_[j] *= b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator*= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : sub 
// DESC   : point-by-point subtraction, return in specified Vector:
//            ret = a*this + b*obj
// ARGS   : ret: GTSTLVec<T> return vector
//          obj: const GTSTLVec<T> summand
// RETURNS: GTSTLVec & 
//**********************************************************************************
template<class T> void GTSTLVec<T>::sub(GTSTLVec<T> &ret, const GTSTLVec<T> &obj, const T a, const T b) 
{

  #if !defined(_G_USE_OPENACC)
  if ( obj.size() < n_ || ret.size() < n_ )
  {
    std::cout << "GTSTLVec<T>::sub: " << "incompatible size" << std::endl;
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
  GSIZET j;
  const T *dobj=obj.data();
  T *dret=ret.data();
  for ( j=0; j<n_; j++ ) 
  {
    dret[j] = a*data_[j] - b*dobj[j];
  }
  #else
  GDOUBLE c = -b;
  szaxpby(dret, data_, &a, dobj, &c, &n_, icsz_);
  #endif

} // end, sub


//**********************************************************************************
//**********************************************************************************
// METHOD : SetCacheBlockingFactor
// DESC   : Set blocking factor
// ARGS   : GINT number of T type members in cache
// RETURNS: none.
//**********************************************************************************
template<class T> void GTSTLVec<T>::SetCacheBlockingFactor(GINT icsz)
{ 
  icsz_ = MAX(icsz,0);

} // end of method SetCacheBlocking

//************************************************************************************
//************************************************************************************
// METHOD : Updatehost
// DESC   : update data from device to host
// ARGS   : none.
// RETURNS: none.

template<class T> void GTSTLVec<T>::updatehost()
{
  #pragma acc update self( data_[0:n_-1] )
} // end of method updatehost


//************************************************************************************
//************************************************************************************
// METHOD : updatedev
// DESC   : Update data from host to device
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<class T> void GTSTLVec<T>::updatedev()
{
  #pragma acc update device(data_[0:n_-1])
} // end of method updatedev


//template class GTSTLVec <GFLOAT>;
template class GTSTLVec<GDOUBLE>;
//template class GTSTLVec   <GINT>;
//template class GTSTLVec <GSIZET>;

