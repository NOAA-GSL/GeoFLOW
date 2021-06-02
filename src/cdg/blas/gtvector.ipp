//==================================================================================
// Module      : gtvector.ipp
// Date        : 1/1/2018 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types.
//                NOTE: Partial support for generalized 'index' that
//                      yields, start, stop, stride, pad indices.
//                The implementation file, gtvector.ipp, is not included
//                in this file, just the declarations.
//                NOTE: Math operations that accept a GTVector object
//                      argument may be treated as constant, by having
//                      size of input argument = 1, and taking as the
//                      operand only the first element.
// Copyright   : Copyright 2018. Colorado State University. All rights reserved
// Derived from: none.
//==================================================================================

#if !defined(_G_VEC_CACHE_SIZE)
  # define _G_VEC_CACHE_SIZE 16
#endif

#if defined(GEOFLOW_USE_CUBLAS)
#include "cublas.h"
#endif

#include "tbox/assert.hpp"
#include "tbox/tracer.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Basic
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<class T>
GTVector<T>::GTVector():
data_           (NULLPTR),
n_                    (0),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_        (TRUE)
{
  GEOFLOW_TRACE();
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;
}


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with array size
// ARGS   : GSIZET n: array size
// RETURNS: none
//**********************************************************************************
template<class T>
GTVector<T>::GTVector(GSIZET n):
data_           (NULLPTR),
n_                    (n),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_        (TRUE)
{
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
  data_ = new T [n_];
#endif

  assert(this->data_!= NULLPTR );
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with general index object
// ARGS   : GIndex gi: contains ibeg, iend, stride, etc.
// RETURNS: none
//**********************************************************************************
template<class T>
GTVector<T>::GTVector(GIndex &gi):
data_           (NULLPTR),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_        (TRUE)
{
  GEOFLOW_TRACE();
  gindex_ = gi;
  gindex_keep_ = gindex_;
  n_=gindex_.end()+1+gindex_.pad();

#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
  data_ = new T [n_];
#endif
  assert(this->data_!= NULLPTR );
}


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: GTVector argument
// DESC   : Instantitate with GTVector argument
// ARGS   : GTVector<T> argument
// RETURNS: 
//**********************************************************************************
template<class T>
GTVector<T>::GTVector(GTVector<T> &obj):
data_           (NULLPTR),
n_           (obj.size()),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_        (TRUE)
{
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
  data_ = new T [n_];
#endif
  assert(this->data_!= NULLPTR );
  
  for ( auto j=0; j<n_; j++ ) {
    this->data_[j] = obj[j];
  }
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;
  updatedev();
}


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: raw pointer and size args
// DESC   : Instantitate with data block and size
// ARGS   : indata: pointer to external data block
//          n      : size of external data block
//          istride: constant stride through external data
// RETURNS: 
//**********************************************************************************
template<class T>
GTVector<T>::GTVector(T *indata, GSIZET n, GSIZET istride):
data_           (NULLPTR),
n_            (n/istride),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_        (TRUE)
{
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
  data_ = new T [n_];
#endif
  assert(this->data_!= NULLPTR );

  GLLONG k=0;
  for ( auto j=0; j<n_; j++ ) {
    data_[j] = indata[k];
    k += istride;
  }
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;
  updatedev();
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method: raw pointer and size args + manage flag
// DESC   : Instantitate with data block and size
// ARGS   : indata: pointer to external data block
//          n      : size of external data block
//          istride: constant stride through external data
//          blocmgd: tell this class that the data should be managed locally 
//                   (is it owned by this class (TRUE) or by caller (FALSE))
// RETURNS: 
//**********************************************************************************
template<class T>
GTVector<T>::GTVector(T *indata, GSIZET n, GSIZET istride, GBOOL blocmgd):
data_           (NULLPTR),
n_            (n/istride),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_     (blocmgd)
{
  GEOFLOW_TRACE();
  if ( bdatalocal_ ) {
#if defined(GEOFLOW_USE_CUBLAS)
    if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
      cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
    }
    else {
      data_ = new T [n_];
    }
#else
    data_ = new T [n_];
#endif
    assert(this->data_!= NULLPTR );
    GLLONG k=0;
    for ( auto j=0; j<n_; j++ ) {
      data_[j] = indata[k];
      k += istride;
    }
    gindex_(n_, n_, 0, n_-1, 1,  0);
  }
  else {
    data_ = indata;
    gindex_(n_, n_, 0, n_-1, istride,  0);
  }
  gindex_keep_ = gindex_;
  updatedev();
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method
// DESC   : Override default copy constructor
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<class T>
GTVector<T>::GTVector(const GTVector<T> &obj):
data_           (NULLPTR),
n_           (obj.size()),
icsz_ (_G_VEC_CACHE_SIZE),
bdatalocal_        (TRUE)
{
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
  data_ = new T [n_];
#endif
  assert(this->data_!= NULLPTR );
  for ( auto j=0; j<n_; j++ ) {
    data_[j] = obj[j];
  }
  gindex_ = obj.gindex_;
  gindex_keep_ = gindex_;
  updatedev();
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   :
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T>
GTVector<T>::~GTVector()
{
  GEOFLOW_TRACE();
#if defined(GEOFLOW_USE_CUBLAS)
    cudaFree(data_);
#else
  if ( data_  != NULLPTR  && bdatalocal_ ) delete [] data_;
#endif
  data_ = NULLPTR;
}


//**********************************************************************************
//**********************************************************************************
// METHOD : size
// DESC   : Get working (current member) size of data block (not including stride)
// ARGS   : 
// RETURNS: GSIZET total size of data block
//**********************************************************************************
template<class T>
GSIZET GTVector<T>::size() const
{
  GEOFLOW_TRACE();
  return gindex_.end() - gindex_.beg() + 1;
} // end, method size


//**********************************************************************************
//**********************************************************************************
// METHOD : capacity
// DESC   : Get total current capacity of block 
// ARGS   : 
// RETURNS: GSIZET total size of data block
//**********************************************************************************
template<class T>
GSIZET GTVector<T>::capacity() const
{
  return n_;
} // end, method capacity


//**********************************************************************************
//**********************************************************************************
// METHOD : getIndex
// DESC   : Get generalized index 
// ARGS   : 
// RETURNS: GIndex object
//**********************************************************************************
template<class T> 
GIndex &GTVector<T>::getIndex() 
{
  return gindex_;
} // end, method getIndex


//**********************************************************************************
//**********************************************************************************
// METHOD : resize 
// DESC   : Modify gindex if capacity is sufficient; else reallocate
//          to instantiate specified number of objects. No attempt
//          is made here to retain data.
// ARGS   : GSIZET n: new vector size
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTVector<T>::resize(GSIZET nnew)
{
  GEOFLOW_TRACE();
  assert(bdatalocal_ && "Data not local; cannot resize");
  
  GLLONG ibeg    = gindex_.beg();
  GLLONG iend    = ibeg + nnew - 1;
  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  if ( (iend-ibeg+1+ipad) > n_ ) {      // must reallocate; change capacity
    if ( this->data_ != NULLPTR ) { 
#if defined(GEOFLOW_USE_CUBLAS)
      cudaFree(this->data_);
#else
      delete [] this->data_;
#endif
      this->data_ = NULLPTR; 
    }
    this->n_ = iend-ibeg+1+ipad;
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
    this->data_ = new T [this->n_];
#endif
  }
  gindex_(nnew, nnew, ibeg, iend, istride, ipad);
  gindex_keep_ = gindex_;
} // end, method resize


//**********************************************************************************
//**********************************************************************************
// METHOD : resize 
// DESC   : Resize vector/datablock, using default allocators
//          to instantiate specified number (via GIndex object) of usable objects
// ARGS   : GSIZET n: new vector size
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTVector<T>::resize(GIndex &gi)
{
  GEOFLOW_TRACE();
  assert(bdatalocal_ && "Data not local; cannot resize");

  GLLONG nnew    = gi.end()+gi.pad()+1;
  GLLONG ibeg    = gi.beg();
  GLLONG iend    = gi.end();
  GLLONG istride = gi.stride();
  GLLONG ipad    = gi.pad();

  if ( (iend-ibeg+1+ipad) > n_ ) {      // must reallocate; change capacity
    if ( this->data_ != NULLPTR ) { 
#if defined(GEOFLOW_USE_CUBLAS)
      cudaFree(this->data_);
#else
      delete [] this->data_; 
#endif
      this->data_ = NULLPTR; 
    }
    this->n_ = iend-ibeg+1+ipad;
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
    this->data_ = new T [this->n_];
#endif
  }
  gindex_(nnew, nnew, ibeg, iend, istride, ipad);
  gindex_keep_ = gindex_;
} // end, method resize


//**********************************************************************************
//**********************************************************************************
// METHOD : resizem 
// DESC   : Resize only if specified number is > current size
// ARGS   : GSIZET n: new vector size
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTVector<T>::resizem(GSIZET nnew)
{
  GEOFLOW_TRACE();
  if ( nnew > n_ ) {

    assert(bdatalocal_ && "Data not local; cannot resize");

#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [n_];
  }
#else
    resize(nnew);
#endif
  }

} // end, method resizem


//**********************************************************************************
//**********************************************************************************
// METHOD : reserve
// DESC   : Change capacity of data block. Instantiate elements using only
//          default allocator of T. Growing or shrinking, the data buffer contains
//          what was already there (or as much of it as possible).
// ARGS   : GSIZET nnew: new capacity
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTVector<T>::reserve(GSIZET nnew)
{
  GEOFLOW_TRACE();
  assert(bdatalocal_ && "Can reserve only on data local vector");

  T *ttmp=NULLPTR;
  GLLONG ibeg    = gindex_.beg();
  GLLONG iend    = gindex_.end();
  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  // Check: is following exception-safe? No....
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&ttmp, (ibeg+nnew+ipad)*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    ttmp = new T [ibeg+nnew+ipad];
  }
#else
  ttmp  = new T [ibeg+nnew+ipad];
#endif
  assert(ttmp != NULLPTR );

  // Copy old data to temp buffer:
  if ( nnew > n_ ) { // growing
    for ( auto j=0; j<n_; j++ ) ttmp[j] = this->data_[j];
    if ( this->data_ != NULLPTR ){
#if defined(GEOFLOW_USE_CUBLAS)
      cudaFree(this->data_);
#else
      delete [] this->data_;
#endif
      this->data_ = NULLPTR; 
    }
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, (ibeg+nnew+ipad)*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [ibeg+nnew+ipad];
  }
#else
    this->data_ = new T [ibeg+nnew+ipad];
#endif
    assert(this->data_ != NULLPTR );

    // Copy only what was there already to expanded buffer,
    // leaving remainder 'uninitialized':
    for ( auto j=0; j<n_; j++ ) this->data_[j] = ttmp[j];
    gindex_(nnew, nnew, ibeg, iend, istride, ipad);
    n_ = nnew;
  }
  else if ( nnew < n_ ) { // shrinking
    for ( auto j=0; j<nnew; j++ ) ttmp[j] = this->data_[j];
    if ( this->data_ != NULLPTR ) {
#if defined(GEOFLOW_USE_CUBLAS)
      cudaFree(this->data_);
#else
      delete [] this->data_;
#endif
      this->data_ = NULLPTR; 
    }
#if defined(GEOFLOW_USE_CUBLAS)
  if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
    cudaMallocManaged(&data_, (ibeg+nnew+ipad)*sizeof(T), cudaMemAttachGlobal);
  }
  else {
    data_ = new T [ibeg+nnew+ipad];
  }
#else
    this->data_ = new T [ibeg+nnew+ipad];
#endif
    assert(this->data_ != NULLPTR );

    // Copy only what of the original fills fills new buffer:
    n_ = nnew;
    for ( auto j=0; j<n_; j++ ) this->data_[j] = ttmp[j];
    gindex_(n_, n_, ibeg, MIN(n_-1,iend), istride, ipad);
  }

  if ( ttmp != NULLPTR ) delete [] ttmp;
} // end, method reserve


//**********************************************************************************
//**********************************************************************************
// METHOD : clear
// DESC   : Remove all elements, leave with a size/capacity of 0
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTVector<T>::clear()
{
  GEOFLOW_TRACE();
  if ( data_ != NULLPTR ) {
#if defined(GEOFLOW_USE_CUBLAS)
      cudaFree(this->data_);
#else
    delete [] data_;
#endif
    this->data_ = NULLPTR; 
  }
  n_ = 0;
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;
} // end, method clear


//**********************************************************************************
//**********************************************************************************
// METHOD : push_back
// DESC   : add element to end of data block.
// ARGS   : T const &t
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTVector<T>::push_back(const T &obj)
{
  GEOFLOW_TRACE();
  GSIZET nnew = gindex_.end() + gindex_.pad() + 2;

  if ( nnew > n_ ) { // reallocate if required
    reserve(nnew);
  }


  GIndex  gi      = gindex_; 
  GLLONG  iglob   = gindex_.szglobal();
  GLLONG  iloc    = gindex_.szlocal();
  GLLONG  istride = gindex_.stride();
  GLLONG  ipad    = gindex_.pad();
  GLLONG  ibeg    = gindex_.beg();
  GLLONG  iend    = gindex_.end() + 1;
  gindex_(iglob, iloc, ibeg, iend, istride,  ipad);

  data_[gindex_.end()] = obj;

  #if defined(_G_AUTO_CREATE_DEV)
    updatedev();
  #endif

} // end, method push_back


//**********************************************************************************
//**********************************************************************************
// METHOD : back
// DESC   : Get reference to last element
// ARGS   : none
// RETURNS: T reference
//**********************************************************************************
template<class T>
T &GTVector<T>::back()
{
  return data_[n_-1];
} // end, method back


//**********************************************************************************
//**********************************************************************************
// METHOD : back
// DESC   : Get reference to last element
// ARGS   : none
// RETURNS: T reference
//**********************************************************************************
template<class T>
T &GTVector<T>::back() const
{
  return data_[n_-1];
} // end, method back


//**********************************************************************************
//**********************************************************************************
// METHOD : data
// DESC   : Get pointer to data block. 
// ARGS   : 
// RETURNS: T* pointer to data
//**********************************************************************************
template<class T> 
T *GTVector<T>::data() 
{
  return data_ + gindex_.beg();
} // end, method data 


//**********************************************************************************
//**********************************************************************************
// METHOD : data (const)
// DESC   : Get (const) pointer to data block. 
// ARGS   : 
// RETURNS: T* pointer to data
//**********************************************************************************
template<class T> 
const T *GTVector<T>::data() const 
{
  return data_ + gindex_.beg();
} // end, method data (const)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator ==
// DESC   : Equivalence with another GTVector
// ARGS   : GTVector<T> & right-hand side arg 
// RETURNS: TRUE or FALSE
//**********************************************************************************
template<class T> 
GBOOL GTVector<T>::operator==(const GTVector<T> &obj) 
{
  GEOFLOW_TRACE();
  if ( !(this->gindex_ == obj.gindex_) ) return FALSE;

  GBOOL bret = TRUE;
  for ( auto j=gindex_.beg(); j<=gindex_.end() && bret; j++ ) {
    bret = bret && data_[j] == obj[j];
  }

  return bret;
} // end, operator==(GTVector &)


//**********************************************************************************
//**********************************************************************************
// METHOD : assignment operator= GTVector
// DESC   : Equate to another GTVector. If incomping vector is constant,
//          s.t. size == 1, the assign this to entire vector
// ARGS   : GTVector<T> & right-hand side arg 
// RETURNS: GTVector & 
//**********************************************************************************
template<class T>
GTVector<T> &GTVector<T>::operator=(const GTVector<T> &obj)
{
  GEOFLOW_TRACE();
  if ( this == &obj) {
    return *this;
  }

  if ( n_ > 0 ) { // If not allocated, allocate; else it's an error
    assert( (this->n_ < obj.capacity() || this->data_ != NULLPTR) && "L-vector has insufficient size");
  }

  if ( data_ == NULLPTR ) {
    n_ = obj.capacity();
#if defined(GEOFLOW_USE_CUBLAS)
    if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
      cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
    }
    else {
      data_ = new T [n_];
    }
#else
    data_ = new T [n_];
#endif
    assert(this->data_ != NULLPTR );
    gindex_ = obj.gindex_;
    gindex_keep_ = gindex_;
  }

  if ( obj.size() > 1 ) {
    for ( auto j=gindex_.beg(); j<=gindex_.end(); j++ ) {
      data_[j] = obj[j-gindex_.beg()];
    }
  }
  else {
    for ( auto j=gindex_.beg(); j<=gindex_.end(); j++ ) {
      data_[j] = obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

  return *this;
} // end, operator=(GTVector &)


//**********************************************************************************
//**********************************************************************************
// METHOD : assignment operator= std::vector
// DESC   : Equate to std::vector of same template type
// ARGS   : GTVector<T> & right-hand side arg 
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> 
GTVector<T> &GTVector<T>::operator=(const std::vector<T> &obj) 
{
  GEOFLOW_TRACE();

  if ( n_ > 0 ) { // If not allocated, allocate; else it's an error
    assert( (this->n_ < obj.capacity() || this->data_ != NULLPTR) && "L-vector has insufficient size");
  }

  if ( data_ == NULLPTR ) {
    n_ = obj.capacity();
#if defined(GEOFLOW_USE_CUBLAS)
    if ( !std::is_class<T>::value && !std::is_pointer<T>::value && std::is_floating_point<T>::value ) {
      cudaMallocManaged(&data_, n_*sizeof(T), cudaMemAttachGlobal);
    }
    else {
      data_ = new T [n_];
    }
#else
    data_ = new T [n_];
#endif
    assert(this->data_ != NULLPTR );
    gindex_(n_, n_, 0, n_-1, 1,  0);
    gindex_keep_ = gindex_;
  }

  for ( auto j=gindex_.beg(); j<=gindex_.end(); j++ ) {
    data_[j] = obj[j-gindex_.beg()];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
  
  return *this;
} // end, operator=(std::vector &)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator=
// DESC   : Set vector to a constant
// ARGS   : T aarg
// RETURNS: none
//**********************************************************************************
template<class T>
void GTVector<T>::operator=(T a)
{
  GEOFLOW_TRACE();
  for ( auto j=gindex_.beg(); j<=gindex_.end(); j+=gindex_.stride() ) {
    data_[j] = a;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end operator=


//************************************************************************************
//************************************************************************************
// METHOD : range 
// DESC   : Sets GTVector<T> range ibeg, iend, by
//          re-defining variable gindex_ object.
// ARGS   : ibeg : starting index
//          iend : ending index
// RETURNS: none.
//************************************************************************************
template<class T> 
void  GTVector<T>::range(GLONG ibeg, GLONG iend)
{
  GEOFLOW_TRACE();
//assert(iend < n_ && ibeg <= iend && "Invalid range specification");
  if (  ibeg <= iend
   &&  (ibeg >= static_cast<GLONG>(n_) 
   ||   iend >= static_cast<GLONG>(n_)) ) {
    std::cout << "GTVector::range: invalid range specification: n_=" << n_ << " ibeg=" << ibeg << " iend=" << iend << std::endl;
    while (1);
  }

  gindex_.beg() = ibeg;
  gindex_.end() = iend;
  
} // end of method range


//************************************************************************************
//************************************************************************************
// METHOD : range_reset
// DESC   : Resets range to that stored in gindex_keep object.
//          Capacity, etc does not change.
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<class T> 
void  GTVector<T>::range_reset()
{
  GEOFLOW_TRACE();
  gindex_ = gindex_keep_;
} // end of method range_reset


//**********************************************************************************
//**********************************************************************************
// METHOD : set
// DESC   : Set vector to input constant
// ARGS   : T-type 
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTVector<T>::set(T a)
{ 
  GEOFLOW_TRACE();
  for ( auto j=gindex_.beg(); j<=gindex_.end(); j+=gindex_.stride() ) {
    data_[j] = a;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
} // end of method set


//**********************************************************************************
//**********************************************************************************
// METHOD : set
// DESC   : Set vector to input elements up to max specified
// ARGS   : T-type array 
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTVector<T>::set(T *a, GSIZET n)
{ 
  GEOFLOW_TRACE();
#if !defined(GEOFLOW_USE_OPENACC)
  ASSERT(!(gindex_.beg()+n > gindex_.end()+1));
#endif
  
  GLLONG m=0;
  for ( auto j=gindex_.beg(); j<MIN(gindex_.beg()+n,gindex_.end()+1) && m < n; j+=gindex_.stride() ) {
    data_[j] = a[m++];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
} // end of method set


//**********************************************************************************
//**********************************************************************************
// METHOD : floor
// DESC   : Set values < floor to specified floor:
//               element = MAX(element, floor)
// ARGS   : T-type 
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTVector<T>::floor(T a)
{ 
  GEOFLOW_TRACE();
  for ( auto j=gindex_.beg(); j<=gindex_.end(); j+=gindex_.stride() ) {
    data_[j] = MAX(data_[j],a);
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
} // end of method floor


//**********************************************************************************
//**********************************************************************************
// METHOD : updatehost
// DESC   : Update data from device to host
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<class T> 
void GTVector<T>::updatehost()
{
  GEOFLOW_TRACE();
} // end of method updatehost


//************************************************************************************
//************************************************************************************
// METHOD : updatedev
// DESC   : Update data from host to device
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<class T> 
void GTVector<T>::updatedev()
{
  GEOFLOW_TRACE();
} // end of method updatedev

//**********************************************************************************
//**********************************************************************************
// METHOD : transpose 
// DESC   : Imposes additional structure on array data, by
//          computing in-place 'transpose', looking at 
//          only first n^2 elements as though they're arraged
//          in a matrix. It is assumed that the matrix size
//          is nxn.
// ARGS   : n : row and column length imposted as matrix
// RETURNS: none
//**********************************************************************************
template<class T> void GTVector<T>::transpose(GSIZET n)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: GTVector<T>::transpose()");

  T       tmp;

  for ( auto j=0; j<n; j++ ) {
    for ( auto i=j; i<n; i++ ) {
       tmp = (*this)[i+n*j];
       (*this)[i+j*n] = (*this)[j+i*n];
       (*this)[j+i*n] = tmp;
    }
  }
} // end, method transpose

//**********************************************************************************
//**********************************************************************************
// METHOD : operator +
// DESC   : addition of vector and constant, create new vector
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::operator+(const T b)
{
  GEOFLOW_TRACE();
    GTVector ans(*this);
    ans += b;
    return ans;
} // end, operator+

//**********************************************************************************
//**********************************************************************************
// METHOD : operator -
// DESC   : subtraction of vector and constant, create new vector
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::operator-(const T b)
{
  GEOFLOW_TRACE();
    GTVector ans(*this);
    ans -= b;
    return ans;
} // end, operator-

//**********************************************************************************
//**********************************************************************************
// METHOD : operator *
// DESC   : product of vector and constant, create new vector
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::operator*(const T b)
{
  GEOFLOW_TRACE();
    GTVector ans(*this);
    ans *= b;
    return ans;
} // end, operator*

//**********************************************************************************
//**********************************************************************************
// METHOD : operator +
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector<T>&
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::operator+(const GTVector &obj)
{
  GEOFLOW_TRACE();
  return this->add_impl_(obj, typename std::is_floating_point<T>::type());
} // end, operator+


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector<T>&
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::operator-(const GTVector &obj)
{
  GEOFLOW_TRACE();
  return this->sub_impl_(obj, typename std::is_floating_point<T>::type());
} // end, operator-

//**********************************************************************************
//**********************************************************************************
// METHOD : operator *
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector<T>&
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::operator*(const GTVector &obj)
{
  GEOFLOW_TRACE();
  return this->mul_impl_(obj, typename std::is_floating_point<T>::type());
} // end, operator*


//**********************************************************************************
//**********************************************************************************
// METHOD : dot
// DESC   :
// ARGS   : GTVector &
// RETURNS: typename T
//**********************************************************************************
template<class T>
T GTVector<T>::dot(const GTVector &obj)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value &&
    "Invalid template type: GVector<T>::dot()");

  T ret = 0;
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    ret += this->data_[j] * obj[j];
  }

  return ret;

} // end, dot


//**********************************************************************************
//**********************************************************************************
// METHOD : gdot (1)
// DESC   : Compute dot product over all ranks
// ARGS   : obj  : input vector
//          comm : communicator
// RETURNS: typename T
//**********************************************************************************
template<class T>
T GTVector<T>::gdot(const GTVector &obj, GC_COMM comm)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value &&
    "Invalid template type: GVector<T>::gdot()");

  T lret = this->dot(obj); 
  T gret;

  GComm::Allreduce(&lret, &gret, 1, T2GCDatatype<T>() , GC_OP_SUM, comm);

  return gret;

} // end, gdot (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : gdot (2)
// DESC   : Compute dot product over all ranks:
//              (this:b)^T . c,
//          where : represents pointProd of operands
// ARGS   : b : input vector
//          c : input vector
//          comm : communicator
// RETURNS: typename T
//**********************************************************************************
template<class T>
T GTVector<T>::gdot(const GTVector &b, const GTVector &c, GC_COMM comm)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value &&
    "Invalid template type: GVector<T>::gdot()");

  T lret=0.0; 
  T gret;

  for ( auto j=gindex_.beg(); j<=gindex_.end(); j+=gindex_.stride() ) {
    lret += data_[j]*b[j]*c[j];
  }
  
  GComm::Allreduce(&lret, &gret, 1, T2GCDatatype<T>() , GC_OP_SUM, comm);

  return gret;

} // end, gdot(2)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator += (2)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T>
void
GTVector<T>::operator+=(const T b)
{
  GEOFLOW_TRACE();
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    this->data_[j] += b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator+= (2)

//**********************************************************************************
//**********************************************************************************
// METHOD : operator -= (2)
// DESC   :
// ARGS   : T arg
// RETURNS: void
//**********************************************************************************
template<class T>
void
GTVector<T>::operator-=(const T b)
{
  GEOFLOW_TRACE();
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    this->data_[j] -= b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operatori= (2)

//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (2)
// DESC   : product of vector and constant
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
void
GTVector<T>::operator*=(const T b)
{
  GEOFLOW_TRACE();
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    this->data_[j] *= b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator*= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator += 
// DESC   : Self addition, for non-float types
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::operator+=(const GTVector &obj)
{
  GEOFLOW_TRACE();
  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] += obj[j-this->gindex_.beg()];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] += obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator+= 


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -= 
// DESC   : Self-subtraction, for non-float types
// ARGS   : GTVector & arg
// RETURNS: void
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::operator-=(const GTVector &obj)
{
  GEOFLOW_TRACE();
  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] -= obj[j-this->gindex_.beg()];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] -= obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator-= 


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (1)
// DESC   : point-by-point product. If incoming factor is constant,
//          s.t size==1, then only first element is used
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
void
GTVector<T>::operator*=(const GTVector &obj)
{
  GEOFLOW_TRACE();
  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] *= obj[j-this->gindex_.beg()];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] *= obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator*= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : operator /= (1)
// DESC   : point-by-point division: this = this/b
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
void
GTVector<T>::operator/=(const GTVector &obj)
{
  GEOFLOW_TRACE();
  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] /= obj[j-this->gindex_.beg()];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      this->data_[j] /= obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator/= (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : isfinite (1)
// DESC   : Find if a nan or inf exists in vector
// ARGS   : 
// RETURNS: TRUE if ok (all members are finite);  else FALSE
//**********************************************************************************
template<class T>
GBOOL GTVector<T>::isfinite()
{
  GEOFLOW_TRACE();
  GBOOL bret = TRUE;

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end() && j<=this->gindex_.end() && bret; j+=this->gindex_.stride() ) {
    bret = std::isfinite(this->data_[j]);
  }
 
  return bret;

} // end isfinite (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : isfinite (2)
// DESC   : Find if a nan or inf exists in vector, and provide index
//          of first occurrence.
// ARGS   : iwhere: index of first occurrence; set only if isfinite == FALSE
// RETURNS: TRUE if ok (all members are finite);  else FALSE
//**********************************************************************************
template<class T>
GBOOL GTVector<T>::isfinite(GSIZET &iwhere)
{
  GEOFLOW_TRACE();
  GBOOL bret = TRUE;

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end() && j<=this->gindex_.end() && bret; j+=this->gindex_.stride() ) {
    bret = std::isfinite(this->data_[j]);
    if ( !bret ) iwhere = j;
  }
 
  return bret;

} // end isfinite (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : maxn
// DESC   : Find max of first n elements
// ARGS   : n : num elemets past gindex.beg() to check
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::maxn(GSIZET n)
{
  GEOFLOW_TRACE();
  T fm = std::numeric_limits<T>::min();

  for ( auto j=this->gindex_.beg(); j<this->gindex_.beg()+n && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MAX(fm,this->data_[j]);
  }
 
  return fm;

} // end maxn


//**********************************************************************************
//**********************************************************************************
// METHOD : imax
// DESC   : Find position of max of vector
// ARGS   : none.
// RETURNS: GSIZET position of max
//**********************************************************************************
template<class T>
GSIZET
GTVector<T>::imax()
{
  GEOFLOW_TRACE();
  GSIZET imax;
  T fm = std::numeric_limits<T>::min();

  for ( auto j=this->gindex_.beg(); j<this->gindex_.end() && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    if ( this->data_[j] > fm ) {
      imax= j;
      fm = this->data_[j];
    }
  }
 
  return imax;

} // end imax


//**********************************************************************************
//**********************************************************************************
// METHOD : max
// DESC   : Find max of all members
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::max()
{
  GEOFLOW_TRACE();
  T fm = std::numeric_limits<T>::min();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MAX(fm,this->data_[j]);
  }
 
  return fm;

} // end max


//**********************************************************************************
//**********************************************************************************
// METHOD : gmax
// DESC   : Find global max of all members (over all MPI ranks)
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::gmax(GC_COMM comm)
{
  GEOFLOW_TRACE();
  T lm = std::numeric_limits<T>::min();
  T gm;

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    lm = MAX(lm,this->data_[j]);
  }
  GComm::Allreduce(&lm, &gm, 1, T2GCDatatype<T>() , GC_OP_MAX, comm);
 
  return gm;

} // end gmax


//**********************************************************************************
//**********************************************************************************
// METHOD : amax
// DESC   : Find absolute max of all members
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::amax()
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T fm = std::numeric_limits<T>::min();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MAX(fm,std::fabs(this->data_[j]));
  }

  return fm;

} // end amax


//**********************************************************************************
//**********************************************************************************
// METHOD : gamax
// DESC   : Find global |max| of all members (over all MPI ranks)
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::gamax(GC_COMM comm)
{
  GEOFLOW_TRACE();
  T lm = std::numeric_limits<T>::min();
  T gm;

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    lm = MAX(lm,fabs(this->data_[j]));
  }
  GComm::Allreduce(&lm, &gm, 1, T2GCDatatype<T>() , GC_OP_MAX, comm);
 
  return gm;

} // end gamax


//**********************************************************************************
//**********************************************************************************
// METHOD : amaxdiff
// DESC   : Find max of absolute sequential differences > 0
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::amaxdiff(T tiny)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T diff;
  T fm = std::numeric_limits<T>::min();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end()-1; j+=this->gindex_.stride() ) {
    diff = fabs(this->data_[j+1] - this->data_[j]);
    if ( diff > tiny ) fm = MAX(fm,diff);
  }

  return fm;

} // end amaxdiff


//**********************************************************************************
//**********************************************************************************
// METHOD : minn
// DESC   : Find max of first n elements
// ARGS   : n : num elements past gindex.beg to check
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::minn(GSIZET n)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T fm = std::numeric_limits<T>::max();

  for ( auto j=this->gindex_.beg(); j<this->gindex_.beg()+n && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MIN(fm,this->data_[j]);
  }

  return fm;

} // end minn


//**********************************************************************************
//**********************************************************************************
// METHOD : imin
// DESC   : Find position of min of vector
// ARGS   : none.
// RETURNS: GSIZET position of min
//**********************************************************************************
template<class T>
GSIZET
GTVector<T>::imin()
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  GSIZET imin;
  T fm = std::numeric_limits<T>::max();

  for ( auto j=this->gindex_.beg(); j<this->gindex_.end() && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    if ( this->data_[j] < fm ) {
      imin= j;
      fm = this->data_[j];
    }
  }

  return imin;

} // end imin


//**********************************************************************************
//**********************************************************************************
// METHOD : min
// DESC   : Find max of all members
// ARGS   : none.
// RETURNS: T-type min
//**********************************************************************************
template<class T>
T
GTVector<T>::min()
{
  GEOFLOW_TRACE();
  T fm = std::numeric_limits<T>::max();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MIN(fm,this->data_[j]);
  }

  return fm;

} // end min


//**********************************************************************************
//**********************************************************************************
// METHOD : gmin
// DESC   : Find global min of all members (over all MPI ranks)
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::gmin(GC_COMM comm)
{
  GEOFLOW_TRACE();
  T lm = std::numeric_limits<T>::max();
  T gm;

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    lm = MIN(lm,this->data_[j]);
  }
  GComm::Allreduce(&lm, &gm, 1, T2GCDatatype<T>() , GC_OP_MIN, comm);
 
  return gm;

} // end gmin


//**********************************************************************************
//**********************************************************************************
// METHOD : gamin
// DESC   : Find global |min| of all members (over all MPI ranks)
// ARGS   : none.
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::gamin(GC_COMM comm)
{
  GEOFLOW_TRACE();
  T lm = std::numeric_limits<T>::max();
  T gm;

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    lm = MIN(lm,fabs(this->data_[j]));
  }
  GComm::Allreduce(&lm, &gm, 1, T2GCDatatype<T>() , GC_OP_MIN, comm);
 
  return gm;

} // end gamin


//**********************************************************************************
//**********************************************************************************
// METHOD : amin
// DESC   : Find absolute min of all members
// ARGS   : none.
// RETURNS: T-type min
//**********************************************************************************
template<class T>
T
GTVector<T>::amin()
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T fm = std::numeric_limits<T>::max();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MIN(fm,fabs(this->data_[j]));
  }

  return fm;

} // end amin


//**********************************************************************************
//**********************************************************************************
// METHOD : amindiff
// DESC   : Find min of absolute sequential differences > 0
// ARGS   : none.
// RETURNS: T-type min
//**********************************************************************************
template<class T>
T
GTVector<T>::amindiff(T tiny)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T diff;
  T fm = std::numeric_limits<T>::max();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end()-1; j+=this->gindex_.stride() ) {
    diff = fabs(this->data_[j+1] - this->data_[j]);
    if ( diff > tiny ) fm = MIN(fm,diff);
  }

  return fm;

} // end amindiff


//**********************************************************************************
//**********************************************************************************
// METHOD : pointProd (1)
// DESC   : point-by-point multiplication, returned in specified GTVector
// ARGS   : obj: const GTVector<>  factor
//          ret: GTVector & ret
// RETURNS: GTVector & 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::pointProd(const GTVector<T> &obj, GTVector<T> &ret ) 
{
  GEOFLOW_TRACE();
  ASSERT(!( (obj.size() > 1 && obj.size() < this->size()) || ret.size() < this->size() ));

  if ( obj.size() > 1 ) {
#pragma acc parallel loop
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      ret[j] = this->data_[j-gindex_.beg()] * obj[j-gindex_.beg()];
    }
  }
  else {
#pragma acc parallel loop
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      ret[j] = this->data_[j-gindex_.beg()] * obj[0];
    }
  }

} // end, pointProd (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : pointProd (2)
// DESC   : point-by-point multiplication, and additional constant
//          factor returned in specified GTVector
//             ret = a * (*this)*obj
// ARGS   : a  : constant factor
//          obj: const GTVector<>  factor
//          ret: GTVector & ret
// RETURNS: GTVector & 
//**********************************************************************************
template<class T>
void
GTVector<T>::pointProd(const T a, const GTVector<T> &obj, GTVector<T> &ret ) 
{
  GEOFLOW_TRACE();
  ASSERT( (ret.size() >= this->size()) && (obj.size() >= this->size()) );

  if ( obj.size() > 1 ) {
#pragma acc parallel loop
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      ret[j] = a * this->data_[j-gindex_.beg()] * obj[j-gindex_.beg()];
    }
  }
  else {
#pragma acc parallel loop
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      ret[j] = a * this->data_[j-gindex_.beg()] * obj[0];
    }
  }

} // end, pointProd (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : pointProd (3)
// DESC   : point-by-point multiplication, 
//            this_j *= obj_j
//          returned in *this
// ARGS   : obj: const GTVector<>  factor
// RETURNS: none
//**********************************************************************************
template<class T>
void
GTVector<T>::pointProd(const GTVector<T> &obj)
{
  GEOFLOW_TRACE();
  ASSERT(!( obj.size() > 1 && obj.size() < this->size() ));

  if ( obj.size() > 1 ) {

#pragma acc parallel loop
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      data_[j] *= obj[j-gindex_.beg()];
    }
  }
  else {
#pragma acc parallel loop
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      data_[j] *= obj[0];
    }
  }

} // end, pointProd (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : apointProd 
// DESC   : point-by-point multiplication, and additional constant
//          factor returned in specified GTVector
//             this *= a * obj
// ARGS   : a  : constant factor
//          obj: const GTVector<>  factor
// RETURNS: GTVector & 
//**********************************************************************************
template<class T>
void
GTVector<T>::apointProd(const T a, const GTVector<T> &obj ) 
{
  GEOFLOW_TRACE();
  ASSERT(!( (obj.size() != this->size()) && obj.size() != 1 ));

  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      (*this)[j] *= a * obj[j];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      (*this)[j] *= a * obj[0];
    }
  }

} // end, apointProd


//**********************************************************************************
//**********************************************************************************
// METHOD : constProd (1)
// DESC   : multiply this by constant, return in ret vector
// ARGS   : b  : T-type constant 
//          ret: GTVector & ret vector
// RETURNS: GTVector & 
//**********************************************************************************
template<class T>
void
GTVector<T>::constProd(const T b, GTVector<T> &ret) 
{
  GEOFLOW_TRACE();
  ASSERT(!(ret.size() < this->size()));

  T *dret=ret.data();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    dret[j] = this->data_[j] * b;
  }

} // end, pointProd


//**********************************************************************************
//**********************************************************************************
// METHOD : constProd (2)
// DESC   : multiply this by constant, at specified indices
// ARGS   : b  : T-type constant 
//          iv : indirection indices at which to do product
//          nv : number of indices in iv
// RETURNS: GTVector & 
//**********************************************************************************
template<class T>
void
GTVector<T>::constProd(const T b, GSIZET *iv, GSIZET nv) 
{
  GEOFLOW_TRACE();
  ASSERT(iv != NULLPTR);

  for ( auto j=0; j<nv; j++ ) {
    this->data_[iv[j]] *= b;
  }

} // end, pointProd (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : sum 
// DESC   : Sum vector contents and return
// ARGS   : none.
// RETURNS: T sum
//**********************************************************************************
template<class T>
T
GTVector<T>::sum()
{
  GEOFLOW_TRACE();
  T      sum=static_cast<T>(0);

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    sum +=  this->data_[j];
  }

  return sum;
} // end, sum


//**********************************************************************************
//**********************************************************************************
// METHOD : sum 
// DESC   : Sum vector contents and return
// ARGS   : ibeg, iend: starting and stopping index
// RETURNS: T sum
//**********************************************************************************
template<class T>
T
GTVector<T>::sum(GLLONG ibeg, GLLONG iend) 
{
  GEOFLOW_TRACE();
  T      sum=static_cast<T>(0);
  assert(ibeg >= this->gindex_.beg() && iend <= this->gindex_.end());
  for ( auto j=ibeg; j<=iend; j+=this->gindex_.stride() ) {
    sum +=  this->data_[j]; 
  }

  return sum;
} // end, sum


//**********************************************************************************
//**********************************************************************************
// METHOD : infnorm
// DESC   : Computes infinity (max)
// ARGS   : none.
// RETURNS: T norm
//**********************************************************************************
template<class T>
T
GTVector<T>::infnorm() 
{
  GEOFLOW_TRACE();
  T xnorm = std::numeric_limits<T>::min();

  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    xnorm = MAX(xnorm,fabs(this->data_[j]));
  }
  
  return xnorm;
} // end, infnorm


//**********************************************************************************
//**********************************************************************************
// METHOD : Eucnorm
// DESC   : Computes Euclidean (RMS) norm
// ARGS   : none.
// RETURNS: T norm
//**********************************************************************************
template<class T>
T
GTVector<T>::Eucnorm() 
{
  GEOFLOW_TRACE();
  GDOUBLE n, xnorm=0.0;

  n = 0.0;
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    xnorm += this->data_[j]*this->data_[j];
    n += 1.0;
  }
  xnorm /= n;
  
  return static_cast<T>(sqrt(xnorm));
} // end, Eucnorm


//**********************************************************************************
//**********************************************************************************
// METHOD : rpow
// DESC   : Raise vector elems to power p, and modify array contents
// ARGS   : p: power
// RETURNS: none
//**********************************************************************************
template<class T>
void
GTVector<T>::rpow(const GDOUBLE p)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  GDOUBLE b;
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    b        = static_cast<GDOUBLE>(data_[j]);
    data_[j] = static_cast<T>(pow(b, p));
  }

} // end, rpow


//**********************************************************************************
//**********************************************************************************
// METHOD : abs
// DESC   : Replace each elem with absolute value
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<class T>
void
GTVector<T>::abs()
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
//  tmp = sqrt( std::pow<GDOUBLE>(static_cast<GDOUBLE>(data_[j]),2) );
    data_[j] = std::fabs(data_[j]);
  }

} // end, abs


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity
// DESC   : Gets multiplicity of specified value
// ARGS   : val : type T
// RETURNS: GSIZET multiplicity
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity(T val)
{
	GEOFLOW_TRACE();
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
//  "Invalid template type: multiplicity(T)");

  GSIZET mult=0;


  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] == val ) {
      mult++;
    }
  }
  
  return mult;

} // multiplicity


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity_s
// DESC   : Gets multiplicity of specified value, assuming data has been sorted.
//          DO NOT use if it has not been!

// ARGS   : val   : type T
//          istart: index at which to start search
//          ifound: index where 'val' is first encountered, if at all
// RETURNS: GSIZET multiplicity
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_s(T val, GSIZET istart, GSIZET &ifound)
{
	GEOFLOW_TRACE();
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
//  "Invalid template type: multiplicity(T)");

  GLLONG i;
  GSIZET mult=0;

  assert(istart >= this->gindex_.beg() && istart <= this->gindex_.end() );

  i = istart;
  while ( val !=  this->data_[i] && i <= this->gindex_.end() ) {
    i += this->gindex_.stride();
  }
  ifound = i;
  while ( val ==  this->data_[i] && i <= this->gindex_.end() ) {
    mult++;
    i += this->gindex_.stride();
  }
  
  return mult;

} // multiplicity_s


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity
// DESC   : Gets multiplicity of specified value, and
//          return indices of occurrences
// ARGS   : val  : type T
//          index: Indices where value occurs in array.
//                 Re-allocated if multiplicity > n
//          n    : new size of index array, if resized. Should be >= multiplicity
// RETURNS: GSIZET multiplicity
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity(T val, GSIZET *&index, GSIZET &n)
{
	GEOFLOW_TRACE();
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
//  "Invalid template type: multiplicity(T,GSIZET*,GSIZET&)");

  GLLONG m=0;
  GSIZET mult=0;

  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] == val ) {
      mult++;
    }
  }

  if ( n < mult ) {
    if ( index != NULLPTR ) delete [] index;
    n = mult;
    index = new GSIZET [n];
  }
  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] == val ) {
      index[m++] = i;
    }
  }
  
  return mult;

} // multiplicity


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity_floor
// DESC   : Gets multiplicity of specified value, for values > floor
// ARGS   : val  : type T
//          floor: check only values > floor
// RETURNS: GSIZET multiplicity
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_floor(T val, T floor)
{
	GEOFLOW_TRACE();
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
//  "Invalid template type: multiplicity_floor(T,T)");

  GSIZET mult=0;

  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] > floor && this->data_[i] == val ) {
      mult++;
    }
  }
  
  return mult;

} // multiplicity_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity_floor
// DESC   : Gets multiplicity of specified value, for values > floor, and
//          return indices of occurrences
// ARGS   : val  : type T
//          index: Re-allocated if multiplicity > n
//          n    : new size of index array, if resized. Should be >= multiplicity
//          floor: check only values > floor
// RETURNS: GSIZET multiplicity
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_floor(T val, GSIZET *&index, GSIZET &n, T floor)
{
	GEOFLOW_TRACE();
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
//  "Invalid template type: multiplicity_floor(T,GSIZET *,GSIZE&,T)");

  GLLONG  m=0;
  GSIZET mult=0;

  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] > floor && this->data_[i] == val ) {
      mult++;
    }
  }

  if ( n < mult ) {
    if ( index != NULLPTR ) delete [] index;
    n = mult;
    index = new GSIZET [n];
  }
  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] > floor && this->data_[i] == val ) {
      index[m++] = i;
    }
  }
  
  return mult;

} // multiplicity_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity_ceil
// DESC   : Gets multiplicity of specified value, for values < ceil
// ARGS   : val  : type T
//          ceil : check only values < ceil
// RETURNS: GSIZET multiplicity
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_ceil(T val, T ceil)
{
	GEOFLOW_TRACE();
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
//  "Invalid template type: multiplicity_ceil(T,T)");

  GSIZET index, mult=1;

  if ( !this->contains(val,index) ) return 0 ;

  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] < ceil && this->data_[i] == val ) {
      mult++;
    }
  }
  
  return mult;

} // multiplicity_ceil


//**********************************************************************************
//**********************************************************************************
// METHOD : onlycontains 
// DESC   : Determines if only candidate value is in the buffer
// ARGS   : val : member to search for in buffer
// RETURNS: TRUE if member is in list, else FALSE. If list is empty, returns TRUE
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::onlycontains(T val)
{
  GEOFLOW_TRACE();
  if ( this->data_ == NULLPTR ) return TRUE;

  GLLONG i=this->gindex_.beg();

  while ( i <= this->gindex_.end() && this->data_[i] == val ) i++;

  if ( i > this->gindex_.end() ) return TRUE;

  return FALSE;

} // end of method onlycontains 


//**********************************************************************************
//**********************************************************************************
// METHOD : findfirst
// DESC   : Find index of first occurrence of val in array, else -1
// ARGS   : val : member to search for in buffer
// RETURNS: index of first occurrence, else -1
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GLONG
GTVector<T>::findfirst(T val)
{
  GEOFLOW_TRACE();
  if ( this->data_ == NULLPTR ) return -1;

  GLLONG i=this->gindex_.beg();

  while ( i <= this->gindex_.end() && this->data_[i] != val ) i+=this->gindex_.stride();

  if ( i > this->gindex_.end() ) return i = -1;

  return i;

} // end of method findfirst

//**********************************************************************************
//**********************************************************************************
// METHOD : findlast
// DESC   : Find index of last occurrence of val in array, else -1
// ARGS   : val : member to search for in buffer
// RETURNS: index of last occurrence, else -1
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GLONG
GTVector<T>::findlast(T val)
{
  GEOFLOW_TRACE();
  if ( this->data_ == NULLPTR ) return -1;

  GLLONG i=this->gindex_.end();

  while ( i >= this->gindex_.beg() && this->data_[i] != val ) i-=this->gindex_.stride();

  if ( i < this->gindex_.beg() ) return i = -1;

  return i;

} // end of method findlast 


//**********************************************************************************
//**********************************************************************************
// METHOD : contains (1)
// DESC   : Determines if candidate value is in the buffer
// ARGS   : val : member to search for in buffer
//          iwhere : first index where val is found first. Valid only
//                   if return is TRUE
// RETURNS: TRUE if member is in list, else FALSE. 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains(T val, GSIZET &iwhere)
{
  GEOFLOW_TRACE();
#if 0
  assert(std::is_arithmetic<T>::value || std::is_string<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: contains(T, GSIZET&)");
#endif

  if ( this->data_ == NULLPTR ) return FALSE;

  GLLONG i=this->gindex_.beg();

  while ( i <= this->gindex_.end() && this->data_[i] != val ) i++;

  if ( i > this->gindex_.end() ) return FALSE;

  iwhere = i;

  return TRUE;

} // end of method contains (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : contains (2)
// DESC   : Determines if candidate value is in the buffer
// ARGS   : val : member to search for in buffer
//          iwhere : array returned with all indices where val was found
//          nw     : current size of iwhere array. If array is too small,
//                   it will be resized to accommodate the new number
//                   of indices where val exsts, and nw will change
// RETURNS: number of indices at which this == val
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::contains(T val, GSIZET *&iwhere, GSIZET &nw)
{
  GEOFLOW_TRACE();
#if 0
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: contains(T, GSIZET&)");
#endif

  if ( this->data_ == NULLPTR ) return FALSE;

  GLLONG n=0;
  
  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i++ ) {
    n += this->data_[i] == val ? 1 : 0;
  }

  if ( n == 0 ) {
     return 0;
  }

  if ( nw < n ) {
    if ( iwhere != NULLPTR ) delete [] iwhere;
    iwhere = new GSIZET [n];
    nw = n;
  }

  n = 0;
  for ( auto i=this->gindex_.beg(); i<=this->gindex_.end(); i++ ) {
    if ( data_[i] == val ) {
      iwhere[n++] = i;
    }
  }

  return n;

} // end of method contains (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : containsn (1)
// DESC   : Determines if candidate value is in the buffer, checking n elements
// ARGS   : val    : member to search for in buffer
//          n      : # elements to check
//          iwhere : first index where val is found; else unchanged
// RETURNS: TRUE if member is in list, else FALSE. 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::containsn(T val, GSIZET n, GSIZET &iwhere)
{
  GEOFLOW_TRACE();
#if 0
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: containsn(T, GSIZET, GSIZET&)");
#endif

  if ( data_ == NULLPTR  ) return FALSE;

  GLLONG i=this->gindex_.beg();

  while ( i < MIN(this->gindex_.beg()+n,n_) &&  data_[i] != val ) i++;

  if ( i >= this->gindex_.beg()+n ) return FALSE;

  iwhere = i;

  return TRUE;

} // end of method containsn (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : contains (2)
// DESC   : Determines if candidate value is in the buffer
// ARGS   : val : member to search for in buffer
// RETURNS: TRUE if member is in list, else FALSE. 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains(T val)
{
  GEOFLOW_TRACE();
#if 0
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: contains(T)");
#endif

  GSIZET nd;

  return contains(val,nd);

} // end of method contains (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : containsn  (2)
// DESC   : Determines if candidate value is in the buffer, searching only the
//          first n elements (after & including gindex.beg)
// ARGS   : val : member to search for in buffer
//          n   : number of indices to check, staring with gindex.beg()
// RETURNS: TRUE if member is in list, else FALSE. 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::containsn(T val, GSIZET n)
{
  GEOFLOW_TRACE();
#if 0
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: containsn(T, GSIZET)");
#endif

  if ( data_ == NULLPTR ) return FALSE;

  GLLONG i=this->gindex_.beg();

  while ( i < MIN(this->gindex_.beg()+n,n_) && data_[i] != val ) i++;

  return i < this->gindex_.beg()+n;

} // end of method containsn (2) 


//**********************************************************************************
//**********************************************************************************
// METHOD : contains_floor (1)
// DESC   : Determines if candidate member is in the buffer, and if so,
//          sets the index where it is found. Search continues until index
//          is found s.t. data >= floor.
// ARGS   : val    : member to search for in buffer
//          iwhere : first index where val is found first; else unchanged
//          floor  : floor value
//          istart : starting index after index.ibeg() (default is 0)
// RETURNS:  TRUE if member is in list, else FALSE. On TRUE, set iwhere to
//           index where member is found first
//************************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains_floor(T val, GSIZET  &iwhere, T floor, GSIZET istart)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value &&
    "Invalid template type: contains_floor");

  if ( data_ == NULLPTR ) return FALSE;

  GSIZET i;
  for ( i=this->gindex_.beg()+istart; i<=this->gindex_.end(); i++ ) {
    if ( data_[i] > floor && data_[i] == val ) break;
  }

  if ( i > this->gindex_.end() || data_[i] <= floor ) return FALSE;

  iwhere = i;

  return TRUE;

} // end of method contains_floor (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : contains_ceil (1)
// DESC   : Determines if candidate member is in the buffer, and if so,
//          sets the index where it is found. Search continues until index
//          is found s.t. data <= ceiling.
// ARGS   : val   : member to search for in buffer
//          iwhere: first index where val is found first; else unchanged
//          ceil  : ceiling value
//          istart : starting index after index.ibeg() (default is 0)
// RETURNS:  TRUE if member is in list, else FALSE. On TRUE, set iwhere to
//           index where member is found first
//************************************************************************************
//#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains_ceil(T val, GSIZET  &iwhere, T ceil, GSIZET istart)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value &&
    "Invalid template type: contains_ceil");

  if ( this->data_ == NULLPTR ) return FALSE;

  GSIZET i;
  for ( i=this->gindex_.beg()+istart; i<=this->gindex_.end(); i++ ) {
    if ( this->data_[i] < ceil && this->data_[i] == val ) break;
  }

  if ( i > this->gindex_.end() || this->data_[i] >= ceil ) return FALSE;

  iwhere = i;

  return TRUE;

} // end of method contains_ceil (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrng 
// DESC   : From buffer, within range given by [ibeg,iend], get indices where there 
//          are distinct values. 
//          Use apariingly, as this is slow, and unoptimized.
// ARG    : 
//          ibeg    : beginning of range (after gindexd_.beg())
//          n       : number of elements to scan
//          is      : stride count for search
//          vals    : array of values of size nd if nd > number of distinct values
//                    found; else = nd;
//          indices : GTVector containing indices into the buffer which contain
//                    the distinct values. Must be dimensioned >= any possible nd,
//                    as no checking is done. This array is deleted and re-allocated.
//          nd      : number of distinct values. nd=0 is not considered an error
//          tunique : tmp array of type T, of size at least of size()
//          itmp    : tmp array of type GSIZET, of size at least of size()
// RETURNS    :  TRUE on success; else FALSE 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng(GSIZET ibeg, GSIZET n, GSIZET is, T *&vals,
                         GSIZET *&indices, GSIZET  &nd, T * const &tunique, GSIZET * const &itmp)
{
  GEOFLOW_TRACE();
  GLLONG nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( auto i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    bcont = FALSE;
    for ( auto j=0; j<nfound && !bcont;j++ ) {
      bcont = this->data_[i] == tunique[j];
    }

    if ( !bcont ) { // if not in unique buff, add it:
      itmp   [nfound] = i;
      tunique[nfound] = this->data_[i];
      nfound++;
    }
  }
  if ( nd < nfound ) {
    if ( vals    != NULLPTR ) delete [] vals;   vals = NULLPTR;
    if ( indices != NULLPTR ) delete [] indices; indices = NULLPTR;
    nd = nfound;
    vals    = new T      [nd];
    assert(vals != NULLPTR );
    indices = new GSIZET [nd];
    assert(indices != NULLPTR );
  }
  for ( auto i=0; i<nfound; i++ ) {
    vals   [i] = tunique[i];
    indices[i] = itmp[i];
  }

  return nfound;

} // end of method distinctrng


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrng 
// DESC   : From buffer, within range given by [ibeg,iend], get indices where there 
//          are distinct values. 
//          Use apariingly, as this is slow, and unoptimized.
// ARG    : 
//          ibeg    : beginning of range (after gindexd_.beg())
//          n       : number of elements to scan
//          is      : stride count for search
//                    found; else = nd;
//          indices : GTVector containing indices into the buffer which contain
//                    the distinct values. Must be dimensioned >= any possible nd,
//                    as no checking is done. This array is deleted and re-allocated.
//          nd      : number of distinct values. nd=0 is not considered an error
//          tunique : tmp array of type T, of size at least of size()
//          itmp    : tmp array of type GSIZET, of size at least of size()
// RETURNS    :  TRUE on success; else FALSE 
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng(GSIZET ibeg, GSIZET n, GSIZET is,
                         GSIZET *&indices, GSIZET  &nd,
                         T * const &tunique, GSIZET * const &itmp)
{
  GEOFLOW_TRACE();
  GLLONG i;
  GLLONG nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    bcont = FALSE;
    for ( auto j=0; j<nfound && !bcont; j++ ) {
      bcont = this->data_[i] == tunique[j];
    }

    if ( !bcont ) { // if not in unique buff, add it:
      itmp   [nfound] = i;
      tunique[nfound] = this->data_[i];
      nfound++;
    }
  }
  if ( nd < nfound ) {
    if ( indices != NULLPTR ) delete [] indices;
    nd = nfound;
    indices = new GSIZET [nd];
  }
  for ( auto i=0; i<nfound; i++ ) {
    indices[i] = itmp[i];
  }

  return nfound;

} // end of method distinctrng


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrng_floor
// DESC   : From buffer, within range given by [ibeg,iend], get indices where there 
//          are distinct values s.t. vals > floor
//          Use apariingly, as this is slow, and unoptimized.
// ARG    : 
//          ibeg    : beginning of range (after gindexd_.beg())
//          n       : number of elements to scan
//          is      : stride count for search
//          vals    : array of values of size nd if nd > number of distinct values
//                    found. Must be deleted by caller
//          indices : array containing indices into the buffer which contain
//                    the distinct values. Must be dimensioned >= any possible nd,
//                    as no checking is done. This array is deleted and re-allocated.
//                    Must be deleted by caller.
//          nd      : number of distinct values. nd=0 is not considered an error
//          floor   : s.t. vals > floor
//          tunique : tmp array of type T, of size at least size()
//          itmp    : tmp array of type GSIZET, of size at least size()
// RETURNS    :  number distinct values found
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng_floor(GSIZET ibeg, GSIZET n, GSIZET is, T *&vals,
                               GSIZET *&indices, GSIZET  &nd, T floor, 
                               T * const &tunique, GSIZET * const &itmp)
{
  GEOFLOW_TRACE();
  GLLONG nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( auto i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    if ( this->data_[i] <= floor ) continue;
    bcont = FALSE;
    for ( auto j=0; j<nfound && !bcont; j++ ) {
      bcont = this->data_[i] == tunique[j];
    }

    if ( !bcont ) { // if not in unique buff, add it:
      itmp   [nfound] = i;
      tunique[nfound] = this->data_[i];
      nfound++;
    }
  }

  if ( nd < nfound ) {
    if ( vals    != NULLPTR ) delete [] vals;   vals = NULLPTR;
    if ( indices != NULLPTR ) delete [] indices; indices = NULLPTR;
    nd = nfound;
    vals    = new T      [nd];
    assert(vals != NULLPTR );
    indices = new GSIZET [nd];
    assert(indices != NULLPTR );
  }
  for ( auto i=0; i<nfound; i++ ) {
    vals   [i] = tunique[i];
    indices[i] = itmp[i];
  }

  return nfound;

} // end of method distinctrng_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrng_floor
// DESC   : From buffer, within range given by [ibeg,iend], get indices where there 
//          are distinct values s.t. vals > floor
//          Use apariingly, as this is slow, and unoptimized.
// ARG    : 
//          ibeg    : beginning of range (after gindexd_.beg())
//          n       : number of elements to scan
//          is      : stride count for search
//                    found. Must be deleted by caller
//          indices : array containing indices into the buffer which contain
//                    the distinct values. Must be dimensioned >= any possible nd,
//                    as no checking is done. This array is deleted and re-allocated.
//                    Must be deleted by caller.
//          nd      : number of distinct values. nd=0 is not considered an error
//          floor   : s.t. vals > floor
//          tunique : tmp array of type T, of size at least of size()
//          itmp    : tmp array of type GSIZET, of size at least of size()
// RETURNS    :  number distinct values found
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng_floor(GSIZET ibeg, GSIZET n, GSIZET is, 
                               GSIZET *&indices, GSIZET  &nd, T floor,
                               T * const &tunique, GSIZET * const &itmp)
{
  GEOFLOW_TRACE();
  GLLONG nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( auto i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+n; i+=is+this->gindex_.stride()-1 ) {
    if ( this->data_[i] <= floor ) continue;
    bcont = FALSE;
    for ( auto j=0; j<nfound && !bcont; j++ ) {
      bcont = this->data_[i] == tunique[j];
    }
    
    if ( !bcont ) { // if not in unique buff, add it:
      itmp   [nfound] = i;
      tunique[nfound] = this->data_[i];
      nfound++;
    }
  }

  if ( nd < nfound ) {
    if ( indices != NULLPTR ) delete [] indices;
    nd = nfound;
    indices = new GSIZET [nd];
  }
  for ( auto i=0; i<nfound; i++ ) {
    indices[i] = itmp[i];
  }

  return nfound;

} // end of method distinctrng_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : distinct 
// DESC   : From buffer, get indices where there are distinct values. This
//          is meant to be used for meaningful arithmetic Ts.
//          Use apariingly, as this is slow, and unoptimized.
// ARG    : indices : GTVector containing indices into the buffer which contain
//                    the distinct values. Must be dimensioned >= any possible nd,
//                    as no checking is done. This array is deleted and re-allocated.
//          nd      : number of distinct values. nd=0 is not considered an error
//          tunique : tmp array of type T, of size at least of size()
//          itmp    : tmp array of type GSIZET, of size at least of size()
// RETURNS    :  no. distinct elements
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinct(GSIZET  *&indices, GSIZET  &nd, 
                      T * const &tunique, GSIZET *const &itmp)
{
  GEOFLOW_TRACE();

  GSIZET n = distinctrng(this->gindex_.beg(), this->gindex_.end()-this->gindex_.end()+1, 1, indices, nd, tunique, itmp);

  return n; 

} // end of method distinct 



//**********************************************************************************
//**********************************************************************************
// METHOD : distinct_floor 
// DESC   : From buffer, get indices where there are distinct values >= floor. This
//          is meant to be used for meaningful arithmetic Ts.
//          Use apariingly, as this is slow, and unoptimized.
// ARGS   : indices : GTVector containing indices into the buffer which contain
//                    the distinct values. Must be dimensioned >= any possible nd,
//                    as no checking is done. This array is deleted and re-allocated.
//          nd      : size of 'indices' array; modified as necessary
//          floor   : min value; s.t. distinct value is >= floor
//          tunique : tmp array of type T, of size at least of size()
//          itmp    : tmp array of type GSIZET, of size at least of size()
// RETURNS :  no. distinct elements
//**********************************************************************************
//#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinct_floor(GSIZET  *&indices, GSIZET  &nd, 
                            T floor, T * const &tunique, GSIZET * const &itmp)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
    "Invalid template type: distinct_floor(GSIZET*, GSIZET&, T)");

  GSIZET n = distinctrng_floor(this->gindex_.beg(), this->gindex_.end(), 1, indices, nd, floor, tunique, itmp);

  return n;
} // end of method distinct_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : sortdecreasing (1)
// DESC   : Sorts buffer in decreasing order, in place. 
// ARGS   : none.
// RETURNS: none
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::sortdecreasing()
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortdecreasing(1)");

  T        tmp;

  if ( gindex_.stride() > 1 ) {
    // Perhaps a more efficient algorithm (e.g., quick sort) 
    // would be better, but for now...
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      for ( auto i=j; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
        if ( this->data_[i] >= this->data_[j] ) {
          tmp      = this->data_[j];
          this->data_[j] = this->data_[i];
          this->data_[i] = tmp;
        }
      }
    }
  }
  else {
    quicksortl2s(this->data_, this->gindex_.beg(), this->gindex_.end());
  }


} // sortdecreasing (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : sortdecreasing (2)
// DESC   : Sorts buffer in decreasing order, in place. 
// ARGS   : isort: integer array that contains indices in increasing order s.t.
//                 if B = *this once sorted, and A, its unsored version, 
//                 then B = A(isort). Resized if necessary.
// RETURNS: none
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::sortdecreasing(GTVector<GSIZET> &isort)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortdecreasing(2)");

  GSIZET   itmp;
  T        tmp;

  isort.resize(this->size());
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    isort[j] = j;
  }

  if ( gindex_.stride() > 1 ) {
    // Perhaps a more efficient algorithm (e.g., quick sort) 
    // could be used, but for now...
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      for ( auto i=j+1; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
        if ( this->data_[i] >= this->data_[j] ) {
          tmp      = this->data_[j];
          this->data_[j] = this->data_[i];
          this->data_[i] = tmp;
          itmp     = isort[j];
          isort[j] = isort[i];
          isort[i] = itmp;
        }
      }
    }
  }
  else {
    quicksortl2s(this->data_, isort.data(), this->gindex_.beg(), this->gindex_.end());
  }

} // sortdecreasing(2)


//**********************************************************************************
//**********************************************************************************
// METHOD : sortincreasing (1)
// DESC   : sorts buffer in increasing order
// ARGS   : none.
// RETURNS: nona
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::sortincreasing()
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortincreasing(1)");

  T        tmp;

  if ( gindex_.stride() > 1 ) {
    // Perhaps a more efficient algorithm (e.g., quick sort) 
    // would be better, but for now...
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      for ( auto i=j+1; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
        if ( this->data_[i] <= this->data_[j] ) {
          tmp      = this->data_[j];
          this->data_[j] = this->data_[i];
          this->data_[i] = tmp;
        }
      }
    }

  }
  else {
    quicksorts2l(this->data_, this->gindex_.beg(), this->gindex_.end());
  }

} // sortincreasing (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : sortincreasing (2)
// DESC   : sorts buffer in increasing order, but does not change
//          the data, returning a vector with indices that point
//          to buffer members in increeaseing order
// ARGS   : isort: integer array that contains indices in increasing order s.t.
//                 if B = *this once sorted, and A, its unsored version, 
//                 then B = A(isort). Resized if necessary.
// RETURNS: nona
//**********************************************************************************
//#pragma acc routine vector
template<class T>
void
GTVector<T>::sortincreasing(GTVector<GSIZET> &isort)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortincreasing (2)");

  GSIZET   itmp;
  T        tmp;

  isort.resize(this->size());
  for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    isort[j] = j;
  }

  if ( gindex_.stride() > 1 ) {
    // Perhaps a more efficient algorithm (e.g., quick sort) 
    // could be used, but for now...
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      for ( auto i=j; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
        if ( this->data_[i] <= this->data_[j] ) {
          tmp      = this->data_[j];
          this->data_[j] = this->data_[i];
          this->data_[i] = tmp;
          itmp      = isort[j];
          isort[j]  = isort[i];
          isort[i]  = itmp;
        }
      }
    }
  }
  else {
    quicksorts2l(this->data_, isort.data(), this->gindex_.beg(), this->gindex_.end());
  }

} // sortincreasing (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : partitions2l (1)
// DESC   : partition function for sort small-to-large values
// ARGS   : a    : vector to be sorred
//          start: starting index   
//          end  : ending index   
// RETURNS: pivot value index
//**********************************************************************************
template<class T>
GLLONG
GTVector<T>::partitions2l(T *a, GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: partitions2l");

        GLLONG i, j;  
  const GLLONG P_index= start + (end - start)/2; // P-index indicates the pivot value index
  const T      pivot=a[P_index]; // pivot value
    
  // Move the mid-point value to the front:
  std::swap(a[P_index], a[start]);

  // Check if array value is less than pivot
  // then place it at left side by swapping 
  i = start + 1;
  j = end;
  while ( i <= j ) {
    while ( i <= j && a[i] <= pivot ) i++;
    while ( i <= j && a[j] >  pivot ) j--;
    if ( i < j ) {   
      std::swap(a    [i], a    [j]);
    }
  } 
  std::swap(a[i-1], a[start]);

  // Return the pivot value index
  return i - 1;

} // end, method partitions2l (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : partitions2l (2)
// DESC   : partition function for sort small-to-large values
// ARGS   : a    : vector to be sorred
//          isort: sort indices in original vector, s.t if B = this vector sorted, 
//                 and A the unsorted version, then B = A(isort); Note: isort 
//                 must be filled on entry.
//          start: starting index   
//          end  : ending index   
// RETURNS: pivot value index
//**********************************************************************************
template<class T>
GLLONG
GTVector<T>::partitions2l(T *a, GSIZET *isort,  GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: partitions2l");

        GLLONG i, j;  
  const GLLONG P_index= start + (end - start)/2; // P-index indicates the pivot value index
  const T      pivot=a[P_index]; // pivot value
    
  // Move the mid-point value to the front:
  std::swap(a    [P_index], a    [start]);
  std::swap(isort[P_index], isort[start]);

  // Check if array value is less than pivot
  // then place it at left side by swapping 
  i = start + 1;
  j = end;
  while ( i <= j ) {
    while ( i <= j && a[i] <= pivot ) i++;
    while ( i <= j && a[j] >  pivot ) j--;
    if ( i < j ) {
      std::swap(a    [i], a    [j]);
      std::swap(isort[i], isort[j]);
    }
  } 
  std::swap(a    [i-1], a    [start]);
  std::swap(isort[i-1], isort[start]);
  
  // Return the pivot value index
  return i - 1;
  
} // end, method partitions2l (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : partitionl2s (1)
// DESC   : partition function for sort large-to-small values
// ARGS   : a    : vector to be sorred
//          start: starting index   
//          end  : ending index   
// RETURNS: pivot value index
//**********************************************************************************
template<class T>
GLLONG
GTVector<T>::partitionl2s(T *a, GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: partitionl2s");

        GLLONG i, j;  
  const GLLONG P_index= start + (end - start)/2; // P-index indicates the pivot value index
  const T      pivot=a[P_index]; // pivot value

  // Move the mid-point value to the front:
  std::swap(a[P_index], a[start]);
    
  // Check if array value is less than pivot
  // then place it at left side by swapping 
  i = start + 1;
  j = end;
  while ( i <= j ) {
    while ( i <= j && a[i] >= pivot ) i++;
    while ( i <= j && a[j] <  pivot ) j--;
    if ( i < j ) {
      std::swap(a    [i], a    [j]);
    }
  } 
  std::swap(a[i-1], a[start]);
  
  // Return the pivot value index
  return i - 1;

} // end, method partitionl2s (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : partitionl2s (2)
// DESC   : partition function for sort large-to-small values
// ARGS   : a    : vector to be sorred
//          isort: sort indices in original vector, s.t if B = this vector sorted, 
//                 and A the unsorted version, then B = A(isort); Note: isort 
//                 must be filled on entry.
//          start: starting index   
//          end  : ending index   
// RETURNS: pivot value index
//**********************************************************************************
template<class T>
GLLONG
GTVector<T>::partitionl2s(T *a, GSIZET *isort, GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: partitionl2s");

        GLLONG i, j;  
  const GLLONG P_index= start + (end - start)/2; // P-index indicates the pivot value index
  const T      pivot=a[P_index]; // pivot value

  // Move the mid-point value to the front:
  std::swap(a    [P_index], a    [start]);
  std::swap(isort[P_index], isort[start]);
    
  // Check if array value is less than pivot
  // then place it at left side by swapping 
  i = start + 1;
  j = end;
  while ( i <= j ) {
    while ( i <= j && a[i] >= pivot ) i++;
    while ( i <= j && a[j] <  pivot ) j--;
    if ( i < j ) {
      std::swap(a    [i], a    [j]);
      std::swap(isort[i], isort[j]);
    }
  } 
  std::swap(a    [i-1], a    [start]);
  std::swap(isort[i-1], isort[start]);
  
  // Return the pivot value index
  return i - 1;

} // end, method partitionl2s (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : quicksorts2l (1)
// DESC   : Do quick sort from small to large values
// ARGS   : a    : vector to be sorred
//          start: starting index   
//          end  : ending index   
// RETURNS: none.
//**********************************************************************************
template<class T>
void
GTVector<T>::quicksorts2l(T *a, GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: quicksorts2l");

  GLLONG P_index;

  if ( start >= end ) return;

  P_index = partitions2l(a,start,end);
  quicksorts2l(a,start,P_index-1);
  quicksorts2l(a,P_index+1,end);

} // end, method quicksorts2l (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : quicksorts2l (2)
// DESC   : Do quick sort from small to large values
// ARGS   : a    : vector to be sorred
//          isort: sort indices in original vector, s.t if B = this vector sorted, 
//                 and A the unsorted version, then B = A(isort). Note: isort 
//                 must be filled on entry.

//          start: starting index   
//          end  : ending index   
// RETURNS: none.
//**********************************************************************************
template<class T>
void
GTVector<T>::quicksorts2l(T *a, GSIZET *isort,  GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: quicksorts2l");

  GLLONG P_index;

  if ( start >= end ) return;

  P_index = partitions2l(a,isort,start,end);
  quicksorts2l(a,isort,start,P_index-1);
  quicksorts2l(a,isort,P_index+1,end);

} // end, method quicksorts2l (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : quicksortl2s (1)
// DESC   : Do quick sort from large to small values
// ARGS   : a    : vector to be sorred
//          start: starting index   
//          end  : ending index   
// RETURNS: none.
//**********************************************************************************
template<class T>
void
GTVector<T>::quicksortl2s(T *a, GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: quicksortl2s");

  GLLONG P_index;

  if ( start >= end ) return;

  P_index = partitionl2s(a,start,end);
  quicksortl2s(a,start,P_index-1);
  quicksortl2s(a,P_index+1,end);

} // end, method quicksortl2s (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : quicksortl2s (2)
// DESC   : Do quick sort from large to small values
// ARGS   : a    : vector to be sorred
//          isort: sort indices in original vector, s.t if B = this vector sorted, 
//                 and A the unsorted version, then B = A(isort). Note: isort 
//                 must be filled on entry.
//          start: starting index   
//          end  : ending index   
// RETURNS: none.
//**********************************************************************************
template<class T>
void
GTVector<T>::quicksortl2s(T *a, GSIZET *isort, GLLONG start, GLLONG end)
{
  GEOFLOW_TRACE();
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
  "Invalid template type: quicksortl2s");

  GLLONG P_index;

  if ( start >= end ) return;

  P_index = partitionl2s(a,isort,start,end);
  quicksortl2s(a,isort,start,P_index-1);
  quicksortl2s(a,isort,P_index+1,end);

} // end, method quicksortl2s (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : add_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::add_impl_(const GTVector<T> &obj, std::false_type d)
{
  GEOFLOW_TRACE();
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(1);
  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
        vret[j] = this->data_[j] + obj[j];
    } 
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
        vret[j] = this->data_[j] + obj[0];
    } 
  }
  

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return vret;
} // end, add_impl_

//**********************************************************************************
//**********************************************************************************
// METHOD : add_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::add_impl_(const GTVector &obj, std::true_type d)
{
  GEOFLOW_TRACE();
  GTVector vret(this->gindex_);
  
  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);

  GSIZET nn = vret.getIndex().end() - vret.getIndex().beg() + 1;
  if      ( std::is_same<T,GFLOAT>::value ) {
    fzaxpby((GFLOAT*)vret.data(), (GFLOAT*)(this->data()), (GFLOAT*)&a,
            (GFLOAT*)(obj.data()), (GFLOAT*)&b, &nn, &icsz_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dzaxpby((GDOUBLE*)vret.data(), (GDOUBLE*)(this->data()), (GDOUBLE*)&a,
            (GDOUBLE*)(obj.data()), (GDOUBLE*)&b, &nn, &icsz_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qzaxpby((GQUAD*)vret.data(), (GQUAD*)(this->data()), (GQUAD*)&a,
            (GQUAD*)(obj.data()), (GQUAD*)&b, &nn, &icsz_);
  }
  else {
    assert(FALSE);
  }


  #if defined(_G_AUTO_UPDATE_DEV)
      vret.updatedev();
  #endif

  return vret;
} // end, add_impl_

//**********************************************************************************
//**********************************************************************************
// METHOD : sub_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::sub_impl_(const GTVector &obj, std::false_type d)
{
  GEOFLOW_TRACE();
  GTVector vret(this->gindex_);

  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      vret[j] = this->data_[j] - obj[j];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      vret[j] = this->data_[j] - obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
      vret.updatedev();
  #endif

  return vret;
} // end, sub_impl_


//**********************************************************************************
//**********************************************************************************
// METHOD : sub_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::sub_impl_(const GTVector &obj, std::true_type d)
{
  GEOFLOW_TRACE();
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  
  GSIZET nn = vret.getIndex().end() - vret.getIndex().beg() + 1;
  if      ( std::is_same<T,GFLOAT>::value ) {
    fzaxpby((GFLOAT*)vret.data(), (GFLOAT*)(this->data()), (GFLOAT*)&a,
            (GFLOAT*)(obj.data()), (GFLOAT*)&b, &nn, &icsz_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dzaxpby((GDOUBLE*)vret.data(), (GDOUBLE*)(this->data()), (GDOUBLE*)&a,
            (GDOUBLE*)(obj.data()), (GDOUBLE*)&b, &nn, &icsz_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qzaxpby((GQUAD*)vret.data(), (GQUAD*)(this->data()), (GQUAD*)&a,
            (GQUAD*)(obj.data()), (GQUAD*)&b, &nn, &icsz_);
  }
  else {
    assert(FALSE);
  }

  #if defined(_G_AUTO_UPDATE_DEV)
      vret.updatedev();
  #endif

  return vret;
} // end, sub_impl_

//**********************************************************************************
//**********************************************************************************
// METHOD : mul_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::mul_impl_(const GTVector &obj, std::false_type d)
{
  GEOFLOW_TRACE();
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  if ( obj.size() > 1 ) {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      vret[j] = this->data_[j] * obj[j];
    }
  }
  else {
    for ( auto j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      vret[j] = this->data_[j] * obj[0];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
      vret.updatedev();
  #endif

  return vret;
} // end, mul_impl_


//**********************************************************************************
//**********************************************************************************
// METHOD : mul_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>
GTVector<T>::mul_impl_(const GTVector &obj, std::true_type d)
{
  GEOFLOW_TRACE();
    return mul_impl_(obj,std::false_type());
} // end, mul_impl_

//**********************************************************************************
//**********************************************************************************
// METHOD : concat
// DESC   : Concatenate input array with existing vector. This vector
//          must have stride 1, and input array is also taken to have 
//          no structure/stride.
// ARGS   : arr : input array
//          narr: no. elements in array
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
void GTVector<T>::concat(T *arr, GSIZET narr)
{
  GEOFLOW_TRACE();
  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  assert(istride       == 1 
      && ipad           > 0 
      && gindex_.beg() == 0
      && gindex_.end() == n_-1
      &&  "Global index structure not allowed in 'concat'");

  this->reserve(this->size()+narr);
  
  for ( auto j=0; j<narr; j++ ) {
    data_[n_+j] = arr[j];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
} // end, concat
