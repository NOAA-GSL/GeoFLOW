//==================================================================================
// Module      : gtvector.ipp
// Date        : 1/1/2018 (DLR)
// Description : Basic template vector class, provides
//               access to contiguous memory. This program unit is
//               included after the header in gtvector.hpp, to resolve
//               'normal' template functions (those not requiring SFINAE).
// Copyright   : Copyright 2018. Colorado State University. All rights reserved
// Derived from: none.
//==================================================================================

#include "gmtk.hpp"

#include <cassert>


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Basic
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<class T>
GTVector<T>::GTVector():
data_    (NULLPTR),
n_             (0),
bdatalocal_ (TRUE)
{
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;

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
template<class T>
GTVector<T>::GTVector(GSIZET n):
data_    (NULLPTR),
n_             (n),
bdatalocal_ (TRUE)
{
  data_ = new T [n_];
  assert(this->data_!= NULLPTR );
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
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
data_    (NULLPTR),
bdatalocal_ (TRUE)
{
  gindex_ = gi;
  gindex_keep_ = gindex_;
  n_=gindex_.end()+1+gindex_.pad();

  data_ = new T [n_];
  assert(this->data_!= NULLPTR );

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
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
data_      (NULLPTR),
n_      (obj.size()),
bdatalocal_   (TRUE)
{
  data_ = new T [n_];
  assert(this->data_!= NULLPTR );
  
  for ( GLLONG j=0; j<obj.capacity(); j++ ) {
    this->data_[j] = obj[j];
  }
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
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
data_     (NULLPTR),
n_      (n/istride),
bdatalocal_  (TRUE)
{
  data_ = new T [n_];
  assert(this->data_!= NULLPTR );

  GLLONG k=0;
  for ( GLLONG j=0; j<n_; j++ ) {
    data_[j] = indata[k];
    k += istride;
  }
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
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
data_     (NULLPTR),
n_      (n/istride),
bdatalocal_  (TRUE)
{
  if ( bdatalocal_ ) {
    data_ = new T [n_];
    assert(this->data_!= NULLPTR );
    GLLONG k=0;
    for ( GLLONG j=0; j<n_; j++ ) {
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
template<class T>
GTVector<T>::GTVector(const GTVector<T> &obj):
data_      (NULLPTR),
n_      (obj.size()),
bdatalocal_   (TRUE)
{
  data_ = new T [n_];
  assert(this->data_!= NULLPTR );
  for ( GLLONG j=0; j<obj.capacity(); j++ ) {
    data_[j] = obj[j];
  }
  gindex_ = obj.gindex_;
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
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
  #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
  if ( data_  != NULLPTR  && bdatalocal_ ) delete [] data_;
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
  assert(bdatalocal_ && "Data not local; cannot resize");

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc exit data delete( data_[0:n_-1] )
    #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
  #endif
  
  GLLONG ibeg    = gindex_.beg();
  GLLONG iend    = ibeg + nnew - 1;
  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  if ( (iend-ibeg+1+ipad) > n_ ) {      // must reallocate; change capacity
    if ( this->data_ != NULLPTR ) { 
      delete [] this->data_;
      this->data_ = NULLPTR; 
    }
    this->n_ = iend-ibeg+1+ipad;
    this->data_ = new T [this->n_];
  }
  gindex_(nnew, nnew, ibeg, iend, istride, ipad);
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc enter data create( data_[0:n_-1] )
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
 
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
  assert(bdatalocal_ && "Data not local; cannot resize");

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc exit data delete( data_[0:n_-1] )
    #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
  #endif

  GLLONG nnew    = gi.end()+gi.pad()+1;
  GLLONG ibeg    = gi.beg();
  GLLONG iend    = gi.end();
  GLLONG istride = gi.stride();
  GLLONG ipad    = gi.pad();

  if ( (iend-ibeg+1+ipad) > n_ ) {      // must reallocate; change capacity
    if ( this->data_ != NULLPTR ) { 
      delete [] this->data_; 
      this->data_ = NULLPTR; 
    }
    this->n_ = iend-ibeg+1+ipad;
    this->data_ = new T [this->n_];
  }
  gindex_(nnew, nnew, ibeg, iend, istride, ipad);
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc enter data create( data_[0:n_-1] )
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
 
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
  if ( nnew > n_ ) {

    assert(bdatalocal_ && "Data not local; cannot resize");
    #if defined(_G_AUTO_CREATE_DEV)
//    #pragma acc exit data delete( data_[0:n_-1] )
      #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
    #endif

    resize(nnew);

    #if defined(_G_AUTO_CREATE_DEV)
//    #pragma acc enter data create( data_[0:n_-1] )
      #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
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
  assert(bdatalocal_ && "Can reserve only on data local vector");

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc exit data delete( data_[0:n_-1] )
    #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
  #endif

  T *ttmp=NULLPTR;
  GLLONG ibeg    = gindex_.beg();
  GLLONG iend    = gindex_.end();
  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  // Check: is following exception-safe? No....
  ttmp  = new T [ibeg+nnew+ipad];
  assert(ttmp != NULLPTR );

  // Copy old data to temp buffer:
  if ( nnew > n_ ) { // growing
    for ( GSIZET j=0; j<n_; j++ ) ttmp[j] = this->data_[j];
    if ( this->data_ != NULLPTR ){
      delete [] this->data_;
      this->data_ = NULLPTR; 
    }
    this->data_ = new T [ibeg+nnew+ipad];
    assert(this->data_ != NULLPTR );

    // Copy only what was there already to expanded buffer,
    // leaving remainder 'uninitialized':
    for ( GSIZET j=0; j<n_; j++ ) this->data_[j] = ttmp[j];
    gindex_(nnew, nnew, ibeg, iend, istride, ipad);
    n_ = nnew;
  }
  else if ( nnew < n_ ) { // shrinking
    for ( GSIZET j=0; j<nnew; j++ ) ttmp[j] = this->data_[j];
    if ( this->data_ != NULLPTR ) {
      delete [] this->data_;
      this->data_ = NULLPTR; 
    }
    this->data_ = new T [ibeg+nnew+ipad];
    assert(this->data_ != NULLPTR );

    // Copy only what of the original fills fills new buffer:
    n_ = nnew;
    for ( GSIZET j=0; j<n_; j++ ) this->data_[j] = ttmp[j];
    gindex_(n_, n_, ibeg, MIN(n_-1,iend), istride, ipad);
  }

  if ( ttmp != NULLPTR ) delete [] ttmp;

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc enter data create( data_[0:n_-1] )
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif
 
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

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc exit data delete( data_[0:n_-1] )
    #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
  #endif
  if ( data_ != NULLPTR ) {
    delete [] data_;
    this->data_ = NULLPTR; 
  }
  n_ = 0;
  gindex_(n_, n_, 0, n_-1, 1,  0);
  gindex_keep_ = gindex_;

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc enter data create( data_[0:n_-1] )
    #pragma acc enter data copyin( this[0:1] ) create( data_[0:n_-1] )
  #endif

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
  return data_; //+gindex_.beg();
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
  return data_ ; //gindex_.beg();
} // end, method data (const)


//**********************************************************************************
//**********************************************************************************
// METHOD : assignment operator= GTVector
// DESC   : Equate to another GTVector
// ARGS   : GTVector<T> & right-hand side arg 
// RETURNS: GTVector & 
//**********************************************************************************
template<class T> 
GTVector<T> &GTVector<T>::operator=(const GTVector<T> &obj) 
{
  if ( this == &obj) {
    return *this;
  }
  
  if ( n_ > 0 ) { // If not allocated, allocate; else it's an error
    assert( (this->n_ < obj.capacity() || this->data_ != NULLPTR) && "L-vector has insufficient size");
  }

  if ( data_ == NULLPTR ) {
    n_ = obj.capacity();
    data_ = new T [n_];
    assert(this->data_ != NULLPTR );
    gindex_ = obj.gindex_;
    gindex_keep_ = gindex_;
  }

  for ( GLLONG j=gindex_.beg(); j<=gindex_.end(); j++ ) {
    data_[j] = obj[j-gindex_.beg()];
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
  
  if ( n_ > 0 ) { // If not allocated, allocate; else it's an error
    assert( (this->n_ < obj.capacity() || this->data_ != NULLPTR) && "L-vector has insufficient size");
  }

  if ( data_ == NULLPTR ) {
    n_ = obj.capacity();
    data_ = new T [n_];
    assert(this->data_ != NULLPTR );
    gindex_(n_, n_, 0, n_-1, 1,  0);
    gindex_keep_ = gindex_;
  }

  for ( GLLONG j=gindex_.beg(); j<=gindex_.end(); j++ ) {
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
  for ( GLLONG j=gindex_.beg(); j<=gindex_.end(); j+=gindex_.stride() ) {
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
  for ( GLLONG j=gindex_.beg(); j<=gindex_.end(); j+=gindex_.stride() ) {
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
  #if defined(_G_BOUNDS_CHK) && !defined(_G_USE_OPENACC)
  if ( gindex_.beg()+n > gindex_.end()+1  ) {
    std::cout <<  "GTVector<T>::set: " << "assignment size mismatch: " << std::endl;
while(1){};
    exit(1);
  }
  #endif
  
  GLLONG j, m=0;
  for ( j=gindex_.beg(); j<MIN(gindex_.beg()+n,gindex_.end()+1) && m < n; j+=gindex_.stride() ) {
    data_[j] = a[m++];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
} // end of method set


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
  #pragma acc update self( data_[0:n_-1] )
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
#if defined(_G_AUTO_UPDATE_DEV)
  #pragma acc update device( data_[0:n_-1] )
#endif
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
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: GMatrix<T>::transpose()");

  T       tmp;
  GSIZET  i, j;

  for ( j=0; j<n; j++ ) {
    for ( i=j; i<n; i++ ) {
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
  return this->mul_impl_(obj, typename std::is_floating_point<T>::type());
} // end, operator-

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
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
#pragma acc routine vector
template<class T>
void
GTVector<T>::operator+=(const GTVector &obj)
{
  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    this->data_[j] += obj[j-this->gindex_.beg()];
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
#pragma acc routine vector
template<class T>
void
GTVector<T>::operator-=(const GTVector &b)
{

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    this->data_[j] -= b[j-this->gindex_.beg()];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator-= 


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *= (1)
// DESC   : point-by-point product
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
void
GTVector<T>::operator*=(const GTVector &b)
{
  GLLONG j;
  T *p = b.data();
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    this->data_[j] *= p[j-this->gindex_.beg()];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif

} // end, operator*= (1)


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
  GBOOL bret = TRUE;

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end() && j<=this->gindex_.end() && bret; j+=this->gindex_.stride() ) {
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
  GBOOL bret = TRUE;

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end() && j<=this->gindex_.end() && bret; j+=this->gindex_.stride() ) {
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
  T fm = std::numeric_limits<T>::min();

  for ( GLLONG j=this->gindex_.beg(); j<this->gindex_.beg()+n && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  GSIZET imax;
  T fm = std::numeric_limits<T>::min();

  for ( GLLONG j=this->gindex_.beg(); j<this->gindex_.end() && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  T fm = std::numeric_limits<T>::min();

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MAX(fm,this->data_[j]);
  }
 
  return fm;

} // end max


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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T fm = std::numeric_limits<T>::min();

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MAX(fm,std::fabs(this->data_[j]));
  }

  return fm;

} // end amax


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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T diff;
  T fm = std::numeric_limits<T>::min();

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end()-1; j+=this->gindex_.stride() ) {
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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T fm = std::numeric_limits<T>::max();

  for ( GLLONG j=this->gindex_.beg(); j<this->gindex_.beg()+n && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  GSIZET imin;
  T fm = std::numeric_limits<T>::max();

  for ( GLLONG j=this->gindex_.beg(); j<this->gindex_.end() && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  T fm = std::numeric_limits<T>::max();

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MIN(fm,this->data_[j]);
  }

  return fm;

} // end min


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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T fm = std::numeric_limits<T>::max();

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  T diff;
  T fm = std::numeric_limits<T>::max();

  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end()-1; j+=this->gindex_.stride() ) {
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
template<class T>
void
GTVector<T>::pointProd(const GTVector<T> &obj, GTVector<T> &ret ) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( obj.size() < this->size() || ret.size() < this->size() ) {
    std::cout << "pointProd(1): " << "incompatible size" << std::endl;
while(1){};
    exit(1);
  }
  #endif

  GLLONG j; 
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    ret[j] = this->data_[j-gindex_.beg()] * obj[j-gindex_.beg()];
  }

} // end, pointProd (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : pointProd (2)
// DESC   : point-by-point multiplication, returned in this
// ARGS   : obj: const GTVector<>  factor
// RETURNS: none
//**********************************************************************************
template<class T>
void
GTVector<T>::pointProd(const GTVector<T> &obj)
{
  #if defined(_G_BOUNDS_CHK)
  if ( obj.size() < this->size() ) {
    std::cout << "pointProd(2): " << "incompatible size" << std::endl;
while(1){};
    exit(1);
  }
  #endif

  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    data_[j] *= obj[j-gindex_.beg()];
  }

} // end, pointProd (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : constProd
// DESC   : multiply this by constant, return
// ARGS   : b  : T-type constant 
//          ret: GTVector & ret vector
// RETURNS: GTVector & 
//**********************************************************************************
template<class T>
void
GTVector<T>::constProd(const T b, GTVector<T> &ret) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( ret.size() < this->size() ) {
    std::cout << "constProd: " << "incompatible size" << std::endl;
while(1){};
    exit(1);
  }
  #endif

  GLLONG j;
  T *dret=ret.data();

  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    dret[j] = this->data_[j] * b;
  }

} // end, pointProd

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
  T      sum=static_cast<T>(0);
  GLLONG j;

  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  GDOUBLE xnorm;
  GLLONG j;

  for ( j=this->gindex_.beg(), xnorm=0; j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    xnorm = MAX(xnorm,fabs(this->data_[j]));
  }
  
  return static_cast<T>(xnorm);
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
  GDOUBLE xnorm;
  GLLONG j;

  for ( j=this->gindex_.beg(), xnorm=0; j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    xnorm += this->data_[j]*this->data_[j];
  }
  
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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");

  GLLONG  j;
  GDOUBLE b;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
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
  assert(std::is_arithmetic<T>::value && "Requires arithmetic template parameter");
  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
//  tmp = sqrt( std::pow<GDOUBLE>(static_cast<GDOUBLE>(data_[j]),2) );
    data_[j] = std::fabs(data_[j]);
  }

} // end, abs


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity
// DESC   : Gets multiplicity of specified value
// ARGS   : val : type T
// RETURNS: GSIZT multiplicity
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity(T val)
{
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
//  "Invalid template type: multiplicity(T)");

  GSIZET mult=0;
  for ( GLLONG i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] == val ) {
      mult++;
    }
  }
  
  return mult;

} // multiplicity


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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity(T val, GSIZET *&index, GSIZET &n)
{
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
//  "Invalid template type: multiplicity(T,GSIZET*,GSIZET&)");

  GLLONG m=0;
  GSIZET mult=0;

  for ( GLLONG i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] == val ) {
      mult++;
    }
  }

  if ( n < mult ) {
    if ( index != NULLPTR ) delete [] index;
    n = mult;
    index = new GSIZET [n];
  }
  for ( GLLONG i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_floor(T val, T floor)
{
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
//  "Invalid template type: multiplicity_floor(T,T)");

  GSIZET mult=0;

  for ( GLLONG i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_floor(T val, GSIZET *&index, GSIZET &n, T floor)
{
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
//  "Invalid template type: multiplicity_floor(T,GSIZET *,GSIZE&,T)");

  GLLONG  m=0;
  GSIZET mult=0;

  for ( GLLONG i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] > floor && this->data_[i] == val ) {
      mult++;
    }
  }

  if ( n < mult ) {
    if ( index != NULLPTR ) delete [] index;
    n = mult;
    index = new GSIZET [n];
  }
  for ( GLLONG i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::multiplicity_ceil(T val, T ceil)
{
//assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
//  "Invalid template type: multiplicity_ceil(T,T)");

  GSIZET index, mult=1;

  if ( !this->contains(val,index) ) return 0 ;

  for ( GSIZET i=this->gindex_.beg(); i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
    if ( this->data_[i] < ceil && this->data_[i] == val ) {
      mult++;
    }
  }
  
  return mult;

} // multiplicity_ceil


//**********************************************************************************
//**********************************************************************************
// METHOD : contains (1)
// DESC   : Determines if candidate value is in the buffer
// ARGS   : val : member to search for in buffer
//          iwhere : first index where val is found first. Valid only
//                   if return is TRUE
// RETURNS: TRUE if member is in list, else FALSE. 
//**********************************************************************************
#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains(T val, GSIZET &iwhere)
{

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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::contains(T val, GSIZET *&iwhere, GSIZET &nw)
{
#if 0
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: contains(T, GSIZET&)");
#endif

  if ( this->data_ == NULLPTR ) return FALSE;

  GLLONG i, n;
  
  for ( i=this->gindex_.beg(), n=0; i<=this->gindex_.end(); i++ ) {
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

  for ( i=this->gindex_.beg(), n=0; i<=this->gindex_.end(); i++ ) {
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
#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::containsn(T val, GSIZET n, GSIZET &iwhere)
{
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
#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains(T val)
{
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
#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::containsn(T val, GSIZET n)
{
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
#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains_floor(T val, GSIZET  &iwhere, T floor, GSIZET istart)
{
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value &&
    "Invalid template type: contains_floor");

  if ( data_ == NULLPTR ) return FALSE;

  GLLONG i;
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
#pragma acc routine vector
template<class T>
GBOOL
GTVector<T>::contains_ceil(T val, GSIZET  &iwhere, T ceil, GSIZET istart)
{
  assert(std::is_arithmetic<T>::value || std::is_enum<T>::value &&
    "Invalid template type: contains_ceil");

  if ( this->data_ == NULLPTR ) return FALSE;

  GLLONG i;
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng(GSIZET ibeg, GSIZET n, GSIZET is, T *&vals,
                         GSIZET *&indices, GSIZET  &nd, T * const &tunique, GSIZET * const &itmp)
{

  GLLONG i, j, nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    bcont = FALSE;
    for ( j=0; j<nfound && !bcont;j++ ) {
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
  for ( i=0; i<nfound; i++ ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng(GSIZET ibeg, GSIZET n, GSIZET is,
                         GSIZET *&indices, GSIZET  &nd,
                         T * const &tunique, GSIZET * const &itmp)
{

  GLLONG i, j, nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    bcont = FALSE;
    for ( j=0; j<nfound && !bcont; j++ ) {
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
  for ( i=0; i<nfound; i++ ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng_floor(GSIZET ibeg, GSIZET n, GSIZET is, T *&vals,
                               GSIZET *&indices, GSIZET  &nd, T floor, 
                               T * const &tunique, GSIZET * const &itmp)
{
  GLLONG i, j, nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    if ( this->data_[i] <= floor ) continue;
    bcont = FALSE;
    for ( j=0; j<nfound && !bcont; j++ ) {
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
  for ( i=0; i<nfound; i++ ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng_floor(GSIZET ibeg, GSIZET n, GSIZET is, 
                               GSIZET *&indices, GSIZET  &nd, T floor,
                               T * const &tunique, GSIZET * const &itmp)
{

  GLLONG i, j, nfound;
  GBOOL bcont;

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+n; i+=is+this->gindex_.stride()-1 ) {
    if ( this->data_[i] <= floor ) continue;
    bcont = FALSE;
    for ( j=0; j<nfound && !bcont; j++ ) {
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
  for ( i=0; i<nfound; i++ ) {
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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinct(GSIZET  *&indices, GSIZET  &nd, 
                      T * const &tunique, GSIZET *const &itmp)
{

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
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinct_floor(GSIZET  *&indices, GSIZET  &nd, 
                            T floor, T * const &tunique, GSIZET * const &itmp)
{
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
#pragma acc routine vector
template<class T>
void
GTVector<T>::sortdecreasing()
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortdecreasing(1)");

  GLLONG   i, j;
  T       tmp;

  // Perhaps a more efficient algorithm (e.g., heapsort) 
  // would be better, but for now...
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    for ( i=j; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
      if ( this->data_[i] > this->data_[j] ) {
        tmp      = this->data_[j];
        this->data_[j] = this->data_[i];
        this->data_[i] = tmp;
      }
    }
  }

} // sortdecreasing (1)

//**********************************************************************************
//**********************************************************************************
// METHOD : sortdecreasing (2)
// DESC   : Sorts buffer in decreasing order, in place. 
// ARGS   : isort: integer array that contains indices in decreasing order.
//                 resized if necessary. 
// RETURNS: none
//**********************************************************************************
#pragma acc routine vector
template<class T>
void
GTVector<T>::sortdecreasing(GTVector<GSIZET> &isort)
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortdecreasing(2)");

  GLLONG   i, j;
  GSIZET   tmp;

  isort.resize(this->size());
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    isort[j] = j;
  }

  // Perhaps a more efficient algorithm (e.g., heapsort) 
  // would be better, but for now...
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    for ( i=j; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
      if ( this->data_[isort[i]] > this->data_[isort[j]] ) {
        tmp      = isort[j];
        isort[j] = isort[i];
        isort[i] = tmp;
      }
    }
  }

} // sortdecreasing(2)


//**********************************************************************************
//**********************************************************************************
// METHOD : sortincreasing (1)
// DESC   : sorts buffer in increasing order
// ARGS   : none.
// RETURNS: nona
//**********************************************************************************
#pragma acc routine vector
template<class T>
void
GTVector<T>::sortincreasing()
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortincreasing(1)");

  GLLONG   i, j;
  T       tmp;

  // Perhaps a more efficient algorithm (e.g., heapsort) 
  // would be better, but for now...
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    for ( i=j; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
      if ( this->data_[i] < this->data_[j] ) {
        tmp      = this->data_[j];
        this->data_[j] = this->data_[i];
        this->data_[i] = tmp;
      }
    }
  }

} // sortincreasing (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : sortincreasing (2)
// DESC   : sorts buffer in increasing order, but does not change
//          the data, returning a vector with indices that point
//          to buffer members in increeaseing order
// ARGS   : isort: integer array that contains indices in increasing order.
//                 resized if necessary. 
// RETURNS: nona
//**********************************************************************************
#pragma acc routine vector
template<class T>
void
GTVector<T>::sortincreasing(GTVector<GSIZET> &isort)
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: sortincreasing (2)");

  GLLONG   i, j;
  GSIZET   tmp;

  isort.resize(this->size());
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    isort[j] = j;
  }

  // Perhaps a more efficient algorithm (e.g., heapsort) 
  // would be better, but for now...
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    for ( i=j; i<=this->gindex_.end(); i+=this->gindex_.stride() ) {
      if ( this->data_[isort[i]] < this->data_[isort[j]] ) {
        tmp       = isort[j];
        isort[j]  = isort[i];
        isort[i]  = tmp;
      }
    }
  }

} // sortincreasing (2)


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
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(1);
  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      vret[j] = this->data_[j] + obj[j];
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
  GTVector vret(this->gindex_);
  
  T a = static_cast<T>(1);
  T b = static_cast<T>(1);
  GMTK::add(vret, *this, obj, a, b);

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
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    vret[j] = this->data_[j] - obj[j];
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
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  GMTK::add(vret, *this, obj, a, b);

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
  GTVector vret(this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    vret[j] = this->data_[j] * obj[j];
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
GTVector<T>::mul_impl_(const GTVector &obj, std::true_type d)
{
    return mul_impl_(obj,std::false_type());
} // end, sub_impl_

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

  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  assert(istride       == 1 
      && ipad           > 0 
      && gindex_.beg() == 0
      && gindex_.end() == n_-1
      &&  "Global index structure not allowed in 'concat'");

  this->reserve(this->size()+narr);
  
  for ( GSIZET j=0; j<narr; j++ ) {
    data_[n_+j] = arr[j];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  updatedev();
  #endif
} // end, concat
