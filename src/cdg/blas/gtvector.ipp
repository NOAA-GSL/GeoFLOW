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
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Basic
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<class T>
GTVector<T>::GTVector():
data_   (NULLPTR),
n_      (0),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
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
data_   (NULLPTR),
n_      (n),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
{
  data_ = new T [n_];
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
data_   (NULLPTR),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
{
  gindex_ = gi;
  gindex_keep_ = gindex_;
  n_=gindex_.end()+1+gindex_.pad();

  data_ = new T [n_];

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
data_   (NULLPTR),
n_      (obj.size()),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
{
  data_ = new T [n_];
  
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
data_   (NULLPTR),
n_      (n/istride),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
{
  data_ = new T [n_];

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
data_   (NULLPTR),
n_      (n/istride),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
{
  if ( bdatalocal_ ) {
    data_ = new T [n_];
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
data_   (NULLPTR),
n_      (obj.size()),
bdatalocal_ (TRUE),
bconstdata_ (FALSE)
{
  data_ = new T [n_];
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
}


//**********************************************************************************
//**********************************************************************************
// METHOD : bconstdata
// DESC   : Sets vector to be a 'constant'. Member data is resized to be of
//          capacity 1, and any get of the data (for any index) returns the
//          value in element 0 of data block
// ARGS   : bdoconst : TRUE or false to make vector 'constant' or not.
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTVector<T>::bconstdata(GBOOL doconst) 
{
  bconstdata_ = doconst;

  if ( bconstdata_ ) {
    if ( data_ != NULLPTR ) delete data_;
    n_ = 1;
    data_ = new T [n_];
    gindex_(n_, n_, 0, n_-1, 1,  0);
  }
} // end, method bconstdata


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
  if ( bconstdata_ ) {
    clear();
    resizem(1);
  }

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
    if ( this->data_ != NULLPTR ) delete [] this->data_;
    this->data_ = new T [nnew];
    n_ = nnew;
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
  if ( bconstdata_ ) {
    clear();
    resize(1);
  }

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

  if ( nnew != n_ ) {
    if ( data_ != NULLPTR ) delete [] data_;
    data_ = new T [nnew];
    n_ = nnew;
  }
  gindex_(nnew, nnew, ibeg, iend, istride, ipad);

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
  if ( bconstdata_ ) {
    clear();
    resize(1);
  }

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
  if ( bconstdata_ ) {
    clear();
    resize(1);
  }

  assert(bdatalocal_ && "Can reserve only on data local vector");

  #if defined(_G_AUTO_CREATE_DEV)
//  #pragma acc exit data delete( data_[0:n_-1] )
    #pragma acc exit data delete( data_[0:n_-1], this[0:1] )
  #endif

  T *ttmp;
  GLLONG ibeg    = gindex_.beg();
  GLLONG iend    = gindex_.end();
  GLLONG istride = gindex_.stride();
  GLLONG ipad    = gindex_.pad();

  // Check: is following exception-safe? No....
  ttmp  = new T [nnew];

  // Copy old data to temp buffer:
  if ( nnew > n_ ) { // growing
    for ( GSIZET j=0; j<n_; j++ ) ttmp[j] = this->data_[j];
    if ( data_ != NULLPTR ) delete [] data_;
    data_ = new T [nnew];

    // Copy only what was there already to expanded buffer:
    for ( GSIZET j=0; j<n_; j++ ) data_[j] = ttmp[j];
    gindex_(nnew, nnew, ibeg, iend, istride, ipad);
    n_ = nnew;
  }
  else if ( nnew < n_ ) { // shrinking
    for ( GSIZET j=0; j<nnew; j++ ) ttmp[j] = this->data_[j];
    if ( this->data_ != NULLPTR ) delete [] this->data_;
    this->data_ = new T [nnew];

    // Copy only the new amount, leaving the reminder
    // 'uninitialized':
    n_ = nnew;
    for ( GSIZET j=0; j<n_; j++ ) this->data_[j] = ttmp[j];
    gindex_(n_, n_, ibeg, MIN(n_-1,iend), istride, ipad);
  }

  delete [] ttmp;

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

  if ( data_ != NULLPTR ) delete [] data_;
  gindex_(0, 0, 0, -1, 0,  0);

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
  if ( bconstdata_ ) return;

  GSIZET nnew = gindex_.end() + gindex_.pad() + 2;

  if ( nnew > n_ ) { // reallocate if required
    reserve(nnew);
  }


  GIndex gi      = gindex_; 
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
// METHOD : data
// DESC   : Get pointer to data block. 
// ARGS   : 
// RETURNS: T* pointer to data
//**********************************************************************************
template<class T> 
T *GTVector<T>::data() 
{
  return data_+gindex_.beg();
} // end, method data 


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
void  GTVector<T>::range(GSIZET ibeg, GSIZET iend) 
{

  if ( bconstdata_ ) return;

  assert(ibeg < n_ && iend < n_ && "Invalid range specification");

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
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector<T>&
//**********************************************************************************
template<class T>
GTVector<T>&
GTVector<T>::operator+(GTVector<T> &obj) 
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
GTVector<T>&
GTVector<T>::operator-(GTVector<T> &obj) 
{
  return this->sub_impl_(obj, typename std::is_floating_point<T>::type());
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
GTVector<T>::operator+=(T b)
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
GTVector<T>::operator-=(T b)
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
GTVector<T>::operator*=(T b)
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
// METHOD : operator * 
// DESC   : product of vector and constant, create new vector
// ARGS   : GTVector &
// RETURNS: void
//**********************************************************************************
template<class T>
GTVector<T>&
GTVector<T>::operator*(T b)
{
  GTVector<T> *vret = new GTVector<T>(gindex_);
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    (*vret)[j] = this->data_[j] * b;
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;

} // end, operator* 


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
GTVector<T>::operator+=(GTVector<T> &obj)
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
GTVector<T>::operator-=(GTVector<T> &b)
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
GTVector<T>::operator*=(GTVector<T> &b)
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
// METHOD : minn
// DESC   : Find max of first n elements
// ARGS   : n : num elements past gindex.beg to check
// RETURNS: T-type max
//**********************************************************************************
template<class T>
T
GTVector<T>::minn(GSIZET n)
{
  T fm = std::numeric_limits<T>::max();

  for ( GLLONG j=this->gindex_.beg(); j<this->gindex_.beg()+n && j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    fm = MIN(fm,this->data_[j]);
  }

  return fm;

} // end minn


//**********************************************************************************
//**********************************************************************************
// METHOD : min
// DESC   : Find max of all members
// ARGS   : none.
// RETURNS: T-type max
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

  GLLONG  m=0;
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

  while ( i < this->gindex_.beg()+n &&  data_[i] != val ) i++;

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

  while ( i < this->gindex_.beg()+n && data_[i] != val ) i++;

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
// RETURNS    :  TRUE on success; else FALSE 
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng(GSIZET ibeg, GSIZET n, GSIZET is, T *&vals,
                                GSIZET *&indices, GSIZET  &nd)
{

  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: distinctrng");

  GLLONG i, j, nfound;
  GBOOL bcont;

  T tunique[this->n_];
  GSIZET itmp[this->n_]; 

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    bcont = FALSE;
    for ( j=0; j<nfound;j++ ) {
      if ( this->data_[i] == tunique[j] ) break;
    }
    if ( j < nfound ) bcont= TRUE; // contained in tunique

    if ( !bcont ) { // if not in unique buff, add it:
      itmp   [nfound] = i;
      tunique[nfound] = this->data_[i];
      nfound++;
    }
  }
  if ( nd < nfound ) {
    if ( vals    != NULLPTR ) delete [] vals;
    if ( indices != NULLPTR ) delete [] indices;
    nd = nfound;
    vals    = new T      [nd];
    indices = new GSIZET [nd];
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
// RETURNS    :  TRUE on success; else FALSE 
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng(GSIZET ibeg, GSIZET n, GSIZET is,
                                GSIZET *&indices, GSIZET  &nd)
{

  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: distinctrng");

  GLLONG i, j, nfound;
  GBOOL bcont;

  GSIZET tunique[this->n_];
  GSIZET itmp[this->n_]; 

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    bcont = FALSE;
    for ( j=0; j<nfound; j++ ) {
      if ( this->data_[i] == tunique[j] ) break;
    }
    if ( j < nfound ) bcont= TRUE; // contained in tunique

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
// RETURNS    :  number distinct values found
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng_floor(GSIZET ibeg, GSIZET n, GSIZET is, T *&vals,
                                      GSIZET *&indices, GSIZET  &nd, T floor)
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
    "Invalid template type: distinctrng_floor");

  GLLONG i, j, nfound;
  GBOOL bcont;

  T tunique[this->n_];
  GSIZET itmp[this->n_]; 

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+ibeg+n && i<=this->gindex_.end(); i+=this->gindex_.stride()+is-1 ) {
    if ( this->data_[i] <= floor ) continue;
    bcont = FALSE;
    for ( j=0; j<nfound; j++ ) {
      if ( this->data_[i] == tunique[j] ) break;
    }
    if ( j < nfound ) bcont= TRUE; // contained in tunique

    if ( !bcont ) { // if not in unique buff, add it:
      itmp   [nfound] = i;
      tunique[nfound] = this->data_[i];
      nfound++;
    }
  }

  if ( nd < nfound ) {
    if ( vals    != NULLPTR ) delete [] vals;
    if ( indices != NULLPTR ) delete [] indices;
    nd = nfound;
    indices = new GSIZET [nd];
    vals    = new T      [nd];
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
// RETURNS    :  number distinct values found
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinctrng_floor(GSIZET ibeg, GSIZET n, GSIZET is, 
                                      GSIZET *&indices, GSIZET  &nd, T floor)
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
    "Invalid template type: distinctrng_floor");

  GLLONG i, j, nfound;
  GBOOL bcont;

  GSIZET tunique[this->n_];
  GSIZET itmp[this->n_]; 

  // This is the brute-force method, and is _slow_
  nfound = 0;
  for ( i=this->gindex_.beg()+ibeg; i<this->gindex_.beg()+n; i+=is+this->gindex_.stride()-1 ) {
    if ( this->data_[i] <= floor ) continue;
    bcont = FALSE;
    for ( j=0; j<nfound;j++ ) {
      if ( this->data_[i] == tunique[j] ) break;
    }
    if ( j < nfound ) bcont= TRUE; // contained in tunique

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
// RETURNS    :  no. distinct elements
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinct(GSIZET  *&indices, GSIZET  &nd)
{

  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value || std::is_pointer<T>::value &&
    "Invalid template type: distinct(GSIZET*, GSIZET&)");

  GSIZET n = distinctrng(this->gindex_.beg(), this->gindex_.end()-this->gindex_.end()+1, 1, indices, nd);

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
//          nd      : number of distinct values. nd=0 is not considered an error
//         floor   : min value; s.t. distinct value is >= floor
// RETURNS :  TRUE on success; else FALSE 
//**********************************************************************************
#pragma acc routine vector
template<class T>
GSIZET
GTVector<T>::distinct_floor(GSIZET  *&indices, GSIZET  &nd, T floor)
{
  assert(std::is_arithmetic<T>::value || std::is_arithmetic<T>::value &&
    "Invalid template type: distinct_floor(GSIZET*, GSIZET&, T)");

  GSIZET n = distinctrng_floor(this->gindex_.beg(), this->gindex_.end()-this->gindex_.end()+1, 1, indices, nd, floor);

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
// ARGS   : none.
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
GTVector<T>&
GTVector<T>::add_impl_(GTVector<T> &obj, std::false_type d) 
{
  GTVector<T> *vret = new GTVector<T> (this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(1);
  for ( GLLONG j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
      (*vret)[j] = this->data_[j] + obj[j];
  } 
  

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;
} // end, add_impl_

//**********************************************************************************
//**********************************************************************************
// METHOD : add_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>&
GTVector<T>::add_impl_(GTVector<T> &obj, std::true_type d) 
{
  GTVector<T> *vret = new GTVector<T> (this->gindex_);
  
  T a = static_cast<T>(1);
  T b = static_cast<T>(1);
  GMTK::add(*vret, *this, obj, a, b);

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;
} // end, add_impl_

//**********************************************************************************
//**********************************************************************************
// METHOD : sub_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>&
GTVector<T>::sub_impl_(GTVector<T> &obj, std::false_type d)
{
  GTVector<T> *vret = new GTVector<T> (this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  GLLONG j;
  for ( j=this->gindex_.beg(); j<=this->gindex_.end(); j+=this->gindex_.stride() ) {
    (*vret)[j] = this->data_[j] - obj[j];
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;
} // end, sub_impl_


//**********************************************************************************
//**********************************************************************************
// METHOD : sub_impl_
// DESC   :
// ARGS   : GTVector &
// RETURNS: GTVector &
//**********************************************************************************
template<class T>
GTVector<T>&
GTVector<T>::sub_impl_(GTVector<T> &obj, std::true_type d)
{
  GTVector<T> *vret = new GTVector<T> (this->gindex_);

  T a = static_cast<T>(1);
  T b = static_cast<T>(-1);
  GMTK::add(*vret, *this, obj, a, b);

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;
} // end, sub_impl_
