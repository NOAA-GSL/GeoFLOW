//==================================================================================
// Module       : gtmatrix.ipp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template matrix object, whose data is composed of
//                a regular array of contiguous data, ordered like a Fortran
//                matrix in row-major order (row data changing fastest over column data).
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived from : none.
//==================================================================================
#include <limits>
#include "cff_blas.h"
#include "gtmatrix.hpp"
#include "gmtk.hpp"

#define GTMATRIX_ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
	a(k,l)=h+s*(g-h*tau);

#if !defined(_G_MAT_CACHE_SIZE)
  # define _G_MAT_CACHE_SIZE 16
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Basic
// ARGS   : none
// RETURNS: GTMatrix<T>
//**********************************************************************************
template<class T> 
GTMatrix<T>::GTMatrix()
:
n1_                   (0),
n2_                   (0),
icsz_                 (_G_MAT_CACHE_SIZE),
singzero_             (1e-12)
{
  data_.resize(1);

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) 
  #endif

} // end of constructor 1 


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Construct with matrix dimanions (rows, cols)
// ARGS   : size1: no. rows
//          size2: no. cols
// RETURNS: GTMatrix<T>
//**********************************************************************************
template<class T>  
GTMatrix<T>::GTMatrix(const GSIZET   size1, const GSIZET   size2)
:
n1_                   (size1),
n2_                   (size2),
icsz_                 (_G_MAT_CACHE_SIZE),
singzero_             (1e-12)
{

  data_.resize(n1_*n2_);

  zero();

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) 
  #endif

} // end of constructor 2


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (3)
// DESC   : Square matrix constructor
// ARGS   : size1: no. rows, cols
// RETURNS: GTMatrix<T>
//**********************************************************************************
template<class T>  
GTMatrix<T>::GTMatrix(const GSIZET size1)
:
n1_                   (size1),
n2_                   (size1),
icsz_                 (_G_MAT_CACHE_SIZE),
singzero_             (1e-12)
{

  data_.resize(n1_*n2_);

  zero();

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) 
  #endif

} // end of constructor 3


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (4)
// DESC   : Instantiate with C-style data block
// ARGS   : array: data block
//          m1   : no. rows
//          m2   : no. cols
// RETURNS: GTMatrix<T>
//**********************************************************************************
template<class T>  
GTMatrix<T>::GTMatrix(T *array, GSIZET   m1, GSIZET   m2)
:
n1_                   (m1),
n2_                   (m2),
icsz_                 (_G_MAT_CACHE_SIZE),
singzero_             (1e-12)
{

  // build matrix data_ structure:
  data_.resize(n1_*n2_);
  for ( GSIZET j=0; j< n1_*n2_; j++ ) data_[j] = array[j];

  zero();

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] )
  #endif

} // end of constructor 4


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (5)
// DESC   : Instantiate from specied columns of input matrix
// ARGS   : array: data block
//          m    : input matrix to copy from
//          ind  : column indices of m to copy
//          nn   : no. cols in matrix
// RETURNS: GTMatrix<T>
//**********************************************************************************
template<class T>  
GTMatrix<T>::GTMatrix(GTMatrix<T> &m, GSIZET *ind, GSIZET nn)
:
n1_                   (0),
n2_                   (0),
icsz_                 (_G_MAT_CACHE_SIZE),
singzero_             (1e-12)
{

  // Build matrix data structure:
  if ( ind == NULL ) {
    std::cout << "GTMatrix<T>::GTMatrix (4): NULL reference index set" << std::endl;
    exit(1);
  }
  n1_ = m.dim(1);
  n2_ = nn;
  data_.resize(n1_*n2_);
  for ( GSIZET j=0; j<n2_; j++ ) {
    for ( GSIZET i=0; i<n1_; i++ ) {
      data_[i+j*n1_] = m(i,ind[j]);
    }
  }

  zero();

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] )
  #endif

} // end of constructor 5


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method 
// DESC   : 
// ARGS   : 
// RETURNS: GTMatrix<T>
//**********************************************************************************
template<class T> 
GTMatrix<T>::GTMatrix(const GTMatrix<T> &m)
{

  // copy member data:
  n1_       = m.n1_;
  n2_       = m.n2_;
  icsz_     = m.icsz_;
  singzero_ = m.singzero_;
 
  data_.resize(n1_*n2_);

  data_ = m.data_;

  #if defined(_G_AUTO_CREATE_DEV)
    #pragma acc enter data copyin( this[0:1] ) )
  #endif
  updatedev();

} // end of copy constructor


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<class T> 
GTMatrix<T>::~GTMatrix<T>()
{
  #pragma acc exit data delete( data_[0:1], this[0:1] )
}

//**********************************************************************************
//**********************************************************************************
// METHOD : operator = (1)
// DESC   : 
// ARGS   : From existing matrix
// RETURNS: GTMatrix & this
//**********************************************************************************
template<class T> 
GTMatrix<T> &GTMatrix<T>::operator=(const GTMatrix<T> &m)
{
  if ( &m != this ) {
    if ( m.n1_ != n1_ || m.n2_ != n2_ ) {
      std::cout << "GTMatrix<T>::=: incompatible matrices" << std::endl;
      while(1);
      exit(1);
    }
    // copy member data:
    n1_       = m.n1_;
    n2_       = m.n2_;
    singzero_ = m.singzero_;
    data_     = m.data_;
  }

  return *this;

} // end = operator


//**********************************************************************************
//**********************************************************************************
// METHOD : size
// DESC   : 
// ARGS   : GINT direction
// RETURNS: none
//**********************************************************************************
template<class T> 
GSIZET GTMatrix<T>::size(GINT idir) const
{
  if      ( idir == 1 ) return n1_;
  else if ( idir == 2 ) return n2_;
  else                  {
    GError();
    exit(1);
  }

} // end = operator


//**********************************************************************************
//**********************************************************************************
// METHOD : operator = (2)
// DESC   : 
// ARGS   : from scalar T
// RETURNS: none
//**********************************************************************************
template<class T> 
void  GTMatrix<T>::operator=(T m)
{
  GSIZET   i;

  for ( i=0; i<n1_*n2_; i++ ) { 
      data_[i] = m; 
  }

} // end = operator


//**********************************************************************************
//**********************************************************************************
// METHOD : operator *
// DESC   : multiplies this by constant, and returns
//          result, without destroying *this data
// ARGS   : T scalar
// RETURNS: product matrix
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::operator*(T a) 
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator*(T)");

  GTMatrix<T>   *mret = new GTMatrix<T>(n1_,n2_);

  GTVector<T> &vdata = mret->data();
  for ( GSIZET j=0; j<n1_*n2_; j++ ) vdata[j] = a*data_[j];

  return *mret;

} // end of * operator 

//**********************************************************************************
//**********************************************************************************
// METHOD : operator *=
// DESC   : multiplies this by constant, and returns
//          result, modifying *this data
// ARGS   :
// RETURNS: product matrix
//**********************************************************************************
template<class T>
void GTMatrix<T>::operator*=(T a) 
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator*=(T)");

  data_ *= a;

} // end of *= operator 


//**********************************************************************************
//**********************************************************************************
// METHOD : operator +=
// DESC   : Matrix addition: this += a, modifying member data
// ARGS   : GTMatrix &
// RETURNS: none
//**********************************************************************************
template<class T>
void GTMatrix<T>::operator+=(GTMatrix<T> &a) 
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator+=(GTMatrix<T> &)");

  #if defined(_G_BOUNDS_CHK)
  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    std::cout << "GTMatrix<T>::+: incompatible matrices"<< std::endl;
    exit(1);
  }
  #endif

  data_ += a.data();

}


//**********************************************************************************
//**********************************************************************************
// METHOD : operator +=
// DESC   : Scalar self-addition: this += a, modifying member data
// ARGS   : T scalar
// RETURNS: none
//**********************************************************************************
template<class T>
void GTMatrix<T>::operator+=(T a) 
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator+=(T)");

  data_ += a;

}


//**********************************************************************************
//**********************************************************************************
// METHOD : operator -=
// DESC   : Matrix subtraction : this -= a, modifying member data
// ARGS   : GTMatrix &
// RETURNS: GTMatrix 
//**********************************************************************************
template<class T>
void GTMatrix<T>::operator-=(GTMatrix<T> &a) 
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator-=(GTMatrix<T> &)");

  #if defined(_G_BOUNDS_CHK)
  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    std::cout << "GTMatrix<T>::+: incompatible matrices"<< std::endl;
    exit(1);
  }
  #endif

  data_ -= a.data();

}

//**********************************************************************************
//**********************************************************************************
// METHOD : operator +
// DESC   : Matrix addition: this + a
// ARGS   : GTMatrix &
// RETURNS: GTMatrix 
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::operator+(GTMatrix<T> &a)
{

  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator+(GTMatrix<T> &)");

  #if defined(_G_BOUNDS_CHK)
  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    std::cout << "GTMatrix<T>::+: incompatible matrices"<< std::endl;
    exit(1);
  }
  #endif

  GTMatrix<T>  *mret=new GTMatrix<T>(n1_,n2_);

  mret->data() = this->data() + a.data();

  return *mret;

} // operator+

//**********************************************************************************
//**********************************************************************************
// METHOD : operator -
// DESC   : Matrix subtraction: this - a
// ARGS   : GTMatrix &
// RETURNS: GTMatrix  &
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::operator-(GTMatrix<T> &a)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::operator-(GTMatrix<T> &)");

  
  #if defined(_G_BOUNDS_CHK)
  if ( this->n1_ != a.dim(1) || this->n2_ !=a.dim(2) ) {
    std::cout << "GTMatrix<T>::-: incompatible matrices"<< std::endl;
    exit(1);
  } 
  #endif

  GTMatrix<T>  *mret=new GTMatrix<T>(n1_,n2_);

  mret->data() = this->data() - a.data();

  return *mret;

} // operator-


//**********************************************************************************
//**********************************************************************************
// METHOD : resize
// DESC   : resizes dynamically allocated quantities
//          if required
// ARGS   :
// RETURNS:  TRUE on success, else FALSE
//**********************************************************************************
template<class T> 
GBOOL GTMatrix<T>::resize(GSIZET new1, GSIZET new2)
{
  if ( n1_*n2_ == new1*new2 ) return TRUE;

  n1_ = new1;
  n2_ = new2;

  data_.resize(n1_*n2_);

  updatedev(); // Update data on device if necessary

  return TRUE;

} // end, method resize


//**********************************************************************************
//**********************************************************************************
// METHOD : resizem
// DESC   : resizes dynamically allocated quantities
//          if required, with new size an upper limit
// ARGS   :
// RETURNS:  TRUE on success, else FALSE
//**********************************************************************************
template<class T> 
GBOOL GTMatrix<T>::resizem(GSIZET new1, GSIZET new2)
{
  

  if ( new1*new2 > n1_*n2_ ) {
    data_.resizem(new1*new2);
  }

  n1_ = new1;
  n2_ = new2;


  updatedev(); // Update data on device if necessary

  return TRUE;

} // end, method resizeM


//**********************************************************************************
//**********************************************************************************
// METHOD : dim
// DESC   : Array dimension (usable)
//          in direction idir 
// ARGS   :
// RETURNS: GSIZET   size
//**********************************************************************************
template<class T> 
GSIZET GTMatrix<T>::dim(GINT idir) const
{
  if      ( idir == 1 ) return n1_;
  else if ( idir == 2 ) return n2_;
  else    {
    GError();
    exit(1);
  }
} // end, method dim


//**********************************************************************************
//**********************************************************************************
// METHOD : zero
// DESC   : Zeros out data elemnts
// ARGS   :
// RETURNS: none
//**********************************************************************************
template<class T> 
void GTMatrix<T>::zero()
{ 
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::zero");

  data_ = 0;
}  // end, method zero



//**********************************************************************************
//**********************************************************************************
// METHOD : transpose (1)
// DESC   : Computes transpose of *this, but
//          does not destroy data.
// ARGS   :
// RETURNS: transpose of this
//**********************************************************************************
template<class T> 
GBOOL  GTMatrix<T>::transpose(GTMatrix<T> &trans)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::transpose(GTMatrix<T> &)");

  GSIZET  i, j;

  if ( trans.dim(2) !=  n1_ || trans.dim(1) != n2_ ) {
    std::cout << "GTMatrix<T>::transpose(1): incompatible matrix"<< std::endl;
    exit(1);
  }

  for ( j=0; j<n1_; j++ ) {
    for ( i=0; i<n2_; i++ ) {
       trans(i,j) = (*this)(j,i);
    }
  }
  return TRUE;
 
} // end, method transpose (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : transpose (2)
// DESC   : Computes transpose of *this, but
//          does not destroy data. Computes in box of size nx x ny
// ARGS   :
// RETURNS: transpose of this
//**********************************************************************************
template<class T> 
GBOOL  GTMatrix<T>::transpose(T *&trans, GSIZET nx, GSIZET ny)
{

  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::transpose(GTMatrix<T> &,GSIZET,GSIZET)");

  GSIZET  i, j;

  if ( n1_ !=  ny || n2_ != nx ) {
    std::cout << "GTMatrix<T>::transpose(2): incompatible matrix"<< std::endl;
    exit(1);
  }

  for ( j=0; j<nx; j++ ) {
    for ( i=0; i<ny; i++ ) {
       trans[j+i*ny] = (*this)(j,i);
    }
  }
  return TRUE;
 
} // end, method transpose (2)


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : transpose (3)
// DESC   : Computes in-place transpose of this, 
//          so member data is altered.
// ARGS   :
// RETURNS: transpose of this
//**********************************************************************************
template<class T> void GTMatrix<T>::transpose()
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::transpose()");

  T       tmp;
  GSIZET  i, j;

  for ( j=0; j<n2_; j++ ) {
    for ( i=j; i<n1_; i++ ) {
       tmp = (*this)(i+j*n1_);
       (*this)(i+j*n1_) = (*this)(j+i*n2_);
       (*this)(j+i*n2_) = tmp;
    }
  }
  j = n2_;
  n2_ = n1_;
  n1_ = j;

} // end, method transpose (3)
#endif

//**********************************************************************************
//**********************************************************************************
// METHOD : inverse (1)
// DESC   : Computes inverse of this, copying the
//          result to mret
// ARGS   : mret: return matrix for inverse
// RETURNS: inverse of this
//**********************************************************************************
template<class T>
GBOOL GTMatrix<T>::inverse(GTMatrix<T> &mret)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::inverse(GTMatrix<T> &)");

  GString    serr = "GTMatrix<T>::inverse(1): ";
  GLLONG     i, j;
  GBOOL      bRet=TRUE;
  

  if ( mret.dim(1) !=  n1_ || mret.dim(2) != n2_ ) {
    std::cout << serr << "incompatible matrix"<< std::endl;
    exit(1);
  }
  if ( n1_ != n2_ ) {
    std::cout << serr << "matrix not square"<< std::endl;
    exit(1);
  }

 
  GLLONG *indx;
  T      **A0, **A, *col,  d, *b;
  mret = 0.0;

  A0    = new T * [n1_];
  A     = new T * [n1_];
  col   = new T   [n1_];
  b     = new T   [n1_];
  indx  = new GLLONG [n2_];
  for ( i=0; i<n1_; i++ ) {
     A0[i] = new T [n2_];
     A [i] = new T [n2_];
  }
  for ( j=0; j<n2_; j++ ) {
    for ( i=0; i<n1_; i++ ) {
      A0[i][j] = data_[i+j*n1_];
      A [i][j] = A0[i][j];
    }
  }

  if ( !wludcmp(A, n1_, indx, &d) ) {
     std::cout << serr << "wludcmp failed" << std::endl;
     bRet = FALSE;
  }

  for ( j=0; j<n2_ && bRet; j++ ) {
    for ( i=0; i<n1_; i++ ) {  col[i] = 0.0; b[i] = 0.0; }
    col[j] = 1.0; b[j] = col[j];
    if ( !(bRet=lubksb(A, n1_, indx, col)) ) {
       std::cout << serr << "lubjsb failed" << std::endl;
       bRet = FALSE;
       break;
    }
    if ( !(bRet=improve(A0, A, n1_, indx, b, col)) ) {
       std::cout << serr << "improve failed" << std::endl;
       bRet = FALSE;
       break;
    }
    for ( i=0; i<n1_; i++ ) mret(i,j) = col[i];
  }

  for ( i=0; i<n1_; i++ ) {
    delete [] A0[i];
    delete [] A [i];
  }
  delete [] A0;
  delete [] A;
  delete [] col;
  delete [] b;
  delete [] indx;

  return bRet;

} // end, method inverse (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : inverse (2)
// DESC   : Computes inverse of this, copying the result to mret;
//          inverse is made on a box of size nx x ny.
// ARGS   : mret : return matrix
//          nx,ny: num rows, cols to be used in mret
// RETURNS: inverse of this
//**********************************************************************************
template<class T>
GBOOL  GTMatrix<T>::inverse(GTMatrix<T> &mret, GSIZET nx, GSIZET ny)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::inverse(GTMatrix<T> &, GSIZET, GSIZET)");

  GLLONG      i, j, *indx;
  GBOOL       bRet=TRUE;
  T         **A0, **A, *col,  d, *b;
  GString     serr = "GTMatrix<T>::inverse(2): ";
  

  if ( n1_ < nx || n2_ < ny ) {
    std::cout << serr << "incompatible matrix"<< std::endl;
    exit(1);
  }
  if ( nx != ny ) {
    std::cout << serr << "matrix not square"<< std::endl;
    exit(1);
  }

  mret       = 0.0;

  A0    = new T * [nx];
  A     = new T * [nx];
  col   = new T   [nx];
  b     = new T   [nx];
  indx  = new GLLONG [ny];
  for ( i=0; i<nx; i++ ) {
     A0[i] = new T [ny];
     A [i] = new T [ny];
  }
  for ( j=0; j<ny; j++ ) {
    for ( i=0; i<nx; i++ ) {
      A0[i][j] = data_[i+j*nx];
      A [i][j] = A0[i][j];
    }
  }

  if ( !wludcmp(A, nx, indx, &d) ) {
     std::cout << serr << "wludcmp failed" << std::endl;
     bRet = FALSE;
  }

  for ( j=0; j<ny && bRet; j++ ) {
    for ( i=0; i<nx; i++ ) {  col[i] = 0.0; b[i] = 0.0; }
    col[j] = 1.0; b[j] = col[j];
    if ( !(bRet=lubksb(A, nx, indx, col)) ) {
       std::cout << serr << "lubjsb failed" << std::endl;
       bRet = FALSE;
       break;
    }
    if ( !(bRet=improve(A0, A, nx, indx, b, col)) ) {
       std::cout << serr << "improve failed" << std::endl;
       bRet = FALSE;
       break;
    }
    for ( i=0; i<nx; i++ ) mret(i,j) = col[i];
  }

  for ( i=0; i<nx; i++ ) {
    delete [] A0[i];
    delete [] A [i];
  }
  delete [] A0;
  delete [] A;
  delete [] col;
  delete [] b;
  delete [] indx;

  return bRet;

} // end, method inverse (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : inverse (3)
// DESC   : Computes inverse of this, but
//          does not destroy data. A copy is made and
//          returned.
// ARGS   :
// RETURNS: inverse of this
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::inverse()
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::inverse()");

  GTMatrix<T> *mret = new GTMatrix<T>(n1_,n2_);

  if ( !inverse(*mret) ) {
    std::cout << "GTMatrix<T>::inverse(3): failed" << std::endl;
    exit(1);
  }

  return *mret;

} // end, method inverse (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : isSymmetric 
// DESC   : determines if matrix is symmetric
// ARGS   :
// RETURNS: TRUE or FALSE 
//**********************************************************************************
template<class T> 
GBOOL  GTMatrix<T>::isSymmetric()
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::inverse()");

  GSIZET i, j;
  GBOOL  bRet;

  if ( n1_ != n2_ ) return FALSE;

  // NOTE: should be symmetric w.r.t some tolerance!!!
  for ( j=0; j<n2_; j++ ) {
    for ( i=j+1,bRet=TRUE; i<n1_-1; i++ ) {
      bRet = bRet && (*this)(i,j)  == (*this)(j,i);
    }
  }
  return bRet;
 
} // end, method isSymmetric


//**********************************************************************************
//**********************************************************************************
// METHOD : improve
// DESC   : Taken largely from Numerical Recipes in C++, p. 59: Iterate to improve the
//          solution to Ax = b for x, given A, its LU decomposition, and b, and the 
//          initial solution, x
// ARGS   :
// RETURNS: GBOOL TRUE on success; else FALSE
//**********************************************************************************
template<class T>
GBOOL  GTMatrix<T>::improve(T **&a, T **alud, GSIZET n, GLLONG *&indx, T b[], T x[])
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::improve()");

  T      sdp, *r=new T[n];
  GLLONG i, j, m=0, nloop=1;
  if ( a == NULL || alud == NULL || indx == NULL ) {
    delete [] r; 
    return FALSE;
  }
  while ( m < nloop ) {
    for ( i=0; i<n; i++ ) {
      sdp = -b[i];
      for ( j=0; j<n; j++ ) 
        sdp += a[i][j] * x[j];
      r[i] = sdp;
    }  
    lubksb(alud, n, indx, r);
    for ( i=0; i<n; i++ ) x[i] -= r[i];
   m++;
  }
    
  delete [] r; 

  return TRUE;
} // end, method improve


//**********************************************************************************
//**********************************************************************************
// METHOD : wludcmp
// DESC   : Taken largely from Numerical Recipes
// ARGS   :
// RETURNS: GBOOL flag
//**********************************************************************************
template<class T>
GBOOL  GTMatrix<T>::wludcmp(T **&a, GSIZET   n, GLLONG   *&indx, T *d)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::wludcmp()");

  if ( a == NULL || indx == NULL ) return FALSE;

  GLLONG i, imax, j, k;
  GBOOL  bRet=TRUE;
  T      big, dum, sum, temp, gtiny=std::numeric_limits<T>::min();
  T     *vv;

  vv = new T [n]; 
  *d=1.0;
  for ( i=0; i<n && bRet; i++ ) {
    big=0.0;
    for ( j=0; j<n; j++ )
      if ( (temp=fabs(a[i][j])) > big ) big=temp;
    if ( big == 0.0 ) {
      bRet = FALSE; 
      break;
    }
    vv[i]=1.0 / big;
  }

  for ( j=0; j<n && bRet; j++ ) {
    for ( i=0; i<j; i++ ) {
      sum=a[i][j];
      for ( k=0; k<i; k++ ) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for ( i=j; i<n; i++ ) {
      sum=a[i][j];
      for ( k=0; k<j; k++ ) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big ) {
        big=dum;
        imax=i;
      }
    }
    if ( j != imax ) {
      for ( k=0; k<n; k++ ) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
     *d = -(*d);
     vv[imax]=vv[j];
    }
    indx[j]=imax;
    if ( fabs(a[j][j]) < gtiny ) a[j][j]=gtiny;
//  if ( a[j][j] == 0.0 ) a[j][j]=gtiny;
    if ( j != n-1 ) {
      dum=1.0/(a[j][j]);
      for ( i=j+1; i<n; i++ ) a[i][j] *= dum;
    }
  }
  delete [] vv;

  return bRet;

} // end, method wludcmp


//**********************************************************************************
//**********************************************************************************
// METHOD : lubksb
// DESC   : Taken largely from Numerical Recipes
// ARGS   :
// RETURNS: GBOOL flag
//**********************************************************************************
template<class T>
GBOOL  GTMatrix<T>::lubksb(T **&a, GSIZET   nd, GLLONG   *&indx, T b[])
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::lubksb()");

  GSIZET n=nd-1; 
  GLLONG i, ii=-1, ip, j;
  T      sum;

#if 0
  for ( i=0; i <= n; i++ ) {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii > -1) {
          for (j = ii; j < i; j++) sum -= a[i][j]*b[j];
      }
      else { 
          if (sum) ii = i;
      }
      b[i] = sum;
  }

  for ( i=n; i>=0; i-- ) {
      sum=b[i];
      if ( i < n )
      for (j = i+1; j <= n; j++) sum -= a[i][j]*b[j];
      b[i] = sum /(a[i][i]);
  }
#else
  n = nd; 
  ii = 0;
  for ( i=0; i<n; i++ ) {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if ( ii != 0 ) {
          for ( j=ii-1; j<i; j++ ) {
            sum -= a[i][j]*b[j];
          }
      }
      else  if ( sum != 0.0 ) ii = i+1;
      b[i] = sum;
  }

  for ( i=n-1; i>=0; i-- ) {
      sum=b[i];
      for ( j=i+1; j<n; j++ ) {
        sum -= a[i][j]*b[j];
      }
      b[i] = sum /(a[i][i]);
  }
#endif

  return TRUE;

} // end, method lubksb


//**********************************************************************************
//**********************************************************************************
// METHOD : ludcmp
// DESC   : LU decomposition method provided by
//          Warren Jasper, NC State Univ.
// ARGS   :
// RETURNS: GBOOL flag
//**********************************************************************************
template<class T>
GBOOL  GTMatrix<T>::ludcmp(T **&a, GSIZET   nd, GLLONG *&indx, T *d)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::ludcmp()");

   GLLONG  i, imax, j, k;
   GSIZET  n=nd-1;
   GBOOL  bRet = TRUE;
   T      big, dum, sum, temp, gtiny=1.0e-20;
   T     *vv = new T [nd];

   *d = 1.0;
   imax = -1;

   for ( i=0; i<=n; i++ ) {
       big = 0.0;
       for ( j=0; j<=n; j++ )
           if ( (temp = fabs( a[i][j] ) ) > big) big = temp;
           if ( big == 0.0 ) {
             bRet = FALSE;  
             std::cout << "GTMatrix::ludcmp: big = 0" << std::endl;
             break;
           }
       vv[i] = 1.0 / big;
   }
   for ( j=0; j<=n; j++ ) {
       for ( i=0; i<j; i++ ) {
           sum = a[i][j];
               for ( k=0; k<i; k++ ) sum -= a[i][k]*a[k][j];
               a[i][j] = sum;
       }
       big = 0.0;
       for ( i=j; i<=n; i++ ) {
           sum = a[i][j];
           for ( k=0; k<j; k++ )
               sum -= a[i][k]*a[k][j];
           a[i][j] = sum;
           if ( (dum = vv[i]*fabs(sum)) >= big ) {
               big = dum;
               imax = i;
           }
       }
       if (j != imax) {
           if ( imax < 0 ) { std::cout << "GTMatrix::ludcmp: imax < 0 " << std::endl; return FALSE;}
           for ( k=0; k<=n; k++ ) {
               dum = a[imax][k];
               a[imax][k] = a[j][k];
               a[j][k] = dum;
           }
           *d = -(*d);
           vv[imax] = vv[j];
       }
       indx[j] = imax;
       if ( a[j][j] == 0.0 ) a[j][j] = gtiny;
       if ( j != n ) {
           dum = 1.0 / (a[j][j]);
           for ( i=j+1; i<=n; i++ ) a[i][j] *= dum;
       }
   }
   delete [] vv;

   return bRet;

} // end, method ludcmp

//**********************************************************************************
//**********************************************************************************
// METHOD : isamax
// DESC   : Routine that finds the index of element having max.
//          absolute value.
// ARGS   : n   : Number of elements to check.
//          sx  : Vector to be checked.
//          incx: Every incx-th element is checked.
// RETURNS: int
//**********************************************************************************
template<class T> 
GSIZET  GTMatrix<T>::isamax(GSIZET n, T *sx, GSIZET incx)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::isamax()");

  T       smax = std::numeric_limits<T>::min();
  GLLONG  i, istmp = 0;

  if( n <= 1 ) return( istmp );
  if( incx != 1 ) {
    // Code for increment not equal to 1. 
//  if( incx < 0 ) sx = sx + ((-n+1)*incx + 1);
    istmp = 0;
    smax  = fabs( *sx );
    sx += incx;
    for( i=1; i<n; i++, sx+=incx )
      if( fabs( *sx ) > smax ) {
        istmp = i;
        smax  = fabs( *sx );
      }
    return( istmp );
  }
  // Code for increment equal to 1.
  istmp = 0;
  smax  = fabs(*sx);
  sx++;
  for( i=1; i<n; i++, sx++ )
    if( fabs( *sx ) > smax ) {
      istmp = i;
      smax  = fabs( *sx );
    }
  return( istmp );
}  // end, method isamax


//**********************************************************************************
//**********************************************************************************
// METHOD : svdcmp (1)
// DESC   : Taken from Numerical Recipes, p.70-72, find 
//          solution to Ax = b for x, given A, using Singular Value Decomp,
//          s.t. A = U W V^T. U replaces A on output. W is a diagonal
//          matrix of singular values represented as a vector in [0, ..., n-1].
//          V (not V^T) is output as an array v[0,...,n-1][0,...,n-1]
// ARGS   : 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<class T>
GBOOL GTMatrix<T>::svdcmp(T **a, GSIZET m, GSIZET n, T w[], T **v)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::svdcmp(T **,GSIZET,GSIZET,T *,T **)");

  GLLONG  i,its,j,jj,k,l, niter=50, nm;
  GBOOL   bRet=TRUE, flag;
  T       anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  GString serr = " GBOOL GTMatrix<T>::svdcmp (1): ";

   rv1 = new T [n];
   g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -GSIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -GSIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
//
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
//
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<niter;its++) {
//-----
			flag=TRUE;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag=FALSE;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == (niter-1)) {
                           std::cout << serr << "no convergence in " << niter << " iterations" << std::endl;
                           return FALSE;
                        }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+GSIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
        delete [] rv1;

  return bRet;
} // end, method svdcmp (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : svdcmp (2)
// DESC   : Taken from Numerical Recipes, p.70-72, find 
//          solution to Ax = b for x, given A, using Singular Value Decomp,
//          s.t. A = U W V^T. U replaces A on output. W is a diagonal
//          matrix of singular values represented as a vector in [0, ..., n-1].
//          V (not V^T) is output as an array v[0,...,n-1][0,...,n-1]
// ARGS   : 
//          w   : vector in which to store e-values; size of at least dim(diag(this))
//          v   : matrix of same size as a containing values of V
//          rv1 : GVector of dimension diag(this) 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<class T>
GBOOL GTMatrix<T>::svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::svdcmp(GTVector<T> &,GTMatrix<T> &,GTVector<T>&)");

	GLLONG  i,its,j,jj,k,l, m, n, niter=50, nm;
        GBOOL   bRet=TRUE, flag;
	T       anorm,c,f,g,h,s,scale,x,y,z;
        GString serr = " GBOOL GTMatrix<T>::svdcmp (2): ";

        m = this->dim(1);
        n = this->dim(2);
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs((*this)(k,i)); 
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					(*this)(k,i) /= scale;
					s += (*this)(k,i)*(*this)(k,i);
				}
				f=(*this)(i,i);
				g = -GSIGN(sqrt(s),f);
				h=f*g-s;
				(*this)(i,i)=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
					f=s/h;
					for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
				}
				for (k=i;k<m;k++) (*this)(k,i) *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs((*this)(i,k));
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					(*this)(i,k) /= scale;
					s += (*this)(i,k)*(*this)(i,k);
				}
				f=(*this)(i,l-1);
				g = -GSIGN(sqrt(s),f);
				h=f*g-s;
				(*this)(i,l-1)=f-g;
				for (k=l-1;k<n;k++) rv1[k]=(*this)(i,k)/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += (*this)(j,k)*(*this)(i,k);
					for (k=l-1;k<n;k++) (*this)(j,k) += s*rv1[k];
				}
				for (k=l-1;k<n;k++) (*this)(i,k) *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
//
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) v(j,i)=((*this)(i,j)/(*this)(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += (*this)(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) (*this)(i,j)=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
//
				for (s=0.0,k=l;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
				f=(s/(*this)(i,i))*g;
				for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
			}
			for (j=i;j<m;j++) (*this)(j,i) *= g;
		} else for (j=i;j<m;j++) (*this)(j,i)=0.0;
		++(*this)(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<niter;its++) {
//-----
			flag=TRUE;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag=FALSE;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=(*this)(j,nm);
						z=(*this)(j,i);
						(*this)(j,nm)=y*c+z*s;
						(*this)(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == (niter-1)) {
                           std::cout << serr << "no convergence in " << niter << " iterations" << std::endl;
                           return FALSE;
                        }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+GSIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=(*this)(jj,j);
					z=(*this)(jj,i);
				        (*this)(jj,j)=y*c+z*s;
					(*this)(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

  return bRet;
} // end, method svdcmp (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : svdcmp (3)
// DESC   : Taken from Numerical Recipes, p.70-72, find 
//          solution to Ax = b for x, given A, using Singular Value Decomp,
//          s.t. A = U W V^T. U replaces A on output. W is a diagonal
//          matrix of singular values represented as a vector in [0, ..., n-1].
//          V (not V^T) is output as an array v[0,...,n-1][0,...,n-1]
// ARGS   : 
//          w    : vector in which to store e-values; size of at least dim(diag(this))
//          v    : matrix of same size as a containing values of V
//          rv1  : GVector of dimension diag(this) 
//          nx,ny: valid dims of matrix on which to operate; may be less than nominal
//                 matrix dimensions. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<class T>
GBOOL GTMatrix<T>::svdcmp(GTVector<T> &w, GTMatrix<T> &v,  GTVector<T> &rv1, GSIZET nx, GSIZET ny)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::svdcmp(GTVector<T> &,GTMatrix<T> &,GTVector<T>&,GSIZET,GSIZET)");

	GLLONG  i,its,j,jj,k,l, m, n, niter=50, nm;
        GBOOL   bRet=TRUE, flag;
	T       anorm,c,f,g,h,s,scale,x,y,z;
        GString serr = " GBOOL GTMatrix<T>::svdcmp (2): ";

        m = nx;
        n = ny;
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs((*this)(k,i)); 
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					(*this)(k,i) /= scale;
					s += (*this)(k,i)*(*this)(k,i);
				}
				f=(*this)(i,i);
				g = -GSIGN(sqrt(s),f);
				h=f*g-s;
				(*this)(i,i)=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
					f=s/h;
					for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
				}
				for (k=i;k<m;k++) (*this)(k,i) *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs((*this)(i,k));
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					(*this)(i,k) /= scale;
					s += (*this)(i,k)*(*this)(i,k);
				}
				f=(*this)(i,l-1);
				g = -GSIGN(sqrt(s),f);
				h=f*g-s;
				(*this)(i,l-1)=f-g;
				for (k=l-1;k<n;k++) rv1[k]=(*this)(i,k)/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += (*this)(j,k)*(*this)(i,k);
					for (k=l-1;k<n;k++) (*this)(j,k) += s*rv1[k];
				}
				for (k=l-1;k<n;k++) (*this)(i,k) *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
//
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) v(j,i)=((*this)(i,j)/(*this)(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += (*this)(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) (*this)(i,j)=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
//
				for (s=0.0,k=l;k<m;k++) s += (*this)(k,i)*(*this)(k,j);
				f=(s/(*this)(i,i))*g;
				for (k=i;k<m;k++) (*this)(k,j) += f*(*this)(k,i);
			}
			for (j=i;j<m;j++) (*this)(j,i) *= g;
		} else for (j=i;j<m;j++) (*this)(j,i)=0.0;
		++(*this)(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<niter;its++) {
//-----
			flag=TRUE;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag=FALSE;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=(*this)(j,nm);
						z=(*this)(j,i);
						(*this)(j,nm)=y*c+z*s;
						(*this)(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == (niter-1)) {
                           std::cout << serr << "no convergence in " << niter << " iterations" << std::endl;
                           return FALSE;
                        }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+GSIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=(*this)(jj,j);
					z=(*this)(jj,i);
				        (*this)(jj,j)=y*c+z*s;
					(*this)(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

  return bRet;
} // end, method svdcmp (3)



//**********************************************************************************
//**********************************************************************************
// METHOD : jacobi (1)
// DESC   : Taken from Numerical Recipes, p.346-348, find 
//          eigen-decomposition for square matrix. Operates on *this.
// ARGS   : 
//          d   : vector of eigenvalues, returned
//          v   : matrix of same size as *this whose columns contain normalized e-vectors, returned
//                nrot: no. Jacobi rotations used, returned 
//          a   : tmp matrix of same size as *this
//          b   : tmp vector of same size as diagnonal of *this
//          z   : tmp vector of same size as diagnonal of *this
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<class T>
void GTMatrix<T>::jacobi(GTVector<T> &d, GTMatrix<T> &v, GSIZET &nrot, GTMatrix<T> &a, GTVector<T> &b, GTVector<T> &z)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::jacobi");

  GLLONG  j,iq,ip,i,m,n;
  GDOUBLE tresh,theta,tau,t,sm,s,h,g,c;
  GString serr = "GTMatrix<T>::jacobi (1): ";

  m = this->dim(1);
  n = this->dim(2);
  if ( m != n ) {
    std::cout << serr << "matrix must be square" << std::endl;
    exit(1);
  }

  a = *this; // make copy

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v(ip,iq)=0.0;
    v(ip,ip)=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a(ip,ip);
    z[ip]=0.0;
  }
  nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a(ip,iq));
      }
      if (sm == 0.0) {
        return;
      }
      if (i < 4)
        tresh=0.2*sm/(n*n);
      else
        tresh=0.0;
      for (ip=0;ip<n-1;ip++) {
        for (iq=ip+1;iq<n;iq++) {
	  g=100.0*fabs(a(ip,iq));
          if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
                    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
            a(ip,iq)=0.0;
          else if (fabs(a(ip,iq)) > tresh) {
            h=d[iq]-d[ip];
	    if ((float)(fabs(h)+g) == (float)fabs(h))
              t=(a(ip,iq))/h;
            else {
              theta=0.5*h/(a(ip,iq));
              t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
              if (theta < 0.0) t = -t;
            }
            c=1.0/sqrt(1+t*t);
            s=t*c;
            tau=s/(1.0+c);
            h=t*a(ip,iq);
            z[ip] -= h;
            z[iq] += h;
            d[ip] -= h;
            d[iq] += h;
            a(ip,iq)=0.0;
            for (j=0;j<ip;j++) {
              GTMATRIX_ROTATE(a,j,ip,j,iq)
            }
            for (j=ip+1;j<iq;j++) {
              GTMATRIX_ROTATE(a,ip,j,j,iq)
            }
            for (j=iq+1;j<n;j++) {
              GTMATRIX_ROTATE(a,ip,j,iq,j)
            }
            for (j=0;j<n;j++) {
              GTMATRIX_ROTATE(v,j,ip,j,iq)
            }
            ++nrot;
          }
        }
      }
      for (ip=0;ip<n;ip++) {
        b[ip] += z[ip];
        d[ip]=b[ip];
        z[ip]=0.0;
      }
  }
  std::cout << serr << "Too many iterations" << std::endl;

} // end, method jacobi (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : setSingularZero
// DESC   : Serves as public method set tolerance at which the singular (eigen-) values
//          from SVD method are zeroed out.
// ARGS   : T  tol 
// RETURNS: none
//**********************************************************************************
template<class T>
void GTMatrix<T>::setSingularZero(GDOUBLE tol)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::setSingularZero");

  singzero_ = tol;

} // end, method setSingularZero


//**********************************************************************************
//**********************************************************************************
// METHOD : set
// DESC   : assings matrix to constant value
// ARGS   : T  value
// RETURNS: none
//**********************************************************************************
template<class T> void GTMatrix<T>::set(T a)
{
  data_.set(a);

} // end, method set


//**********************************************************************************
//**********************************************************************************
// METHOD : set
// DESC   : assings matrix data from input array, and reshapesaccording
//          to internal data format (row or column major ordering). Matrix
//          will be resized to accomodate this data.
// ARGS   : T *: input array of data
//          n1 : number of rows
//          n2 : number of columns
// RETURNS: none
//**********************************************************************************
template<class T> void GTMatrix<T>::set(T *array, GSIZET n1, GSIZET n2)
{
  
  n1_ = n1;
  n2_ = n2;
  data_.resize(n1_*n2_);
  for ( GSIZET j=0; j< n1_*n2_; j++ ) data_[j] = array[j];

} // end, method set


//**********************************************************************************
//**********************************************************************************
// METHOD : dpythag
// DESC   : Taken from Numerical Recipes and serves as a utility for svdcmp method
//          Computes (a^2 + b^2)^1/2 without underflow or overflow.
// ARGS   : 
// RETURNS: T result
//**********************************************************************************
template<class T>
T GTMatrix<T>::dpythag(T a, T b)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::dpythag");

  T    absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+pow((GDOUBLE)(absb/absa),2.0));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow((GDOUBLE)(absa/absb),2.0)));
} // end, method dpythag


//**********************************************************************************
//**********************************************************************************
// METHOD : svbksub
// DESC   : Taken from Numerical Recipes and serves back substitution 
//              method for SVD
// ARGS   : 
// RETURNS: T result
//**********************************************************************************
template<class T>
void GTMatrix<T>::svbksub(T **u, T w[], T **v, GSIZET n, GSIZET m, T b[], T x[])
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::svbksub");

  GLLONG  jj,j,i;
  T       s,*tmp;

  tmp = new T [n];
  for (j=0;j<n;j++) {
    s=0.0;
    if (w[j] != 0.0) {
      for (i=0;i<m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  delete [] tmp;
} // end, method svbksub


//**********************************************************************************
//**********************************************************************************
// METHOD : choldc (1)
// DESC   : Taken from Numerical Recipes: serves as method to perform
//          efficient Cholesky decomposition of an input matrix, a, s.t.
//          a = L L^T. Upper diag. of a isn't touched, and L is
//          stored in the lower diag. part of a. The diag elements
//          are stored in p, and set in a. 
// ARGS   : 
// RETURNS:  T result
//**********************************************************************************
template<class T>
void GTMatrix<T>::choldc(T **a, T p[], GSIZET n)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::choldc(T **,T *, GSIZET)");

  GLLONG i, j, k;
  T      sum;
  T      gtiny = std::numeric_limits<T>::min();

  for (i=0;i<n;i++) {
    for ( j=i; j<n;j++ ) {
      for ( sum=a[i][j], k=i-1; k>=0; k--) sum -= a[i][k]*a[j][k];
        if ( i == j ) {
          if ( sum <= 0.0 ) { std::cout << "GTMatrix<T>::choldc (1): failed" << std::endl; exit(1); }
          p[i] = sqrt(sum);
        } else a[j][i] = sum / (p[i]+gtiny);
    }
  }

  // By default, add the diagonal back into matrix:
  for (i=0;i<n;i++) a[i][i] = p[i];

} // end, method choldc (1)



//**********************************************************************************
//**********************************************************************************
// METHOD   : choldc (2)
// DESC     : Taken from Numerical Recipes: serves as method to perform
//            efficient Cholesky decomposition of an input matrix, that is
//            symmetric positive definite. That is, 
//            B (=this)  = L L^T.  The result is returned in the 
//            lower triangular part of A; the diagonal part is returned in p, and
//            this is set in A. The upper diagonal part of A is zeroed.
// ARGS     : 
// RETURNS  : T result
//**********************************************************************************
template<class T>
void GTMatrix<T>::choldc(GTMatrix<T> &a, GTVector<T> &p)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::choldc(GTMatrix<T>&,GTVector<T>&)");

  GLLONG  i, j, k, n;
  GDOUBLE sum;
  T       gtiny = std::numeric_limits<T>::min();

  a = 0.0;
  n = this->dim(1);
  for (i=0;i<n;i++) {
    for ( j=i;j<n;j++ ) {
      for ( sum=(*this)(i,j), k=i-1; k>=0; k--) sum -= (*this)(i,k) * (*this)(j,k); 
      if ( i == j ) {
        if ( sum <= 0.0 ) { 
          std::cout << "GTMatrix<T>::choldc (2): failed" << std::endl;  
          std::cout << "a=" << a << std::endl; 
          exit(1); }
      p[i] = sqrt(sum);
      } else a(j,i) = sum / (p[i]+gtiny);
    }
  }
 
  // By default, add the diagonal back into matrix:
  for (i=0;i<n;i++) a(i,i) = p[i];
 
} // end, method choldc (2)
 
 
//**********************************************************************************
//**********************************************************************************
// METHOD : createIdentity
// DESC   : Creates identity from the matrix
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTMatrix<T>::createIdentity()
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::createIdentity");

  GString serr = "GTMatrix<T>::createIdentity: ";
 
  if ( n1_ != n2_ ) {
    std::cout << serr << "matrix not square" << std::endl;
     exit(1);
  }

  (*this) = 0;
  for ( GSIZET i=0; i<n1_; i++ ) (*this)(i,i) = 1;
} // end, method createIdentity


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity
// DESC   : Find multiplicity of specified value
// ARGS   : val  : value to find multiplicity of. If val no found in matrix, then
//                 return is 0
// RETURNS: multiplicity of val
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::multiplicity(T val)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::multiplicity(T)");
 
  return data_.multiplicity(val);

} // end, method multiplicity


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity
// DESC   : Find multiplicity of specified value, and return the 
//          indices of the occurrences in the index array. 
// ARGS   : val  : value to find multiplicity of. If val no found in matrix, then
//                 return is 0
//          irow : Resized if mult > n to contain rows where 'val' appears
//          icol : Resized if mult > n to contain columns where 'val' appears.
//          n    : size of irow/icol arrays on output if n < multiplicity; else
//                 n = multiplicity.
// RETURNS: multiplicity of val, which may be 0
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::multiplicity(T val, GSIZET *&irow, GSIZET *&icol, GSIZET &n)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::multiplicity(T,GSIZET*,GZSIZET*,GSIZET&)");

  GSIZET mult = data_.multiplicity(val);

  if ( mult == 0 ) return mult;

  if ( n < mult ) {
    delete [] irow;
    delete [] icol;
    n    = mult;
    irow = new GSIZET [n];
    icol = new GSIZET [n];
  }
 
  GSIZET m;
  for ( GSIZET i=0,m=0; i<n1_; i++ ) {
    for ( GSIZET j=0; j<n2_; j++ ) {
    
      if ( (*this)(i,j) == val ) {
        irow[m] = i ; icol[m] = j; 
        m++;
      }
    }
  }

} // end, method multiplicity


//**********************************************************************************
//**********************************************************************************
// METHOD : multiplicity
// DESC   : Find multiplicity of specified value, and return the 
//          indices of the occurrences in the index array. 
// ARGS   : val  : value to find multiplicity of. If val no found in matrix, then
//                 return is 0
//          irow : Resized if mult > n to contain rows where 'val' appears
//          icol : Resized if mult > n to contain columns where 'val' appears.
// RETURNS: multiplicity of val, which may be 0
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::multiplicity(T val, GSZBuffer &irow, GSZBuffer &icol)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::multiplicity(T,GSIZEBuffer&,GSZBuffer&)");
 
  GSIZET mult = data_.multiplicity(val);

  if ( mult == 0 ) return mult;
 
  if ( irow.size() < mult ) {
    irow.resize(mult);
  }
  if ( icol.size() < mult ) {
    icol.resize(mult);
  }
 
  GSIZET m;
  for ( GSIZET i=0,m=0; i<n1_; i++ ) {
    for ( GSIZET j=0; j<n2_; j++ ) {
    
      if ( (*this)(i,j) == val ) {
        irow[m] = i ; icol[m] = j; 
        m++;
      }
    }
  }
 
  return mult;
 
} // end, method multiplicity


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrow
// DESC   : Find for a given row, the distinct values in that row. 
// ARGS   :
//          irow : Specified row
//          vals : distinct values found for row irow. Will be resized if
//                 number distinct vals found > n.
//          icol : Resized if mult > n to contain columns where 'vals' appear.
//          n    : size of vals and icol arrays on output if n < number distinct
//                 values found, nv, else, n = nv
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctrow( GSIZET irow, T *&val, GSIZET *&icol, GSIZET &n)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::distinctrow(GSIZET,T*,GSIZET*,GSIZET&)");

  GSIZET nd = data_.distinctrng(irow*n2_, n1_*n2_, n1_, val, icol, n); 

  return nd;

} // end, method distinctrow


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrow
// DESC   : Find for a given row, the distinct values in that row. 
// ARGS   :
//          irow : Specified row
//          vals : distinct values found for row irow. Will be resized if
//                 number distinct vals found > n.
//          icol : Resized if mult > current size to contain columns where 'vals' appear.
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctrow( GSIZET irow, GTVector<T> &val, GSZBuffer &icol)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::distinctrow(GSIZET,T*,GTVector<T>&,GSZBuffer&)");

  GSIZET *index=0, nd, nv=0;
  T      *vvals=0;
 
  nd = data_.distinctrng(irow*n2_, n1_*n2_, n1_, vvals, index, nv); 
  
  if ( val .size() < nd || icol.size() < nd ) {
    val .resize(nd); 
    icol.resize(nd); 
  }
 
  for ( GINT j=0; j<nd; j++ ) {
    val [j] = vvals[j];
    icol[j] = index[j];
  }
  
  if ( vvals != NULLPTR ) delete [] vvals;
  if ( index != NULLPTR ) delete [] index;
 
  return nd;

} // end, method distinctrow


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrow_floor
// DESC   : Find for a given row, the distinct values in that row, considering
//          only values > floor
// ARGS   :
//          irow : Specified row
//          vals : distinct values found for row irow. Will be resized if
//                 number distinct vals found > n.
//          icol : Resized if mult > n to contain columns where 'vals' appear.
//          n    : size of vals and icol arrays on output if n < number distinct
//                 values found, nv, else, n = nv
//          floor: all distinct vals must be > floor
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctrow_floor( GSIZET irow, T *&val, GSIZET *&icol, GSIZET &n, T floor)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::distinctrow_floor(GSIZET,T*,GSIZET*,GSIZET&)");

  GSIZET nd = data_.distinctrng_floor(irow*n2_, n1_*n2_, n1_, val, icol, n, floor); 

  return nd;

} // end, method distinctrow_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctrow_floor
// DESC   : Find for a given row, the distinct values in that row, considering
//          only values > floor
// ARGS   :
//          irow : Specified row
//          vals : distinct values found for row irow. Will be resized if
//                 number distinct vals found > n.
//          icol : Resized if mult > icol size to contain columns where 'vals' appear.
//          floor: all distinct vals must be > floor
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctrow_floor( GSIZET irow, GTVector<T> &val, GSZBuffer &icol, T floor)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::distinctrow_floor(GSIZET,T*,GSIZET*,GSIZET&)");

  T      *vval=0;
  GSIZET  nm=0;
  GSIZET *vicol=0;
  GSIZET nd = data_.distinctrng_floor(irow*n2_, n1_*n2_, n1_, vval, vicol, nm, floor); 
  if ( val.size() < nd || icol.size() < nd ) {
    val .resize(nd); 
    icol.resize(nd); 
  }
  val .set(vval ,nd);
  icol.set(vicol,nd);
  if ( vval  != NULLPTR ) delete [] vval;
  if ( vicol != NULLPTR ) delete [] vicol;

  return nd;

} // end, method distinctrow_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctcol
// DESC   : Find for a given column, the distinct values in that column. 
// ARGS   :
//          icol : Specified col
//          vals : distinct values found for col icol. Will be resized if
//                 number distinct vals found > n.
//          irow : Resized if mult > n to contain rows where 'vals' appear.
//          n    : size of vals and irow arrays on output if n < number distinct
//                 values found, nv, else, n = nv
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctcol( GSIZET icol, T *&val, GSIZET *&irow, GSIZET &n)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::distinctcol(GSIZET,T*,GSIZET*,GSIZET&)");

  // Since data is col-major, can locate contiguous data and use GTVector function:
  GSIZET nd = data_.distinctrng(icol*n1_, n1_, 1, val, irow, n); 

  return nd;

} // end, method distinctcol


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctcol
// DESC   : Find for a given column, the distinct values in that column. 
// ARGS   :
//          icol : Specified column
//          vals : distinct values found for col icol. Will be resized if
//                 number distinct vals found > n.
//          irow : Resized if mult > current size to contain rows where 'vals' appear.
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctcol( GSIZET icol, GTVector<T> &val, GSZBuffer &irow)
{
  static_assert(std::is_arithmetic<T>::value || std::is_pointer<T>::value,
    "Invalid template type: GTMatrix<T>::distinctcol(GSIZET,T*,GTVector<T>&,GSZBuffer&)");

  GSIZET *index=0, nd, nv=0;
  T      *vvals=0;
 
  // Since data is col-major, can locate contiguous data and use GTVector function:
  nd = data_.distinctrng(icol*n1_, n1_, 1, vvals, index, nv); 
  
  if ( val .size() < nd || irow.size() < nd ) {
    val .resize(nd); 
    irow.resize(nd); 
  }
 
  for ( GSIZET j=0; j<nd; j++ ) {
    val [j] = vvals[j];
    irow[j] = index[j] - icol*n1_; // make a row index out of linearized index
  }
  
  if ( vvals != NULLPTR ) delete [] vvals;
  if ( index != NULLPTR ) delete [] index;
 
  return nd;

} // end, method distinctcol


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctcol_floor
// DESC   : Find for a given column, the distinct values in that column, considering
//          only values > floor
// ARGS   :
//          icol : Specified col 
//          vals : distinct values found for col icol. Will be resized if
//                 number distinct vals found > n.
//          irow : Resized if mult > n to contain rows where 'vals' appear.
//          n    : size of vals and irow arrays on output if n < number distinct
//                 values found, nv, else, n = nv
//          floor: all distinct vals must be > floor
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctcol_floor( GSIZET icol, T *&vals, GSIZET *&irow, GSIZET &n, T floor)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::distinctcol_floor(GSIZET,T*,GSIZET*,GSIZET&)");

  // Since data is col-major, can locate contiguous data and use GTVector function:
  GSIZET nd = data_.distinctrng_floor(icol*n1_, n1_, 1, vals, irow, n, floor); 

  return nd;

} // end, method distinctcol_floor


//**********************************************************************************
//**********************************************************************************
// METHOD : distinctcol_floor
// DESC   : Find for a given column , the distinct values in that column, considering
//          only values > floor
// ARGS   :
//          icol : Specified column
//          vals : distinct values found for row irow. Will be resized if
//                 number distinct vals found > n.
//          irow : Resized if mult > irow size to contain rows where 'vals' appear.
//          floor: all distinct vals must be > floor
// RETURNS: nv = number of distinct values found (may not be 0)
//**********************************************************************************
template<class T>
GSIZET GTMatrix<T>::distinctcol_floor( GSIZET icol, GTVector<T> &vals, GSZBuffer &irow, T floor)
{
  static_assert(std::is_arithmetic<T>::value,
    "Invalid template type: GTMatrix<T>::distinctcol_floor(GSIZET,T*,GSIZET*,GSIZET&)");

  T      *vvals=0;
  GSIZET  nm=0;
  GSIZET *vind=0;
  
  // Since data is row-major, can locate contiguous data and use GTVector function:
  GSIZET nd = data_.distinctrng_floor(icol*n1_, n1_, 1, vvals, vind, nm, floor); 
  if ( vals.size() < nd || irow.size() < nd ) {
    vals.resize(nd); 
    irow.resize(nd); 
  }
  vals .set(vvals ,nd);
  for ( GSIZET i=0; i<nd; i++ ) irow[i] = vind[i] - icol*n1_; // make a row index out of linearized index
  if ( vvals != NULLPTR ) delete [] vvals;
  if ( vind  != NULLPTR ) delete [] vind;

  return nd;

} // end, method distinctcol_floor



//**********************************************************************************
//**********************************************************************************
// METHOD : setCacheBlockingFactor
// DESC   : Set blocking factor
// ARGS   : GINT number of T type members in cache
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTMatrix<T>::setCacheBlockingFactor(GINT icsz)
{
  icsz_ = MAX(icsz,0);

} // end, method SetCacheBlocking


//**********************************************************************************
//**********************************************************************************
// METHOD : updatehost
// DESC   : Update data from device to host
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T>
void GTMatrix<T>::updatehost()
{
  data_.updatehost();
} // end, method updatehost


//**********************************************************************************
//**********************************************************************************
// METHOD : updatedev
// DESC   : Update data from host to device
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<class T> 
void GTMatrix<T>::updatedev()
{
  data_.updatedev();
} //nd of method updatedev


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * mat-vec
// DESC   : matrix-vector product, returns product,
//          without destroying *this data, for GFLOAT type
// ARGS   : GTVector &
// RETURNS: GTVector vector
//**********************************************************************************
template<class T>
GTVector<T> &GTMatrix<T>::operator*(GTVector<T> &a)
{
  return this->matvec_impl_(a, typename std::is_floating_point<T>::type());
} // end, method operator*


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * mat-mat
// DESC   : Multiplies this by matrix m, and returns
//          result, for GFLOAT types
// ARGS   : GTMatrix a factor
// RETURNS: GTMatrix product matrix
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::operator*(GTMatrix<T> &a)
{
  return this->matmat_impl_(a, typename std::is_floating_point<T>::type());

} // end, method operator*


//**********************************************************************************
//**********************************************************************************
// METHOD : matvec_impl_ (1)
// DESC   : 
// ARGS   : GTVector &
//          std::type telling us this is for a 'false type' (so, T is
//                a floating point type)
// RETURNS: GTMatrix &
//**********************************************************************************
template<class T>
GTVector<T> &GTMatrix<T>::matvec_impl_(GTVector<T> &obj, std::true_type)
{
  GTVector<T> *vret = new GTVector<T> (this->size(1));

  GMTK::matvec_prod(*vret, *this, obj);

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;
} // end, matvec_impl_ (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : matvec_impl_ (2)
// DESC   : 
// ARGS   : obj      : GTVector & operand
//          std::type: telling us this is for a 'false type' (so, T is
//                not a floating point type)
// RETURNS: GTMatrix &
//**********************************************************************************
template<class T>
GTVector<T> &GTMatrix<T>::matvec_impl_(GTVector<T> &obj, std::false_type)
{
  GTVector<T> *vret = new GTVector<T> (this->size(1));

  for ( GSIZET i=0; i<this->size(1); i++ ) {
    (*vret)[i] = static_cast<T>(0);
    for ( GSIZET j=0; j<this->size(2); j++ ) {
      (*vret)[i] += (*this)(i,j)*obj[j];
    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *vret;
} // end, matvec_impl_ (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : matmat_impl_ (1)
// DESC   : 
// ARGS   : obj      : GTMatrix & mat operatand
//          std::type: telling us this is for a 'false type' (so, T is
//                     a floating point type)
// RETURNS: GTMatrix &
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::matmat_impl_(GTMatrix<T> &obj, std::true_type)
{
  GTMatrix<T> *mret = new GTMatrix<T> (this->size(1),obj.size(2));

  GMTK::matmat_prod(*mret, *this, obj);

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *mret;
} // end, matmat_impl_ (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : matmat_impl_ (2)
// DESC   : 
// ARGS   : obj      : GTMatrix & mat operatand
//          std::type: telling us this is for a 'false type' (so, T is
//                     not a floating point type)
// RETURNS: GTMatrix &
//**********************************************************************************
template<class T>
GTMatrix<T> &GTMatrix<T>::matmat_impl_(GTMatrix<T> &obj, std::false_type)
{
  GTMatrix<T> *mret = new GTMatrix<T> (this->size(1),obj.size(2));

  for ( GSIZET j=0; j<obj.size(2); j++ ) {
    for ( GSIZET i=0; i<this->size(1); i++ ) {

      (*mret)(i,j) = 0.0;
      for ( GSIZET k=0; k<this->size(2); k++ ) {
        (*mret)(i,j) += (*this)(i,k) * obj(k,j);
      }

    }
  }

  #if defined(_G_AUTO_UPDATE_DEV)
  vret->updatedev();
  #endif

  return *mret;
} // end, matmat_impl_ (2)



