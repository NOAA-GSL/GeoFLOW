//==================================================================================
// Module       : gmtk.hpp
// Date         : 1/31/18 (DLR)
// Description  : Math TooKit: namespace of C routines for various
//                math functions
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GMTK_HPP)
#define _GMTK_HPP

#include "gtypes.h"
#include "cff_blas.h"


template<typename T> class GTVector;
template<typename T> class GTMatrix;
                     class GGrid;

extern GINT szMatCache_;
extern GINT szVecCache_;

namespace GMTK
{


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod_k (1)
// DESC   : Compute cross/vector product with hat(k)
//             C = A X hat(k)
//          
// ARGS   : A    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          isgn : multiplies product by sign(isgn)
//          C    : array of pointers to cross produc vectors; must have
//                 at least 3 elements for 3-d vector products. All vector
//                 elements must have the same length, of at least length nind if
//                 iind != NULLPTR, and of length of x[?], y[?] if iind==NULLPTR.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod_k(GTVector<GTVector<T>*> &A, GINT *iind, GINT nind, GINT isgn, GTVector<GTVector<T>*> &C)
{
  assert( A.size() >= 2 && C.size() >= 2 &&  "Incompatible dimensionality");

  GSIZET n;
  T      fact = isgn < 0 ? -1.0 : 1.0;
  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      C[0][k] =  (*A[1])[n]*fact;
      C[1][k] = -(*A[0])[n]*fact;
      if ( C.size() > 2 ) C[2][n] = 0.0;
    }
  }
  else {
    for ( GSIZET n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
      C[0][n] =  (*A[1])[n]*fact;
      C[1][n] = -(*A[0])[n]*fact;
      if ( C.size() > 2 ) C[2][n] = 0.0;
    }
  }

} // end of method cross_prod_k (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod_k (2)
// DESC   : Compute cross/vector product with hat(k)
//             C = A X hat(k)
//          
// ARGS   : Ai   : Vector components x, y, z; must each have same no. elements.
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          isgn : multiplies product by sign(isgn)
//          Ci   : Vector components for solution, must each have the same no. elements
//                 as in Ai, Bi, unless iind != NULLPTR, in which case, Ci must each
//                 have at least nind elements.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod_k(GTVector<T> &Ax, GTVector<T> &Ay, 
                  GINT *iind, GINT nind, GINT isgn, 
                  GTVector<T> &Cx, GTVector<T> &Cy)
{

  GSIZET n;
  T      fact = isgn < 0 ? -1.0 : 1.0;
  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      Cx[k] =  Ay[n]*fact;
      Cy[k] = -Ax[n]*fact;
    }
  }
  else {
    for ( GSIZET n=0; n<Ax.size(); n++ ) { // cycle over all coord pairs
      Cx[n] =  Ay[n]*fact;
      Cy[n] = -Ax[n]*fact;
    }
  }

} // end of method cross_prod_k (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod (1)
// DESC   : compute cross/vector product 
//             C = A X B
//          
// ARGS   : A    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          B    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          C    : array of pointers to cross produc vectors; must have
//                 at least 3 elements for 3-d vector products. All vector
//                 elements must have the same length, of at least length nind if
//                 iind != NULLPTR, and of length of x[?], y[?] if iind==NULLPTR.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod(GTVector<GTVector<T>*> &A, GTVector<GTVector<T>*> &B, 
                GINT *iind, GINT nind, GTVector<GTVector<T>*> &C)
{
  assert( A.size() >= 3 && B.size() && C.size() >= 3 && "Incompatible dimensionality");

  GSIZET n;
  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      x1 = (*A[0])[n]; y1 = (*A[1])[n]; z1 = (*A[2])[n];
      x2 = (*B[0])[n]; y2 = (*B[1])[n]; z2 = (*B[2])[n];
      C[0][k] = y1*z2 - z1*y2; 
      C[1][k] = z1*x2 - z2*x1; 
      C[2][k] = x1*y2 - x2*y1;
    }
  }
  else {
    for ( GSIZET n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
      x1 = (*A[0])[n]; y1 = (*A[1])[n]; z1 = (*A[2])[n];
      x2 = (*B[0])[n]; y2 = (*B[1])[n]; z2 = (*B[2])[n];
      C[0][n] = y1*z2 - z1*y2; 
      C[1][n] = z1*x2 - z2*x1; 
      C[2][n] = x1*y2 - x2*y1;
    }
  }

} // end of method cross_prod (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod (2)
// DESC   : compute cross/vector product 
//             C = A X B
//          
// ARGS   : Ai   : Vector components x, y, z; must each have same no. elements.
// ARGS   : Bi   : Vector components x, y, z; must each have same no. elements, as Ai
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          Ci   : Vector components for solution, must each have the same no. elements
//                 as in Ai, Bi, unless iind != NULLPTR, in which case, Ci must each
//                 have at least nind elements.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod(GTVector<T> &Ax, GTVector<T> &Ay, GTVector<T> &Az,
                GTVector<T> &Bx, GTVector<T> &By, GTVector<T> &Bz,
                GINT *iind, GINT nind, 
                GTVector<T> &Cx, GTVector<T> &Cy, GTVector<T> &Cz)
{

  GSIZET n;
  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      x1 = Ax[n]; y1 = Ay[n]; z1 = Az[n];
      x2 = Bx[n]; y2 = By[n]; z2 = Bz[n];
      Cx[k] = y1*z2 - z1*y2; 
      Cy[k] = z1*x2 - z2*x1; 
      Cz[k] = x1*y2 - x2*y1;
    }
  }
  else {
    for ( GSIZET n=0; n<Ax.size(); n++ ) { // cycle over all coord pairs
      x1 = Ax[n]; y1 = Ay[n]; z1 = Az[n];
      x2 = Bx[n]; y2 = By[n]; z2 = Bz[n];
      Cx[n] = y1*z2 - z1*y2; 
      Cy[n] = z1*x2 - z2*x1; 
      Cz[n] = x1*y2 - x2*y1;
    }
  }

} // end of method cross_prod (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : normalize_euclidean
// DESC   : 
//             Compute Euclidean norm of each 'point', and return, overriding
//             input vectors
//          
// ARGS   : x    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void normalize_euclidean(GTVector<GTVector<T>*> &x, GINT *iind, GINT nind)
{
  GSIZET n;
  T      xn(3);

  // NOTE: The better way to do this, especially for long lists is use
  //       outer loop over coord dimensions, and store progressive sum
  //       over each for each tuple. But this requires vector storage for
  //       each coordinate for all tuples.
  xn = 0.0;
  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all n-tuples
      n = iind[k];
      for ( GSIZET l=0, xn=0.0; l<x.size(); l++ ) xn += (*x[l])[n]*(*x[l])[n];
      xn = 1.0/xn;
      for ( GSIZET l=0, xn=0.0; l<x.size(); l++ ) (*x[l])[n] *= xn;
    }
  }
  else {
    for ( GSIZET n=0; n<x[0]->size(); n++ ) { // cycle over all n-tuples
      for ( GSIZET l=0, xn=0.0; l<x.size(); l++ ) xn += (*x[l])[n]*x[l][n];
      xn = 1.0/xn;
      for ( GSIZET l=0, xn=0.0; l<x.size(); l++ ) (*x[l])[n] *= xn;
    }
  }

} // end of method normalize_euclidean


//**********************************************************************************
//**********************************************************************************
// METHOD :    saxpby
// DESC   : 
//             compute x = ax + by
//          
// ARGS   : x : vector , updated
//          a : const multiplying x
//          y : vector, must be same size as x
//          b : const multiplying y
//          
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void saxpby(GTVector<T> &x, T a, GTVector<T> &y, T b) 
{
  assert(x.size() == y.size() && "Incompatible array sizes");
  for ( GSIZET j=0; j<x.size(); j++ ) { 
    x[j] = a*x[j] + b*y[j];
  }
} // end of method saxpby


  template<typename T>     
  void D2_X_D1(GTMatrix<T> &D1, GTMatrix<T> &D2T, GTVector<T> &u, 
               GTVector<T> &tmp, GTVector<T> &y);
  template<typename T>     
  void D3_X_D2_X_D1(GTMatrix<T> &D1, GTMatrix<T> &D2T, GTMatrix<T> &D3T,  
                    GTVector<T> &u , GTVector<T> &tmp, GTVector<T> &y  );
  template<typename T>     
  void I2_X_D1(GTMatrix<T> &D1, GTVector<T> &u, GSIZET N1, GSIZET N2, GTVector<T> &y);
  template<typename T>     
  void D2_X_I1(GTMatrix<T> &D2T, GTVector<T> &u, GSIZET N1, GSIZET N2, GTVector<T> &y);
  template<typename T>     
  void I3_X_I2_X_D1(GTMatrix<T> &D1, GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET N3,
                    GTVector<T> &y);
  template<typename T>     
  void I3_X_D2_X_I1(GTMatrix<T> &D2T, GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET N3,
                    GTVector<T> &y);
  template<typename T>     
  void D3_X_I2_X_I1(GTMatrix<T> &D3T, GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET N3,
                    GTVector<T> &y);
  template<typename T>     
  void matvec_prod(GTVector<T> &vret, GTMatrix<T> &A, GTVector<T> &b);

  template<typename T>     
  void matmat_prod(GTMatrix<T> &C, GTMatrix<T> &A, GTMatrix<T> &B);

  template<typename T>  
  void    cross_prod_k(GTVector<GTVector<T>*> &x, GINT *iind, GINT nind, GINT isgn, 
                       GTVector<GTVector<T>*> &cross);
  template<typename T>  
  void    cross_prod  (GTVector<GTVector<T>*> &x, GTVector<GTVector<T>*> &y, 
                       GINT *iind, GINT nind, GTVector<GTVector<T>*> &cross);
  template<typename T>  
  void    normalize_euclidean(GTVector<GTVector<T>*> &x, GINT *iind, GINT nind);

  template<typename T>  
  void     add(GTVector<T> &vret, GTVector<T> &va, GTVector<T> &vb, T a, T b);

  template<typename T>  
  void    Dg2_X_D1           (GTMatrix<T> &D1, GTVector<T> &Dg2, GTVector<T> &x, GTVector<T> &tmp, GTVector<T> &y);

  template<typename T>  
  void    D2_X_Dg1           (GTVector<T> &Dg1, GTMatrix<T> &D2T, GTVector<T> &x, GTVector<T> &tmp, GTVector<T> &y);

  template<typename T>  
  void    Dg3_X_Dg2_X_D1     (GTMatrix<T> &D1 , GTVector<T> &Dg2, GTVector<T> &Dg3, GTVector<T> &x, GTVector<T> &tmp, GTVector<T> &y);

  template<typename T>  
  void    Dg3_X_D2_X_Dg1     (GTVector<T> &Dg1, GTMatrix<T> &D2T,  GTVector<T> &Dg3, GTVector<T> &x, GTVector<T> &tmp, GTVector<T> &y);

  template<typename T>  
  void    D3_X_Dg2_X_Dg1     (GTVector<T> &Dg1, GTVector<T> &Dg2, GTMatrix<T> &D3T, GTVector<T> &x, GTVector<T> &tmp, GTVector<T> &y);

  template<typename T>  
  void compute_grefderiv(GGrid &grid, GTVector<T> &u, GTVector<T> &etmp,
                         GINT idir, GBOOL dotrans, GTVector<T> &du);

  template<typename T>  
  void    compute_grefderivs(GGrid &grid, GTVector<T> &u, GTVector<T> &etmp,
                             GBOOL btrans, GTVector<GTVector<T>*> &du);
  template<typename T>  
  void    compute_grefderivsW(GGrid &grid, GTVector<T> &u, GTVector<T> &etmp,
                              GBOOL btrans, GTVector<GTVector<T>*> &du);
  template<typename T>  
  void    compute_grefdiv(GGrid &grid, GTVector<GTVector<T>*> &u, GTVector<T> &etmp,
                          GBOOL btrans, GTVector<T> &divu);
  template<typename T>  
  void    compute_grefdiviW(GGrid &grid, GTVector<GTVector<T>*> &u, GTVector<T> &etmp,
                           GBOOL btrans, GTVector<T> &divu);


};


#endif
