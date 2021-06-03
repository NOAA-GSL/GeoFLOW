//==================================================================================
// Module       : gshapefcn_hostd
// Date         : 9/19/18 (DLR)
// Description  : Forms class for iso-parametric shape functions, N_i, of 
//                high order type, to be used in cases where a 2d surface
//                is not embedded in a 3d space. This is the shape function
//                use in the 'standard' way. Note: 'hostd' stands for 'high
//                order standard'.
//
//                Shape functions define element locations in terms of
//                the (1d, 2d, 3d) reference interval, s.t.:
//                  x^j = Sum_i v^j_i N_i,
//                where v_i is the ith vertex of the element, and v^j_i
//                represents the jth component of the ith vertex and, 
//
//                  Ni = Psi_I(xi,eta,zeta, ...),
//
//                and Psi_I is the I-th tensor product 2d GL or GLL basis
//                function, computed as:
//
//                  Psi_I = h_i(xi) h_j(eta) h_k(zeta) ....
//
//                where I = I(i,j,k) is the tensor product index computed
//                from the coordinate indices, i, j, k.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gshapefcn_hostd.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GShapeFcn_hostd<T>::GShapeFcn_hostd(GINT dim):
GShapeFcn_base<T>(dim)
{
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GShapeFcn_hostd<T>::GShapeFcn_hostd(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT dim):
GShapeFcn_base<T>(b, dim)
{
} // end of constructor method (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GShapeFcn_hostd<T>::~GShapeFcn_hostd()
{
} // end, destructor



//**********************************************************************************
//**********************************************************************************
// METHOD : Ni
// DESC   : Compute ishape-th shape function. ishape is a tensor produc index
//          (array of indices), so that 
//              Psi_ishape = h_ishape[0] h_ishape[1] ...
//          Shape functions are based on these high-order basis functions: 
//          In 1d:
//              N_i = Psi_i(xi)
//          In 2d:
//              N_i = Psi_i(xi, eta)
//          In 3d:
//              N_i =  Psi_i(xi, eta, zeta)
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction.
//          N     : computation of coordinate at each reference node point, 
//                  corresponding to direction icoord
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::Ni(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{

  switch ( this->dim_ ) {
    case 1:
      Ni_1d(ishape, xi, N);
      break;
    case 2:
      Ni_2d(ishape, xi, N);
      break;
    case 3:
      Ni_3d(ishape, xi, N);
      break;
    default:
      assert(FALSE);
  }

} // end of method Ni


//**********************************************************************************
//**********************************************************************************
// METHOD : Ni_1d
// DESC   : Compute ishape-th shape function. ishape is a tensor produc index
//          (array of indices), so that 
//              Psi_ishape = h_ishape[0] h_ishape[1] ...
//          Shape functions are based on these high-order basis functions:
//              N_i = Psi_i(xi)
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::Ni_1d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  
  this->gbasis_[0]->evalBasis(ishape[0], *xi[0], N);

} // end of method Ni_1d


//**********************************************************************************
//**********************************************************************************
// METHOD : Ni_2d
// DESC   : Compute ishape-th shape function. ishape is a tensor produc index
//          (array of indices), so that 
//              Psi_ishape = h_ishape[0] h_ishape[1] ...
//          Shape functions are based on these high-order basis functions:
//              N_i = Psi_i(xi, eta)
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::Ni_2d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  
  d_.resizem(xi.size());

  for ( auto j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    this->gbasis_[j]->evalBasis(ishape[j], *xi[j], d_[j]);
  }

  GSIZET n = 0;
  for ( auto j=0; j<xi[1]->size(); j++ ) {
    for ( auto i=0; i<xi[0]->size(); i++ ) {
      N[n++] = d_[0][i]*d_[1][j];
    }
  }

} // end of method Ni_2d


//**********************************************************************************
//**********************************************************************************
// METHOD : Ni_3d
// DESC   : Compute ishape-th shape function. ishape is a tensor produc index
//          (array of indices), so that 
//              Psi_ishape = h_ishape[0] h_ishape[1] ...
//          Shape functions are based on these high-order basis functions:
//              N_i = Psi_i(xi, eta, zeta)
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::Ni_3d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  d_.resizem(xi.size());

  for ( auto j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    this->gbasis_[0]->evalBasis(ishape[j], *xi[j], d_[j]);
  }

  GSIZET n=0;
  for ( auto k=0; k<xi[1]->size(); k++ ) {
    for ( auto j=0; j<xi[1]->size(); j++ ) {
      for ( auto i=0; i<xi[0]->size(); i++ ) {
        N[n++] = d_[0][i]*d_[1][j]*d_[2][k];
      }
    }
  }

} // end of method Ni



//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi
// DESC   : Compute jth derivative in reference space for ishape-th shape
//          function
//
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^j; may be 1, 2, or 3
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::dNdXi(GTVector<GINT> &ishape, GINT jder, 
                            GTVector<GTVector<T>*> &xi, 
                            GTVector<T> &dNdxi)
{

  switch ( this->dim_ ) {
    case 1:
      dNdXi_1d(ishape, jder, xi, dNdxi);
      break;
    case 2:
      dNdXi_2d(ishape, jder, xi, dNdxi);
      break;
    case 3:
      dNdXi_3d(ishape, jder, xi, dNdxi);
      break;
    default:
      assert(FALSE);
  } 

} // end of method dNdXi


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_1d
// DESC   : Compute jth derivative in reference space for ishape-th shape
//          function in 1d
//
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^j  (1-deriv only)
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::dNdXi_1d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<T>*> &xi, 
                               GTVector<T> &dNdxi)
{
  assert(jder==1 && "Invalid matrix element");
  this->gbasis_[0]->evalDBasis(ishape[0], *xi[0], dNdxi);

  
} // end of method dNdXi_1d


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_2d
// DESC   : Compute jth derivative in reference space for ishape-th shape
//          function in 2d. The shape fcn in 2d is
//                 N = h0(xi) * h1(eta)
//
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^jder: 1, or 2
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::dNdXi_2d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<T>*> &xi, 
                               GTVector<T> &dNdxi)
{
  assert(jder>0 && jder<=2 && "Invalid matrix element");

  d_.resizem(xi.size());

  // Compute 1d components of derivative of full basis:
  for ( auto j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    if ( j == (jder-1)) {
      this->gbasis_[j]->evalDBasis(ishape[j], *xi[j], d_[j]);
    }
    else {
      this->gbasis_[j]->evalBasis (ishape[j], *xi[j], d_[j]);
    }
  }

  // Compute tensor product function, dNdxi:
  GSIZET n = 0;
  for ( auto j=0; j<xi[1]->size(); j++ ) {
    for ( auto i=0; i<xi[0]->size(); i++ ) {
      dNdxi[n++] = d_[0][i]*d_[1][j];
    }
  }
 
} // end of method dNdXi_2d


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_3d
// DESC   : Compute jth derivative in reference space for ishape-th shape
//          function in 3d
//
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^jder: 1, 2, or 3
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interval points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_hostd<T>::dNdXi_3d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<T>*> &xi, 
                               GTVector<T> &dNdxi)
{
  assert(jder>0 && jder<=3 && "Invalid matrix element");

  d_.resizem(xi.size());

  // Compute 1d components of derivative of full basis:
  for ( auto j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    if ( j == (jder-1) ) {
      this->gbasis_[j]->evalDBasis(ishape[j], *xi[j], d_[j]);
    }
    else {
      this->gbasis_[j]->evalBasis (ishape[j], *xi[j], d_[j]);
    }
  }

  // Compute tensor product function, dNdxi:
  GSIZET n = 0;
  for ( auto k=0; k<xi[2]->size(); k++ ) {
    for ( auto j=0; j<xi[1]->size(); j++ ) {
      for ( auto i=0; i<xi[0]->size(); i++ ) {
        dNdxi[n++] = d_[0][i]*d_[1][j]*d_[2][k];
      }
    }
  }
  
} // end of method dNdXi_3d
