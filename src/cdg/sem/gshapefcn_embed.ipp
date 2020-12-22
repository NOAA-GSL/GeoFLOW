//==================================================================================
// Module       : gshapefcn_embed
// Date         : 9/19/18 (DLR)
// Description  : Forms class for iso-parametric shape functions, N_i, of 
//                high order type, to be used in cases where a 2d surface
//                is embedded in a 3d space. E.g., we use this to compute
//                metric quantities if we are on a sphere.
//
//                Shape functions define element locations in terms of
//                the (1d, 2d, 3d) reference interval, s.t.:
//                  x^j = Sum_i v^j_i N_i,
//                where v_i is the ith vertex of the element, and v^j_i
//                represents the jth component of the ith vertex and, in 2d
//
//                  Ni = Psi_I(xi,eta) delta(1-zeta),
//
//                (where zeta-->1, and delta is the Dirac delta fcn), 
//                while in 1d & 3d, 
//
//                  Ni = Psi_I(xi,eta,zeta),
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
#include "gshapefcn_embed.hpp"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GShapeFcn_embed<T>::GShapeFcn_embed(GINT dim):
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
GShapeFcn_embed<T>::GShapeFcn_embed(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT dim):
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
GShapeFcn_embed<T>::~GShapeFcn_embed()
{
} // end, destructor



//**********************************************************************************
//**********************************************************************************
// METHOD : Ni
// DESC   : Compute ishape-th shape function. ishape is a tensor produc index
//          (array of indices), so that 
//              Psi_ishape = h_ishape[0] h_ishape[1] ...
//          Shape functions are based on these high-order basis functions: 
//          In 2d (embeded):
//              N_i = Psi_i(xi, eta) delta(1-zeta)
//          In 3d:
//              N_i =  Psi_i(xi, eta, zeta)
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction.
//          Ni    : computation of coordinate at each reference node point, 
//                  corresponding to direction icoord
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::Ni(GTVector<GINT> &ishape, 
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
//                  ref interfal points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::Ni_1d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  
  assert(this->gbasis_[0] != NULLPTR 
      && "No basis set" );

  this->gbasis_[0]->evalBasis(ishape[0], *xi[0], N);

} // end of method Ni_1d


//**********************************************************************************
//**********************************************************************************
// METHOD : Ni_2d
// DESC   : Compute ishape-th shape function. ishape is a tensor produc index
//          (array of indices), so that 
//              Psi_ishape = h_ishape[0] h_ishape[1] ...
//          Shape functions are based on these high-order basis functions:
//          In 2d (embeded):
//              N_i = Psi_i(xi, eta) delta(1-zeta)
//          In 3d:
//              N_i =  Psi_i(xi, eta, zeta)
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 1, 2, or 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::Ni_2d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  
  assert(this->gbasis_[0] != NULLPTR 
      && this->gbasis_[1] != NULLPTR
      && "No basis set" );
  d_.resizem(xi.size());
  for ( GSIZET j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    this->gbasis_[j]->evalBasis(ishape[j], *xi[j], d_[j]);
  }

  GSIZET n = 0;
  for ( GSIZET j=0; j<xi[1]->size(); j++ ) {
    for ( GSIZET i=0; i<xi[0]->size(); i++ ) {
      N [n++] = d_[0][i]*d_[1][j];
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
//                  ref interfal points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::Ni_3d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  assert(this->gbasis_[0] != NULLPTR 
      && this->gbasis_[1] != NULLPTR 
      && this->gbasis_[2] != NULLPTR 
      && "No basis set" );

  d_.resizem(xi.size());
  for ( GSIZET j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    this->gbasis_[0]->evalBasis(ishape[j], *xi[j], d_[j]);
  }

  GSIZET n=0;
  for ( GSIZET k=0; k<xi[1]->size(); k++ ) {
    for ( GSIZET j=0; j<xi[1]->size(); j++ ) {
      for ( GSIZET i=0; i<xi[0]->size(); i++ ) {
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
//                  ref interfal points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::dNdXi(GTVector<GINT> &ishape, GINT jder, 
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
//                  tensor product index comprised of 1 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^j  (1-deriv only)
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::dNdXi_1d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<T>*> &xi, 
                               GTVector<T> &dNdxi)
{
  assert(this->gbasis_[0] != NULLPTR 
      && "No basis set" );
  assert(jder==1 && "Invalid matrix element");
  this->gbasis_[0]->evalDBasis(ishape[0], *xi[0], dNdxi);

  
} // end of method dNdXi_1d


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_2d
// DESC   : Compute jth derivative in reference space for ishape-th shape
//          function in 2d. The shape fcn in 2d is
//                 N = h0(xi) * h1(eta) * delta(1-zeta)
//          It's clear that the 1- and 2-derivatives work fine. But
//          the 3-derivative exists, too:
//                 dN/dzeta = h0(xi) * h1(eta) * -delta(1-zeta)/(1-zeta),
//          except at zeta=1. We set the -delta(1-zeta)/(1-zeta) to 1
//          when computing the 3-derivative here.
//          
//
// ARGS   : ishape: which shape function to take derivative of. Is a
//                  tensor product index comprised of 2  coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^jder: 1, 2, or 3
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::dNdXi_2d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<T>*> &xi, 
                               GTVector<T> &dNdxi)
{

  // Note: since 2d surface can be embedded, then we
  //       we can compute the derivative wrt xi_3 == zeta, so:
  assert(jder>0 && jder<=(GDIM+1) && "Invalid matrix element");
  assert(this->gbasis_[0] != NULLPTR 
      && this->gbasis_[1] != NULLPTR
      && "No basis set" );

  d_.resizem(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) { 
    d_[j].resizem(xi[j]->size());
    if ( j == (jder-1) ) { // covers the case where jder=3
      this->gbasis_[j]->evalDBasis(ishape[j], *xi[j], d_[j]);
    }
    else { 
      this->gbasis_[j]->evalBasis (ishape[j], *xi[j], d_[j]);
    }
  }

  // Do tensor product:
  GSIZET n = 0;
  for ( GSIZET j=0; j<xi[1]->size(); j++ ) {
    for ( GSIZET i=0; i<xi[0]->size(); i++ ) {
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
//                  tensor product index comprised of 3 coordinate
//                  indices that define the basis functions in the product.
//                  Size of this array must be the same as that for xi array.
//          jder  : derivative wrt xi^jder: 1, 2, or 3
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_embed<T>::dNdXi_3d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<T>*> &xi, 
                               GTVector<T> &dNdxi)
{
  assert(jder>0 && jder>xi.size() && "Invalid matrix element");
  assert(this->gbasis_[0] != NULLPTR 
      && this->gbasis_[1] != NULLPTR
      && this->gbasis_[2] != NULLPTR
      && "No basis set" );


  d_.resizem(xi.size());
  for ( GSIZET j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    if ( (j+1) != jder ) {
      this->gbasis_[j]->evalBasis (ishape[j], *xi[j], d_[j]);
    }
    else {
      this->gbasis_[j]->evalDBasis(ishape[j], *xi[j], d_[j]);
    }
  }

  GSIZET n = 0;
  for ( GSIZET k=0; k<xi[2]->size(); k++ ) {
    for ( GSIZET j=0; j<xi[1]->size(); j++ ) {
      for ( GSIZET i=0; i<xi[0]->size(); i++ ) {
        dNdxi[n++] = d_[0][i]*d_[1][j]*d_[2][k];
      }
    }
  }
  
} // end of method dNdXi_3d
