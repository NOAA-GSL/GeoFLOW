//==================================================================================
// Module       : gshapefcn_hostd
// Date         : 9/19/18 (DLR)
// Description  : Forms class for iso-parametric shape functions, N_i, of 
//                high order type, to be used in cases where a 2d surface
//                is not embedded in a 3d space. This is the shape function
//                use in the 'standard' way.
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
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
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
GShapeFcn_hostd::GShapeFcn_hostd()
{
  gbasis_.resize(GDIM);
  gbasis_ = NULLPTR;
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GShapeFcn_hostd::GShapeFcn_hostd(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)
{
  gbasis_.resize(GDIM);
  gbasis_ = b;
} // end of constructor method (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GShapeFcn_hostd::~GShapeFcn_hostd()
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
//                  ref interfal points for that coord direction.
//          N     : computation of coordinate at each reference node point, 
//                  corresponding to direction icoord
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_hostd::Ni(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{

#if defined(_G_IS1D)
  Ni_1d(ishape, xi, N);
#elif defined(_G_IS2D)
  Ni_2d(ishape, xi, N);
#elif defined(_G_IS3D)
  Ni_3d(ishape, xi, N);
#endif


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
void GShapeFcn_hostd::Ni_1d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  
  gbasis_[0]->evalBasis(ishape[0], *xi[0], N);

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
//                  ref interfal points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_hostd::Ni_2d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  
  for ( GSIZET j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    gbasis_[j]->evalBasis(ishape[j], *xi[j], d_[j]);
  }

  GSIZET n = 0;
  for ( GSIZET j=0; j<xi[1]->size(); j++ ) {
    for ( GSIZET i=0; i<xi[0]->size(); i++ ) {
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
//                  ref interfal points for that coord direction.
//          N     : computation of shape fcn at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_hostd::Ni_3d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  for ( GSIZET j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    gbasis_[0]->evalBasis(ishape[j], *xi[j], d_[j]);
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
void GShapeFcn_hostd::dNdXi(GTVector<GINT> &ishape, GINT jder, 
                            GTVector<GTVector<GFTYPE>*> &xi, 
                            GTVector<GFTYPE> &dNdxi)
{

#if defined(_G_IS1D)
  dNdXi_1d(ishape, jder, xi, dNdxi);
#elif defined(_G_IS2D)
  dNdXi_2d(ishape, jder, xi, dNdxi);
#elif defined(_G_IS3D)
  dNdXi_3d(ishape, jder, xi, dNdxi);
#endif

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
//                  ref interfal points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_hostd::dNdXi_1d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dNdxi)
{
  assert(jder==1 && "Invalid matrix element");
  gbasis_[0]->evalDBasis(ishape[0], *xi[0], dNdxi);

  
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
//                  ref interfal points for that coord direction
//          dNdxi : computation of derivative at each tensor prod node point
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_hostd::dNdXi_2d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dNdxi)
{
  assert(jder>0 && jder<=(GDIM+1) && "Invalid matrix element");

  for ( GSIZET j=0; j<GDIM+1; j++ ) { 
    d_[j].resize(xi[j]->size());
    if ( j == (jder-1)) {
      gbasis_[j]->evalDBasis(ishape[j], *xi[j], d_[j]);
    }
    else {
      gbasis_[j]->evalBasis (ishape[j], *xi[j], d_[j]);
    }
  }

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
//                  tensor product index comprised of 1, 2, or 3 coordinate
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
void GShapeFcn_hostd::dNdXi_3d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dNdxi)
{
  assert(jder>0 && jder>xi.size() && "Invalid matrix element");

  for ( GSIZET j=0; j<xi.size(); j++ ) {
    d_[j].resize(xi[j]->size());
    if ( j == jder ) {
      gbasis_[j]->evalDBasis(ishape[j], *xi[j], d_[j]);
    }
    else {
      gbasis_[j]->evalBasis (ishape[j], *xi[j], d_[j]);
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
