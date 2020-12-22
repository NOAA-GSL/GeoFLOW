//==================================================================================
// Module       : gshapefcn_linear
// Date         : 9/19/18 (DLR)
// Description  : Forms class for iso-parametric shape functions, N_i, of 
//                bi- or tri-linear type. This class is restrictive, as it
//                requires (2, 4, 8) shape functions for (1, 2, 3)d. This can
//                be generalized for elements of different type.
//
//                Shape functions define element locations in terms of
//                the (1d, 2d, 3d) reference interval, s.t.:
//                  x^j = Sum_i v^j_i N_i,
//                where v_i is the ith vertex of the element, and v^j_i
//                represents the jth component of the ith vertex. The N_i
//                should be provided in the order that the bounding vertices
//                of the element take (see gelement_base.hpp) so that this 
//                sum can be computed efficiently.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gshapefcn_linear.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GShapeFcn_linear<T>::GShapeFcn_linear(GINT dim):
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
GShapeFcn_linear<T>::GShapeFcn_linear(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT dim):
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
GShapeFcn_linear<T>::~GShapeFcn_linear()
{
} // end, destructor



//**********************************************************************************
//**********************************************************************************
// METHOD : Ni
// DESC   : Compute x_j = sum_i v^j_i N_i: to compute geometric position, given 
//          vertices v_i correspoinding to shape function N_i.
//          Shape functions are bi-linear:
//              N_i = (1 +/- xi ) ( (1 +/- eta)
//          where +/- are taken depending on the position of the vertex.
// ARGS   : ishape: which shape function to compute. Is a 1-index array, with 
//                  ishape[0] going from 0 - 1 in 1d , 0-3 in 2d, and 0-7 in 3d.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction.
//          Ni    : computation of coordinate at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::Ni(GTVector<GINT> &ishape, 
                             GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{

  switch ( this->dim_ ) {
    case 1: 
      Ni_1d(ishape, xi, N);
      break;
    case 2: 
      Ni_1d(ishape, xi, N);
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
// DESC   : Compute x_j = sum_i v^j_i N_i: to compute geometric position, given 
//          vertices v_i correspoinding to shape function N_i.
//          Shape functions are bi-linear:
//              N_i = (1 +/- xi ) ( (1 +/- eta)
//          where +/- are taken depending on the position of the vertex, for 1d
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
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::Ni_1d(GTVector<GINT> &ishape, 
                                GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  GSIZET n=0;

  GTVector<T>  *xi0 = xi[0];
  T rm, rp;
  if ( ishape[0] == 0 ) {
    for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
      rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
      N   [n++] = rm;
    }
  }
  else if ( ishape[0] == 1 ) {
    for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
      rp = 1.0 + (*xi0)[i];
      N   [n++] = rp;
    }
  }
} // end of method Ni_1d


//**********************************************************************************
//**********************************************************************************
// METHOD : Ni_2d
// DESC   : Compute x_j = sum_i v^j_i N_i: to compute geometric position, given 
//          vertices v_i correspoinding to shape function N_i.
//          Shape functions are bi-linear:
//              N_i = (1 +/- xi ) ( (1 +/- eta)
//          where +/- are taken depending on the position of the vertex, for 2d
// ARGS   : ishape: which shape function to compute. Is a 1-index array, with 
//                  ishape[0] going from 0 - 3 in 2d.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction.
//          N     : shape function computed, defined at all node points
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::Ni_2d(GTVector<GINT> &ishape, 
                                GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  GSIZET n=0;

  GTVector<T>  *xi0 = xi[0];
  GTVector<T>  *xi1 = xi[1];
  T rm, rp, sm, sp;
  switch ( ishape[0] ) {
    case 0:
      for ( GSIZET j=0; j<xi1->size(); j++ ) {
        sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
        for ( GSIZET i=0; i<xi0->size(); i++ ) {
          rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
          N   [n++] = rm*sm; 
        }
      }
      break;
    case 1:
      for ( GSIZET j=0; j<xi1->size(); j++ ) {
        sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
        for ( GSIZET i=0; i<xi0->size(); i++ ) {
          rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
          N   [n++] = rp*sm;
        }
      }
      break;
    case 2:
      for ( GSIZET j=0; j<xi1->size(); j++ ) {
        sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
        for ( GSIZET i=0; i<xi0->size(); i++ ) {
          rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
          N   [n++] = rp*sp;
        }
      }
      break;
    case 3:
      for ( GSIZET j=0; j<xi1->size(); j++ ) {
        sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
        for ( GSIZET i=0; i<xi0->size(); i++ ) {
          rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
          N   [n++] = rm*sp; 
        }
      }
      break;
  }

} // end of method Ni_2d


//**********************************************************************************
//**********************************************************************************
// METHOD : Ni_3d
// DESC   : Compute ishape shape function. 
//          Shape functions are bi-linear:
//              N_i = (1 +/- xi ) ( (1 +/- eta) (1 +/- zeta)
//          where +/- are taken depending on the position of the vertex, in 3d
// ARGS   : ishape: which shape function to compute. Is a 1-index array, with 
//                  ishape[0] going from 0 - 7 in 3d.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction.
//          N     : shape function computed, defined at all node points
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::Ni_3d(GTVector<GINT> & ishape, 
                                GTVector<GTVector<T>*> &xi, GTVector<T> &N)
{
  GSIZET n=0;

  GTVector<T>  *xi0 = xi[0];
  GTVector<T>  *xi1 = xi[1];
  GTVector<T>  *xi2 = xi[2];
  T rm, rp, sm, sp, tm, tp;

  switch ( ishape[0] ) {
    case 0:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rm*sm*tm; 
          }
        }
      }
      break;
    case 1:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rp*sm*tm;
          }
        }
      }
      break;
    case 2:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rp*sp*tm;
          }
        }
      }
      break;
    case 3:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rm*sp*tm;
          }
        }
      }
      break;
    case 4:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rm*sm*tp;
          }
        }
      }
      break;
    case 5:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rp*sm*tp;
          }
        }
      }
      break;
    case 6:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rp*sp*tp;
          }
        }
      }
      break;
    case 7:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            N   [n++] = rm*sp*tp; 
          }
        }
      }
      break;

  } // end, switch

} // end of method Ni



//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi
// DESC   : Compute jth derivative in reference space for ith coordinate computed 
//          from the sum over shape functions:
//               x^i = sum_k v^i_k N_k(xi):
//
// ARGS   : ishape: which shape function to compute deriv of. Is a 1-index array, with 
//                  ishape[0] going from 0 - 1 in 1d , 0-3 in 2d, and 0-7 in 3d.
//          jder  : derivative wrt xi_jder
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          out   : computation of derivative
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::dNdXi(GTVector<GINT> &ishape, GINT j, 
                                GTVector<GTVector<T>*> &xi, 
                                GTVector<T> &dxdxi)
{

  switch ( this->dim_ ) {
    case 1: 
      dNdXi_1d(ishape, j, xi, dxdxi);
      break;
    case 2: 
      dNdXi_2d(ishape, j, xi, dxdxi);
      break;
    case 3: 
      dNdXi_3d(ishape, j, xi, dxdxi);
      break;
    default:
      assert(FALSE);
  }

} // end of method dNdXi


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_1d
// DESC   : Compute jth derivative of shape function in reference space 
//          in 1d.
//
// ARGS   : ishape: which shape function to compute deriv of. Is a 1-index array, with 
//                  ishape[0] going from 0 - 1 in 1d.
//          jder  : derivative wrt xi^j 
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          dxdxi : computation of derivative
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::dNdXi_1d(GTVector<GINT> &ishape, GINT jder, 
                                   GTVector<GTVector<T>*> &xi, 
                                   GTVector<T> &dxdxi)
{

  assert(jder==1 && "Invalid matrix element");

  GTVector<T>  *xi0 = xi[0];

  switch (jder) {
    case 0:
      for ( GSIZET n=0; n<xi0->size(); n++ ) {
        dxdxi [n] = -1.0;
      }
      break;
    case 1:
      for ( GSIZET n=0; n<xi0->size(); n++ ) {
        dxdxi [n] = 1.0;
      }
      break;
   }
  
} // end of method dNdXi_1d


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_2d
// DESC   : Compute jth derivative of shape function in reference space 
//          in 2d.
//
// ARGS   : ishape: which shape function to compute. Is a 1-index array, with 
//                  ishape[0] going from 0 - 3 in 2d.
//          jder  : derivative wrt xi^jder
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          dxdxi : computation of derivative
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::dNdXi_2d(GTVector<GINT> &ishape, GINT jder, 
                                   GTVector<GTVector<T>*> &xi, 
                                   GTVector<T> &dxdxi)
{

  assert(jder>=1 && jder<=2 && "Invalid matrix element");

  GSIZET n = 0;
  GTVector<T>  *xi0 = xi[0];
  GTVector<T>  *xi1 = xi[1];
  T rm, rp, sm, sp;

  switch (ishape[0]) {
    case 0:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            dxdxi[n++] = -sm;
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
           for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            dxdxi[n++] = -rm ;
          }
        }
      }
      break;
    case 1:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            dxdxi[n++] = sm;
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
           for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            dxdxi[n++] = -rp;   
          }
        }
      }
      break;
    case 2:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            dxdxi[n++] = sp;
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
           for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            dxdxi[n++] = rp;
          }
        }
      }
      break;
    case 3:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++ ) {
            dxdxi[n++] = -sp;
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
           for ( GSIZET i=0; i<xi0->size(); i++ ) {
            rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
            dxdxi[n++] = rm;
          }
        }
      }
      break;
      
    default:
      assert(FALSE && "Unknown shape function");
     
    }
 
} // end of method dNdXi_2d


//**********************************************************************************
//**********************************************************************************
// METHOD : dNdXi_3d
// DESC   : Compute jth derivative of shape function in reference space 
//          in 3d.
//
// ARGS   : ishape: which shape function to compute. Is a 1-index array, with 
//                  ishape[0] going from 0 - 3 in 2d.
//          jder  : derivative wrt xi^jder
//          xi    : [*xi_0[N0], *xi_1[N1], *xi_2[N3]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction
//          dxdxi : computation of derivative
//             
// RETURNS:  none
//**********************************************************************************
template<typename T>
void GShapeFcn_linear<T>::dNdXi_3d(GTVector<GINT> &ishape, GINT jder, 
                                   GTVector<GTVector<T>*> &xi, 
                                   GTVector<T> &dxdxi)
{

  assert(jder>=1 && jder<=3 && "Invalid matrix element");

  GTVector<T>  *xi0 = xi[0];
  GTVector<T>  *xi1 = xi[1];
  GTVector<T>  *xi2 = xi[2];
  T rm, rp, sm, sp, tm, tp;
  GSIZET n;

  n = 0;
  switch (ishape[0]) {
    case 0:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] = -sm*tm;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rm*   tm;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rm*sm;
            }
          }
        }
      }
      break;
    case 1:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] =    sm*tm;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rp*   tm;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rp*sm;
            }
          }
        }
      }
      break;
    case 2:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] =    sp*tm;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rp*   tm;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rp*sp;
            }
          }
        }
      }
      break;
    case 3:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] = -   sp*tm;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rm*   tm;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rm*sp;
            }
          }
        }
      }
      break;
    case 4:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] = -   sm*tp;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rm*   tp;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rm*sm;
            }
          }
        }
      }
      break;
    case 5:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] =   sm*tp;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = -rp*   tp;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rp*sm;
            }
          }
        }
      }
      break;
    case 6:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] =    sp*tp;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rp*   tp;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rp*sp;
            }
          }
        }
      }
      break;
    case 7:
      if ( jder == 1 ) { // 1-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              dxdxi [n++] = -   sp*tp;
            }
          }
        }
      }
      else if ( jder == 2 ) { // 2-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) {
          tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rm*   tp;
            }
          }
        }
      }
      else if ( jder == 3 ) { // 3-derivative
        for ( GSIZET k=0; k<xi2->size(); k++ ) { 
          for ( GSIZET j=0; j<xi1->size(); j++ ) {
            sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
            for ( GSIZET i=0; i<xi0->size(); i++ ) {
              rm = 1.0 - (*xi0)[i]; rp = 1.0 + (*xi0)[i];
              dxdxi [n++] = rm*sp;
            }
          }
        }
      }
      break;
  }
  
} // end of method dNdXi_3d

