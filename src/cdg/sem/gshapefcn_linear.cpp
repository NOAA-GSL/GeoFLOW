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
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gshapefcn_linear.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GShapeFcn_linear::GShapeFcn_linear()
{
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GShapeFcn_linear::~GShapeFcn_linear()
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
void GShapeFcn_linear::Ni(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  GSIZET n=0;

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
void GShapeFcn_linear::Ni_1d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  GSIZET n=0;

  GTVector<GFTYPE>  *xi0 = xi[0];
  GFTYPE rm, rp;
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
//          N     : computation of coordinate at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_linear::Ni_2d(GTVector<GINT> &ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  GSIZET n=0;

  GTVector<GFTYPE>  *xi0 = xi[0];
  GTVector<GFTYPE>  *xi1 = xi[1];
  GFTYPE rm, rp, sm, sp;
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
// DESC   : Compute ishape shape function. Only 
//          Shape functions are bi-linear:
//              N_i = (1 +/- xi ) ( (1 +/- eta) (1 +/- zeta)
//          where +/- are taken depending on the position of the vertex, in 3d
// ARGS   : ishape: which shape function to compute. Is a 1-index array, with 
//                  ishape[0] going from 0 - 7 in 3d.
//          xi    : [*xi_0[N0], *xi_1[N1], ... *xi_2[N_GDIM]]: array of reference interval 
//                  points in each direction. Each array should have the number of 
//                  ref interfal points for that coord direction.
//          N     : computation of coordinate at each reference node point, 
//                  corresponding to direction icoord
// 
//             
// RETURNS:  none
//**********************************************************************************
void GShapeFcn_linear::Ni_3d(GTVector<GINT> & ishape, 
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N)
{
  GSIZET n=0;

  GTVector<GFTYPE>  *xi0 = xi[0];
  GTVector<GFTYPE>  *xi1 = xi[1];
  GTVector<GFTYPE>  *xi2 = xi[2];
  GFTYPE rm, rp, sm, sp, tm, tp;

  switch ( ishape[0] ) {
    case 0:
      for ( GSIZET k=0; k<xi2->size(); k++ ) {
        tm = 1.0 - (*xi2)[k]; tp = 1.0 + (*xi2)[k];
        for ( GSIZET j=0; j<xi1->size(); j++ ) {
          sm = 1.0 - (*xi1)[j]; sp = 1.0 + (*xi1)[j];
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
          for ( GSIZET i=0; i<xi0->size(); i++, n++ ) {
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
void GShapeFcn_linear::dNdXi(GTVector<GINT> &ishape, GINT j, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dxdxi)
{

#if defined(_G_IS1D)
  dNdXi_1d(ishape, j, xi, dxdxi);
#elif defined(_G_IS2D)
  dNdXi_2d(ishape, j, xi, dxdxi);
#elif defined(_G_IS3D)
  dNdXi_3d(ishape, j, xi, dxdxi);
#endif

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
void GShapeFcn_linear::dNdXi_1d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dxdxi)
{

  assert(jder==1 && "Invalid matrix element");

  GTVector<GFTYPE>  *xi0 = xi[0];

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
void GShapeFcn_linear::dNdXi_2d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dxdxi)
{

  assert(jder>=1 && jder<=2 && "Invalid matrix element");

  GSIZET n = 0;
  GTVector<GFTYPE>  *xi0 = xi[0];
  GTVector<GFTYPE>  *xi1 = xi[1];
  GFTYPE rm, rp, sm, sp;

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
void GShapeFcn_linear::dNdXi_3d(GTVector<GINT> &ishape, GINT jder, 
                               GTVector<GTVector<GFTYPE>*> &xi, 
                               GTVector<GFTYPE> &dxdxi)
{

  assert(jder>=1 && jder<=3 && "Invalid matrix element");

  GTVector<GFTYPE>  *xi0 = xi[0];
  GTVector<GFTYPE>  *xi1 = xi[1];
  GTVector<GFTYPE>  *xi2 = xi[2];
  GFTYPE rm, rp, sm, sp, tm, tp;
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

