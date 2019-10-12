//==================================================================================
// Module       : gmtk.cpp
// Date         : 1/31/18 (DLR)
// Description  : Math TooKit: namespace of C routines for various
//                math functions
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gmtk.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "ggrid.hpp"

using namespace std;

//#if !defined(_GMTK_GLOBAL_DATA)
//  #define _GMTK_GLOBAL_DATA
  GINT szMatCache_ = _G_MAT_CACHE_SIZE;
  GINT szVecCache_ = _G_VEC_CACHE_SIZE;
//#endif


namespace GMTK 
{

//**********************************************************************************
//**********************************************************************************
// METHOD : curl 
// DESC   : Compute curl component, idir, of input vector field
//          
// ARGS   : grid : grid
//          u    : input vector field. Must have >= GDIM components.
//          idir : curl component to compute. Must be appropriate for 
//                 problem dimension.
//          tmp  : tmp vector; must be of at least length 2.
//          curlc: result
// RETURNS: none.
//**********************************************************************************
template<>
void curl(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &u, const GINT idir, 
          GTVector<GTVector<GFTYPE>*> &tmp, GTVector<GFTYPE> &curlc)
{

  assert(tmp.size() >= 2 && "Insufficient temp space");

  // Handle 1c cases in 2d or 3d:
  if  ( u.size() < 2 ) {
     curlc = 0.0; 
  }

  // Handle 2.5-d or 2d-3c case:
  else if ( GDIM == 2 && u.size() > GDIM && grid.gtype() != GE_2DEMBEDDED ) {
    switch (idir) {
      case 1:
        grid.deriv(*u[2], 2, *tmp[0], curlc);
        curlc *= -1.0;
        break;
      case 2:
        grid.deriv(*u[2], 1, *tmp[0], curlc);
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
  }

  // Handle 2d-2c regular types:
  else if ( GDIM == 2  && u.size() == 2 && grid.gtype() == GE_REGULAR ) {
    switch (idir) {
      case 1:
      case 2:
        curlc = 0.0;
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
  }

  // Handle 3d-3c or embedded cases:
  else if ( GDIM == 3 || grid.gtype() == GE_2DEMBEDDED ) {
    switch (idir) {
      case 1:
        grid.deriv(*u[1], 3, *tmp[0], curlc);
        grid.deriv(*u[2], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      case 2:
        grid.deriv(*u[2], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 3, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
    }
  }
  else {
    assert(FALSE && "Curl cannot be computed");
  }


  return;

} // end of method curl


//**********************************************************************************
//**********************************************************************************
// METHOD : grad
// DESC   : Compute gradient component, idir, of input vector field.
//          Note: Don't really need this, as it's just another way
//                to refer to the Cartesian 'deriv' method in GGrid
//          
// ARGS   : grid : grid
//          u    : input (scalar) field. 
//          idir : gradient component to compute. Must be appropriate for 
//                 problem dimension.
//          tmp  : tmp vector; must be of at least length 1.
//          gradc: result
// RETURNS: none.
//**********************************************************************************
template<>
void grad(GGrid &grid, GTVector<GFTYPE> &u, const GINT idir, 
          GTVector<GTVector<GFTYPE>*> &tmp, GTVector<GFTYPE> &gradc)
{
  assert ( idir >0 && idir <=3 && "Invalid compoment specified");

  grid.deriv(u, idir, *tmp[0], gradc);

} // end of method grad



//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_D1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template <>
void D2_X_D1<GFLOAT>(GTMatrix<GFLOAT> &D1, GTMatrix<GFLOAT>  &D2T, 
           GTVector<GFLOAT> &u, GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET   N11, N12, N21, N22;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N12*N21 && y.size() >= N11*N22) ) {
    cout << "GMTK::D2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D2_X_D1 u as: y = D1 U D2T, where U is u 
  // considered in matrix form:

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N21);

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  fmxm(tmp.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &N21, &szMatCache_);

  // y = D2_X_I1 * tmp == TMP D2T (in mat form):
  fmxm(y.data(), tmp.data(), &N11, &N12, D2T.data().data(), &N21, &N22, &szMatCache_);

} // end of method D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template <>
void D2_X_D1<GDOUBLE>(GTMatrix<GDOUBLE> &D1, GTMatrix<GDOUBLE>  &D2T, 
           GTVector<GDOUBLE> &u, GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET   N11, N12, N21, N22;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N12*N21 && y.size() >= N11*N22) ) {
    cout << "GMTK::D2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D2_X_D1 u as: y = D1 U D2T, where U is u 
  // considered in matrix form:

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N21);

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  dmxm(tmp.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &N21, &szMatCache_);

  // y = D2_X_I1 * tmp == TMP D2T (in mat form):
  dmxm(y.data(), tmp.data(), &N11, &N12, D2T.data().data(), &N21, &N22, &szMatCache_);

} // end of method D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_D1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_D1<GQUAD>(GTMatrix<GQUAD> &D1, GTMatrix<GQUAD>  &D2T, 
                     GTVector<GQUAD> &u, GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET   N11, N12, N21, N22;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N12*N21 && y.size() >= N11*N22) ) {
    cout << "GMTK::D2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D2_X_D1 u as: y = D1 U D2T, where U is u 
  // considered in matrix form:

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N21);

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  qmxm(tmp.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &N21, &szMatCache_);

  // y = D2_X_I1 * tmp == TMP D2T (in mat form):
  qmxm(y.data(), tmp.data(), &N11, &N12, D2T.data().data(), &N21, &N22, &szMatCache_);

} // end of method D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I2_X_D1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = I2 X D1 u
//          where I2 is the 2-direction's identity, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          u   : operand vector; of size N1 X N2
//          y   : return vector result; must be at least of size
//                N1 X N2
// RETURNS: none
//**********************************************************************************
template<>
void I2_X_D1<GFLOAT>(GTMatrix<GFLOAT> &D1, 
           GTVector<GFLOAT> &u, GSIZET N1, GSIZET N2, GTVector<GFLOAT> &y)
{
  GSIZET ND1, ND2;

  ND1 = D1.size(1);
  ND2 = D1.size(2);

  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N2 && y.size() >= N1*N2) ) {
    cout << "GMTK::I2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = I2_X_D1 u:
  fmxm(y.data(), D1.data().data(), &ND1, &ND2, u.data(), &N1, &N2, &szMatCache_);


} // end of method I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = I2 X D1 u
//          where I2 is the 2-direction's identity, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          u   : operand vector; of size N1 X N2
//          y   : return vector result; must be at least of size
//                N1 X N2
// RETURNS: none
//**********************************************************************************
template<>
void I2_X_D1<GDOUBLE>(GTMatrix<GDOUBLE> &D1, 
           GTVector<GDOUBLE> &u, GSIZET N1, GSIZET N2, GTVector<GDOUBLE> &y)
{
  GSIZET ND1, ND2;

  ND1 = D1.size(1);
  ND2 = D1.size(2);

  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N2 && y.size() >= N1*N2) ) {
    cout << "GMTK::I2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = I2_X_D1 u:
  dmxm(y.data(), D1.data().data(), &ND1, &ND2, u.data(), &N1, &N2, &szMatCache_);


} // end of method I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I2_X_D1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = I2 X D1 u
//          where I2 is the 2-direction's identity, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          u   : operand vector; of size N1 X N2
//          y   : return vector result; must be at least of size
//                N1 X N2
// RETURNS: none
//**********************************************************************************
template<>
void I2_X_D1<GQUAD>(GTMatrix<GQUAD> &D1, 
           GTVector<GQUAD> &u, GSIZET N1, GSIZET N2, GTVector<GQUAD> &y)
{
  GSIZET ND1, ND2;

  ND1 = D1.size(1);
  ND2 = D1.size(2);

  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N2 && y.size() >= N1*N2) ) {
    cout << "GMTK::I2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = I2_X_D1 u:
  qmxm(y.data(), D1.data().data(), &ND1, &ND2, u.data(), &N1, &N2, &szMatCache_);


} // end of method I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg2_X_D1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D2) X D1 u
//          where D2 is specified as a vector, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          Dg2 : diag(D2), as a vector
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void Dg2_X_D1<GFLOAT>(GTMatrix<GFLOAT> &D1, GTVector<GFLOAT> &Dg2, GTVector<GFLOAT> &u, 
                       GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET   N11, N12, N2;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N11*N2 && y.size() >= N11*N2) ) {
    cout << "GMTK::Dg_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N2);

  // Compute y = Dg2_X_D1 u as (Dg2 X I) (I X D1) U:

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  fmxm(tmp.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &N2, &szMatCache_);

  // y = Dg2_X_I1 * tmp == TMP diag(D2T) = TMP diag(D2)  (in mat form):
  fmxDm(y.data(), tmp.data(), &N11, &N2, Dg2.data(), &N2, &szMatCache_);

} // end of method Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D2) X D1 u
//          where D2 is specified as a vector, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          Dg2 : diag(D2), as a vector
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void Dg2_X_D1<GDOUBLE>(GTMatrix<GDOUBLE> &D1, GTVector<GDOUBLE> &Dg2, GTVector<GDOUBLE> &u, 
                       GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET   N11, N12, N2;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N11*N2 && y.size() >= N11*N2) ) {
    cout << "GMTK::Dg_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N2);

  // Compute y = Dg2_X_D1 u as (Dg2 X I) (I X D1) U:

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  dmxm(tmp.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &N2, &szMatCache_);

  // y = Dg2_X_I1 * tmp == TMP diag(D2T) = TMP diag(D2)  (in mat form):
  dmxDm(y.data(), tmp.data(), &N11, &N2, Dg2.data(), &N2, &szMatCache_);

} // end of method Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D2) X D1 u
//          where D2 is specified as a vector, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          Dg2 : diag(D2), as a vector
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void Dg2_X_D1<GQUAD>(GTMatrix<GQUAD> &D1, GTVector<GQUAD> &Dg2, GTVector<GQUAD> &u, 
                       GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET   N11, N12, N2;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N11*N2 && y.size() >= N11*N2) ) {
    cout << "GMTK::Dg_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N2);

  // Compute y = Dg2_X_D1 u as (Dg2 X I) (I X D1) U:

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  qmxm(tmp.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &N2, &szMatCache_);

  // y = Dg2_X_I1 * tmp == TMP diag(D2T) = TMP diag(D2)  (in mat form):
  qmxDm(y.data(), tmp.data(), &N11, &N2, Dg2.data(), &N2, &szMatCache_);

} // end of method Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_I1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X I1 u
//          where I1 is the 1-direction's identity, and D2
//          the deriv matrix in the 2-direction
// ARGS   : D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_I1<GFLOAT>(GTMatrix<GFLOAT> &D2T, 
           GTVector<GFLOAT> &u, GSIZET N1, GSIZET N2, GTVector<GFLOAT> &y)
{
  GSIZET N21, N22;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N2 && y.size() >= N1*N2) ) {
    cout << "GMTK::D2_X_I1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = I2_X_D1 u = u * D2T:
  fmxm(y.data(), u.data(), &N1, &N2, D2T.data().data(), &N21, &N22, &szMatCache_);


} // end of method D2_X_I1

//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_I1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X I1 u
//          where I1 is the 1-direction's identity, and D2
//          the deriv matrix in the 2-direction
// ARGS   : D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_I1<GDOUBLE>(GTMatrix<GDOUBLE> &D2T, 
           GTVector<GDOUBLE> &u, GSIZET N1, GSIZET N2, GTVector<GDOUBLE> &y)
{
  GSIZET N21, N22;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N2 && y.size() >= N1*N2) ) {
    cout << "GMTK::D2_X_I1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = I2_X_D1 u = u * D2T:
  dmxm(y.data(), u.data(), &N1, &N2, D2T.data().data(), &N21, &N22, &szMatCache_);


} // end of method D2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_I1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X I1 u
//          where I1 is the 1-direction's identity, and D2
//          the deriv matrix in the 2-direction
// ARGS   : D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_I1<GQUAD>(GTMatrix<GQUAD> &D2T, 
           GTVector<GQUAD> &u, GSIZET N1, GSIZET N2, GTVector<GQUAD> &y)
{
  GSIZET N21, N22;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N2 && y.size() >= N1*N2) ) {
    cout << "GMTK::D2_X_I1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = I2_X_D1 u = u * D2T:
  qmxm(y.data(), u.data(), &N1, &N2, D2T.data().data(), &N21, &N22, &szMatCache_);


} // end of method I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_Dg1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X diag(D1) u
//          where Dg  is the 1-direction's operator expressed as a vector, and D2
//          the dense matrix in the 2-direction
// ARGS   : Dg1 : diag(D1), as a vector
//          D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_Dg1<GFLOAT>(GTVector<GFLOAT> &Dg1, GTMatrix<GFLOAT> &D2T, GTVector<GFLOAT> &u, 
                      GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET   N1, N21, N22;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N21 && y.size() >= N1*N21) ) {
    cout << "GMTK::D2_X_Dg1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D2_X_Dg u as: y = (D2_X_I)(I_X_Dg) U,
  //  where U is u considered in matrix form:
  tmp.resizem(N1*N21);
  // y = D2_X_I1 * tmp == TMP D2T (in mat form):
  fmxm(tmp.data(), u.data(), &N1, &N21, D2T.data().data(), &N21, &N22, &szMatCache_);

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  fDmxm(y.data(), Dg1.data(), &N1, tmp.data(), &N1, &N21, &szMatCache_);

} // end of method D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_Dg1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X diag(D1) u
//          where Dg  is the 1-direction's operator expressed as a vector, and D2
//          the dense matrix in the 2-direction
// ARGS   : Dg1 : diag(D1), as a vector
//          D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_Dg1<GDOUBLE>(GTVector<GDOUBLE> &Dg1, GTMatrix<GDOUBLE> &D2T, GTVector<GDOUBLE> &u, 
                      GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET   N1, N21, N22;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N21 && y.size() >= N1*N21) ) {
    cout << "GMTK::D2_X_Dg1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D2_X_Dg u as: y = (D2_X_I)(I_X_Dg) U,
  //  where U is u considered in matrix form:
  tmp.resizem(N1*N21);
  // y = D2_X_I1 * tmp == TMP D2T (in mat form):
  dmxm(tmp.data(), u.data(), &N1, &N21, D2T.data().data(), &N21, &N22, &szMatCache_);

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  dDmxm(y.data(), Dg1.data(), &N1, tmp.data(), &N1, &N21, &szMatCache_);

} // end of method D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_Dg1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X diag(D1) u
//          where Dg  is the 1-direction's operator expressed as a vector, and D2
//          the dense matrix in the 2-direction
// ARGS   : Dg1 : diag(D1), as a vector
//          D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D2_X_Dg1<GQUAD>(GTVector<GQUAD> &Dg1, GTMatrix<GQUAD> &D2T, GTVector<GQUAD> &u, 
                      GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET   N1, N21, N22;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N1*N21 && y.size() >= N1*N21) ) {
    cout << "GMTK::D2_X_Dg1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D2_X_Dg u as: y = (D2_X_I)(I_X_Dg) U,
  //  where U is u considered in matrix form:
  tmp.resizem(N1*N21);
  // y = D2_X_I1 * tmp == TMP D2T (in mat form):
  qmxm(tmp.data(), u.data(), &N1, &N21, D2T.data().data(), &N21, &N22, &szMatCache_);

  // tmp = I2_X_D1 * u == D1 U (in mat form):
  qDmxm(y.data(), Dg1.data(), &N1, tmp.data(), &N1, &N21, &szMatCache_);

} // end of method D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_D2_X_D1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          D3T : transpose of 3-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) x D3.dm(2) = D1.size(2) x D2T.size(1) x D3T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2)  = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) x D3.size(1) = D1.size(1) x D2T.size(2) x D3T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_D2_X_D1<GFLOAT>(GTMatrix<GFLOAT> &D1, GTMatrix<GFLOAT>  &D2T, GTMatrix<GFLOAT> &D3T,
           GTVector<GFLOAT> &u, GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET   N11, N12, N21, N22, N31, N32;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N12*N21*N31 && y.size() >= N11*N22*N32) ) {
    cout << "GMTK::D2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D3_X_D2_X_D1 u as (D3 X I2 X I1)(I3 X D2 X I1)(I3 X I2 X D1) U: 


  // y = I3_X_I2_X_D1 * u == D1 U (in mat form):
  GSIZET nxy = N21*N31;

  // Resize tmp only if its current size is less than required:
  tmp.resizem(nxy*N32);
  fmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &nxy, &szMatCache_);

  // tmp = I3_X_D2_X_I1 y:
  for ( GSIZET k=0; k<N32; k++ ) { // do mxm op for each 'plane':
    fmxm(tmp.data()+k*N11*N22, y.data()+k*N11*N22, &N11, &N22, D2T.data().data(), &N21, &N22, &szMatCache_);
  }

  // y = D3 X I X I tmp:
  nxy = N11*N22;
  fmxm(y.data(), tmp.data(), &nxy, &N32, D3T.data().data(), &N31, &N32, &szMatCache_);

} // end of method D3_X_D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_D2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          D3T : transpose of 3-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) x D3.dm(2) = D1.size(2) x D2T.size(1) x D3T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) x D3.size(1) = D1.size(1) x D2T.size(2) x D3T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_D2_X_D1<GDOUBLE>(GTMatrix<GDOUBLE> &D1, GTMatrix<GDOUBLE>  &D2T, GTMatrix<GDOUBLE> &D3T,
           GTVector<GDOUBLE> &u, GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET   N11, N12, N21, N22, N31, N32;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N12*N21*N31 && y.size() >= N11*N22*N32) ) {
    cout << "GMTK::D2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D3_X_D2_X_D1 u as (D3 X I2 X I1)(I3 X D2 X I1)(I3 X I2 X D1) U: 


  // y = I3_X_I2_X_D1 * u == D1 U (in mat form):
  GSIZET nxy = N21*N31;

  // Resize tmp only if its current size is less than required:
  tmp.resizem(nxy*N32);
  dmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &nxy, &szMatCache_);

  // tmp = I3_X_D2_X_I1 y:
  for ( GSIZET k=0; k<N32; k++ ) { // do mxm op for each 'plane':
    dmxm(tmp.data()+k*N11*N22, y.data()+k*N11*N22, &N11, &N22, D2T.data().data(), &N21, &N22, &szMatCache_);
  }

  // y = D3 X I X I tmp:
  nxy = N11*N22;
  dmxm(y.data(), tmp.data(), &nxy, &N32, D3T.data().data(), &N31, &N32, &szMatCache_);

} // end of method D3_X_D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_D2_X_D1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          D3T : transpose of 3-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) x D3.dm(2) = D1.size(2) x D2T.size(1) x D3T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) x D3.size(1) = D1.size(1) x D2T.size(2) x D3T.size(2)
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_D2_X_D1<GQUAD>(GTMatrix<GQUAD> &D1, GTMatrix<GQUAD>  &D2T, GTMatrix<GQUAD> &D3T,
           GTVector<GQUAD> &u, GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET N11, N12, N21, N22, N31, N32;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  #if defined(_G_BOUNDS_CHK)
  if ( !(u.size() >= N12*N21*N31 && y.size() >= N11*N22*N32) ) {
    cout << "GMTK::D2_X_D1" << "incompatible size" << endl;
    exit(1);
  }
  #endif

  // Compute y = D3_X_D2_X_D1 u as (D3 X I2 X I1)(I3 X D2 X I1)(I3 X I2 X D1) U: 


  // y = I3_X_I2_X_D1 * u == D1 U (in mat form):
  GSIZET nxy = N21*N31;

  // Resize tmp only if its current size is less than required:
  tmp.resizem(nxy*N32);
  qmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N12, &nxy, &szMatCache_);

  // tmp = I3_X_D2_X_I1 y:
  for ( GSIZET k=0; k<N32; k++ ) { // do mxm op for each 'plane':
    qmxm(tmp.data()+k*N11*N22, y.data()+k*N11*N22, &N11, &N22, D2T.data().data(), &N21, &N22, &szMatCache_);
  }

  // y = D3 X I X I tmp:
  nxy = N11*N22;
  qmxm(y.data(), tmp.data(), &nxy, &N32, D3T.data().data(), &N31, &N32, &szMatCache_);

} // end of method D3_X_D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_I2_X_D1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X I2 X D1 u
// ARGS   : D1      : 1-direction (dense) operator 
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void I3_X_I2_X_D1<GFLOAT>(GTMatrix<GFLOAT> &D1, GTVector<GFLOAT> &u,
                           GSIZET N1, GSIZET N2, GSIZET N3,
                           GTVector<GFLOAT> &y)
{
  GSIZET  N11, N12, NYZ, NN;

  N11 = D1.size(1);
  N12 = D1.size(2);
  NYZ = N2*N3;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::I3_X_I2_X_D1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

  fmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N1, &NYZ, &szMatCache_);


} // end of method I3_X_I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_I2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X I2 X D1 u
// ARGS   : D1      : 1-direction (dense) operator 
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void I3_X_I2_X_D1<GDOUBLE>(GTMatrix<GDOUBLE> &D1, GTVector<GDOUBLE> &u,
                           GSIZET N1, GSIZET N2, GSIZET N3,
                           GTVector<GDOUBLE> &y)
{
  GSIZET  N11, N12, NYZ, NN;

  N11 = D1.size(1);
  N12 = D1.size(2);
  NYZ = N2*N3;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::I3_X_I2_X_D1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

  dmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N1, &NYZ, &szMatCache_);


} // end of method I3_X_I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_I2_X_D1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X I2 X D1 u
// ARGS   : D1      : 1-direction (dense) operator 
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void I3_X_I2_X_D1<GQUAD>(GTMatrix<GQUAD> &D1, GTVector<GQUAD> &u,
                           GSIZET N1, GSIZET N2, GSIZET N3,
                           GTVector<GQUAD> &y)
{
  GSIZET  N11, N12, NYZ, NN;

  N11 = D1.size(1);
  N12 = D1.size(2);
  NYZ = N2*N3;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::I3_X_I2_X_D1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

  qmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N1, &NYZ, &szMatCache_);


} // end of method I3_X_I2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_D2_X_I1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X D2 X I1 u
// ARGS   : D2T     : 2-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void I3_X_D2_X_I1<GFLOAT>(GTMatrix<GFLOAT> &D2T, GTVector<GFLOAT> &u, 
                           GSIZET N1, GSIZET N2, GSIZET N3,
                           GTVector<GFLOAT> &y)
{
  GSIZET  N21, N22, NXY, NN;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::I3_X_D2_X_I1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

 for ( GSIZET k=0; k<N3; k++ ) {
    fmxm(y.data()+k*NXY, u.data()+k*NXY, &N1, &N2, D2T.data().data(), 
         &N21, &N22, &szMatCache_);
  }

} // end of method I3_X_D2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_D2_X_I1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X D2 X I1 u
// ARGS   : D2T     : 2-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void I3_X_D2_X_I1<GDOUBLE>(GTMatrix<GDOUBLE> &D2T, GTVector<GDOUBLE> &u, 
                           GSIZET N1, GSIZET N2, GSIZET N3,
                           GTVector<GDOUBLE> &y)
{
  GSIZET  N21, N22, NXY, NN;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::I3_X_D2_X_I1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

 for ( GSIZET k=0; k<N3; k++ ) {
    dmxm(y.data()+k*NXY, u.data()+k*NXY, &N1, &N2, D2T.data().data(), 
         &N21, &N22, &szMatCache_);
  }

} // end of method I3_X_D2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_D2_X_I1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X D2 X I1 u
// ARGS   : D2T     : 2-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void I3_X_D2_X_I1<GQUAD>(GTMatrix<GQUAD> &D2T, GTVector<GQUAD> &u, 
                           GSIZET N1, GSIZET N2, GSIZET N3,
                           GTVector<GQUAD> &y)
{
  GSIZET  N21, N22, NXY, NN;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::I3_X_D2_X_I1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

 for ( GSIZET k=0; k<N3; k++ ) {
    qmxm(y.data()+k*NXY, u.data()+k*NXY, &N1, &N2, D2T.data().data(), 
         &N21, &N22, &szMatCache_);
  }

} // end of method I3_X_D2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_I2_X_I1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X I2 X I1 u
// ARGS   : D3T     : 3-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_I2_X_I1<GFLOAT>(GTMatrix<GFLOAT> &D3T, GTVector<GFLOAT> &u,
                           GSIZET N1, GSIZET N2, GSIZET N3, 
                           GTVector<GFLOAT> &y)
{
  GSIZET  N31, N32, NXY, NN;

  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::D3_X_I2_X_I1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

 fmxm(y.data(), u.data(), &NXY, &N3, D3T.data().data(), 
      &N31, &N32, &szMatCache_);

} // end of method D3_X_I2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_I2_X_I1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X I2 X I1 u
// ARGS   : D3T     : 3-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_I2_X_I1<GDOUBLE>(GTMatrix<GDOUBLE> &D3T, GTVector<GDOUBLE> &u,
                           GSIZET N1, GSIZET N2, GSIZET N3, 
                           GTVector<GDOUBLE> &y)
{
  GSIZET  N31, N32, NXY, NN;

  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::D3_X_I2_X_I1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

 dmxm(y.data(), u.data(), &NXY, &N3, D3T.data().data(), 
      &N31, &N32, &szMatCache_);

} // end of method D3_X_I2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_I2_X_I1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X I2 X I1 u
// ARGS   : D3T     : 3-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimentions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_I2_X_I1<GQUAD>(GTMatrix<GQUAD> &D3T, GTVector<GQUAD> &u,
                           GSIZET N1, GSIZET N2, GSIZET N3, 
                           GTVector<GQUAD> &y)
{
  GSIZET  N31, N32, NXY, NN;

  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;

#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN || y.size() < NN ) {
    cout << "GMTK::D3_X_I2_X_I1: incompatible dimensions" << endl;
    exit(1);
  }
#endif

 qmxm(y.data(), u.data(), &NXY, &N3, D3T.data().data(), 
      &N31, &N32, &szMatCache_);

} // end of method D3_X_I2_X_I1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_Dg2_X_D1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X diag(D2) X D1 u
// ARGS   : D1   : GTMatrix of 1-operator
//          Dg1  : diag(D2)
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void Dg3_X_Dg2_X_D1<GFLOAT>(GTMatrix<GFLOAT> &D1, GTVector<GFLOAT> &Dg2, GTVector<GFLOAT> &Dg3,
                           GTVector<GFLOAT> &u, GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET N11, N12, N2, N3, NN, NXY, NYZ;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  N3  = Dg3.size();
  NXY = N11*N2;
  NYZ = N2*N3;
  NN  = N11*N2*N3;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"Dg3_X_D2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N3);

  // tmp = I X I X D1 u:
  fmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N11, &NYZ, &szMatCache_);

  // tmp1 = I X Diag(D2) X I tmp:
  for ( GSIZET k=0; k<N3; k++ ) {
    fmxDm(tmp.data()+k*NXY,  y.data()+k*NXY, &N11, &N12, Dg2.data(), &N2, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  fmxDm(y.data(), tmp.data(), &NXY, &N3, Dg3.data(), &N3, &szMatCache_);


} // end, method Dg3_X_Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_Dg2_X_D1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X diag(D2) X D1 u
// ARGS   : D1   : GTMatrix of 1-operator
//          Dg1  : diag(D2)
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void Dg3_X_Dg2_X_D1<GDOUBLE>(GTMatrix<GDOUBLE> &D1, GTVector<GDOUBLE> &Dg2, GTVector<GDOUBLE> &Dg3,
                           GTVector<GDOUBLE> &u, GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET N11, N12, N2, N3, NN, NXY, NYZ;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  N3  = Dg3.size();
  NXY = N11*N2;
  NYZ = N2*N3;
  NN  = N11*N2*N3;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"Dg3_X_D2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N3);

  // tmp = I X I X D1 u:
  dmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N11, &NYZ, &szMatCache_);

  // tmp1 = I X Diag(D2) X I tmp:
  for ( GSIZET k=0; k<N3; k++ ) {
    dmxDm(tmp.data()+k*NXY,  y.data()+k*NXY, &N11, &N12, Dg2.data(), &N2, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  dmxDm(y.data(), tmp.data(), &NXY, &N3, Dg3.data(), &N3, &szMatCache_);


} // end, method Dg3_X_Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_Dg2_X_D1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X diag(D2) X D1 u
// ARGS   : D1   : GTMatrix of 1-operator
//          Dg1  : diag(D2)
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void Dg3_X_Dg2_X_D1<GQUAD>(GTMatrix<GQUAD> &D1, GTVector<GQUAD> &Dg2, GTVector<GQUAD> &Dg3,
                           GTVector<GQUAD> &u, GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET N11, N12, N2, N3, NN, NXY, NYZ;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  N3  = Dg3.size();
  NXY = N11*N2;
  NYZ = N2*N3;
  NN  = N11*N2*N3;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"Dg3_X_D2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N3);

  // tmp = I X I X D1 u:
  qmxm(y.data(), D1.data().data(), &N11, &N12, u.data(), &N11, &NYZ, &szMatCache_);

  // tmp1 = I X Diag(D2) X I tmp:
  for ( GSIZET k=0; k<N3; k++ ) {
    qmxDm(tmp.data()+k*NXY,  y.data()+k*NXY, &N11, &N12, Dg2.data(), &N2, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  qmxDm(y.data(), tmp.data(), &NXY, &N3, Dg3.data(), &N3, &szMatCache_);


} // end, method Dg3_X_Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_D2_X_Dg1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X D2 X diag(D1) u
// ARGS   : Dg1  : transpose of dense matrix, D1
//          D2T  : GTMatrix for 2-operator, transpose
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void Dg3_X_D2_X_Dg1<GFLOAT>(GTVector<GFLOAT> &Dg1, GTMatrix<GFLOAT> &D2T, GTVector<GFLOAT> &Dg3,
                           GTVector<GFLOAT> &u, GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET N1, N21, N22, N3, NN, NXY, NYZ;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N3  = Dg3.size();
  NXY = N1*N21;
  NYZ = N21*N3;
  NN  = N1*N21*N3;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"Dg3_X_D2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N3);

  // tmp = I X I X Diag(D1) u:
  fDmxm(y.data(), Dg1.data(), &N1, u.data(), &N1, &NYZ, &szMatCache_);

  // tmp1 = I X D2 X I tmp:
  for ( GSIZET k=0; k<N3; k++ ) {
    fmxm(tmp.data()+k*NXY, y.data()+k*NXY, &N1, &N21, D2T.data().data(), &N21, &N22, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  fmxDm(y.data(),  tmp.data(), &NXY, &N3, Dg3.data(), &N3, &szMatCache_);


} // end, method Dg3_X_D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_D2_X_Dg1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X D2 X diag(D1) u
// ARGS   : Dg1  : transpose of dense matrix, D1
//          D2T  : GTMatrix for 2-operator, transpose
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void Dg3_X_D2_X_Dg1<GDOUBLE>(GTVector<GDOUBLE> &Dg1, GTMatrix<GDOUBLE> &D2T, GTVector<GDOUBLE> &Dg3,
                           GTVector<GDOUBLE> &u, GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET N1, N21, N22, N3, NN, NXY, NYZ;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N3  = Dg3.size();
  NXY = N1*N21;
  NYZ = N21*N3;
  NN  = N1*N21*N3;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"Dg3_X_D2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N3);

  // tmp = I X I X Diag(D1) u:
  dDmxm(y.data(), Dg1.data(), &N1, u.data(), &N1, &NYZ, &szMatCache_);

  // tmp1 = I X D2 X I tmp:
  for ( GSIZET k=0; k<N3; k++ ) {
    dmxm(tmp.data()+k*NXY, y.data()+k*NXY, &N1, &N21, D2T.data().data(), &N21, &N22, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  dmxDm(y.data(),  tmp.data(), &NXY, &N3, Dg3.data(), &N3, &szMatCache_);


} // end, method Dg3_X_D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_D2_X_Dg1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X D2 X diag(D1) u
// ARGS   : Dg1  : transpose of dense matrix, D1
//          D2T  : GTMatrix for 2-operator, transpose
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void Dg3_X_D2_X_Dg1<GQUAD>(GTVector<GQUAD> &Dg1, GTMatrix<GQUAD> &D2T, GTVector<GQUAD> &Dg3,
                           GTVector<GQUAD> &u, GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET N1, N21, N22, N3, NN, NXY, NYZ;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N3  = Dg3.size();
  NXY = N1*N21;
  NYZ = N21*N3;
  NN  = N1*N21*N3;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"Dg3_X_D2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N3);

  // tmp = I X I X Diag(D1) u:
  qDmxm(y.data(), Dg1.data(), &N1, u.data(), &N1, &NYZ, &szMatCache_);

  // tmp1 = I X D2 X I tmp:
  for ( GSIZET k=0; k<N3; k++ ) {
    qmxm(tmp.data()+k*NXY, y.data()+k*NXY, &N1, &N21, D2T.data().data(), &N21, &N22, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  qmxDm(y.data(),  tmp.data(), &NXY, &N3, Dg3.data(), &N3, &szMatCache_);


} // end, method Dg3_X_D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_Dg2_X_Dg1 (float)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X diag(D2) X diag(D1) u
// ARGS   : Dg1  : diag(D1)
//          Dg1  : diag(D2)
//          D3T  : GMatrix for 3-operator, transpose
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_Dg2_X_Dg1<GFLOAT>(GTVector<GFLOAT> &Dg1, GTVector<GFLOAT> &Dg2, GTMatrix<GFLOAT> &D3T,
                           GTVector<GFLOAT> &u, GTVector<GFLOAT> &tmp, GTVector<GFLOAT> &y)
{
  GSIZET N1, N2, N31, N32, NN, NXY, NYZ;

  N1  = Dg1.size();
  N2  = Dg2.size();
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NYZ = N2*N31;
  NN  = N1*N2*N31;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"D3_X_Dg2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N31);

  // tmp = I X I X Diag(D1) u:
  fDmxm(y.data(), Dg1.data(), &N1, u.data(), &N1, &NYZ, &szMatCache_);

  // tmp1 = I X D2 X I tmp:
  for ( GSIZET k=0; k<N31; k++ ) {
    fmxDm(tmp.data()+k*NXY, y.data()+k*NXY, &N1, &N2, Dg2.data(), &N2, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  fmxm(y.data(), tmp.data(), &NXY, &N31, D3T.data().data(), &N31, &N32, &szMatCache_);


} // end, method D3_X_Dg2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_Dg2_X_Dg1 (double)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X diag(D2) X diag(D1) u
// ARGS   : Dg1  : diag(D1)
//          Dg1  : diag(D2)
//          D3T  : GMatrix for 3-operator, transpose
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_Dg2_X_Dg1<GDOUBLE>(GTVector<GDOUBLE> &Dg1, GTVector<GDOUBLE> &Dg2, GTMatrix<GDOUBLE> &D3T,
                           GTVector<GDOUBLE> &u, GTVector<GDOUBLE> &tmp, GTVector<GDOUBLE> &y)
{
  GSIZET N1, N2, N31, N32, NN, NXY, NYZ;

  N1  = Dg1.size();
  N2  = Dg2.size();
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NYZ = N2*N31;
  NN  = N1*N2*N31;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"D3_X_Dg2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N31);

  // tmp = I X I X Diag(D1) u:
  dDmxm(y.data(), Dg1.data(), &N1, u.data(), &N1, &NYZ, &szMatCache_);

  // tmp1 = I X D2 X I tmp:
  for ( GSIZET k=0; k<N31; k++ ) {
    dmxDm(tmp.data()+k*NXY, y.data()+k*NXY, &N1, &N2, Dg2.data(), &N2, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  dmxm(y.data(), tmp.data(), &NXY, &N31, D3T.data().data(), &N31, &N32, &szMatCache_);


} // end, method D3_X_Dg2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_Dg2_X_Dg1 (quad)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X diag(D2) X diag(D1) u
// ARGS   : Dg1  : diag(D1)
//          Dg1  : diag(D2)
//          D3T  : GMatrix for 3-operator, transpose
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<>
void D3_X_Dg2_X_Dg1<GQUAD>(GTVector<GQUAD> &Dg1, GTVector<GQUAD> &Dg2, GTMatrix<GQUAD> &D3T,
                           GTVector<GQUAD> &u, GTVector<GQUAD> &tmp, GTVector<GQUAD> &y)
{
  GSIZET N1, N2, N31, N32, NN, NXY, NYZ;

  N1  = Dg1.size();
  N2  = Dg2.size();
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NYZ = N2*N31;
  NN  = N1*N2*N31;
#if defined(GARRAY_BOUNDS)
  if ( u.size() < NN  || y.size() < NN ) {
    cout <<"D3_X_Dg2_X_Dg1: incompatible vectors" << endl;
    exit(1);
  }
#endif

  tmp.resizem(NXY*N31);

  // tmp = I X I X Diag(D1) u:
  qDmxm(y.data(), Dg1.data(), &N1, u.data(), &N1, &NYZ, &szMatCache_);

  // tmp1 = I X D2 X I tmp:
  for ( GSIZET k=0; k<N31; k++ ) {
    qmxDm(tmp.data()+k*NXY, y.data()+k*NXY, &N1, &N2, Dg2.data(), &N2, &szMatCache_);
  }

  // y = Dg3 X I X I tmp1:
  qmxm(y.data(), tmp.data(), &NXY, &N31, D3T.data().data(), &N31, &N32, &szMatCache_);


} // end, method D3_X_Dg2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : add (float)
// DESC   : point-by-point addition, returned in specified GTVector:
//            vret = a*va + b*vb, for GFLOAT types
// ARGS   : vret: GTVector<T> return vector
//          va  : first vector operatnd
//          vb  : second vector operatnd
//          a,b: factors multiplying va and vb, respectively
// RETURNS: GTVector & 
//**********************************************************************************
#pragma acc routine vector
template<>
void add<GFLOAT>(GTVector<GFLOAT> &vret, const GTVector<GFLOAT> &va, const GTVector<GFLOAT> &vb, GFLOAT a, GFLOAT b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( va.size() < vret.size() || vb.size() < vret.size() ) {
    cout << "GTVector<T>::add: " << "incompatible size" << endl;
while(1){};
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
  for ( GLLONG j=vret.getIndex().beg(); j<=vret.getIndex().end(); j+=vret.getIndex().stride() ) {
    vret[j] = a*va[j] + b*vb[j];
  }
  #else
  GSIZET nn = vret.getIndex().end() - vret.getIndex().beg() + 1;
  fzaxpby(vret.data(), const_cast<GFLOAT*>(va.data()), &a, 
                       const_cast<GFLOAT*>(vb.data()), &b, &nn, &szVecCache_);
  #endif

} // end, add 


//**********************************************************************************
//**********************************************************************************
// METHOD : add (double)
// DESC   : point-by-point addition, returned in specified GTVector:
//            vret = a*va + b*vb, for GFLOAT types
// ARGS   : vret: GTVector<T> return vector
//          va  : first vector operatnd
//          vb  : second vector operatnd
//          a,b: factors multiplying va and vb, respectively
// RETURNS: none
//**********************************************************************************
#pragma acc routine vector
template<>
void add<GDOUBLE>(GTVector<GDOUBLE> &vret, const GTVector<GDOUBLE> &va, const GTVector<GDOUBLE> &vb, GDOUBLE a, GDOUBLE b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( va.size() < vret.size() || vb.size() < vret.size() ) {
    cout << "GTVector<T>::add: " << "incompatible size" << endl;
while(1){};
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
  for ( GLLONG j=vret.getIndex().beg(); j<=vret.getIndex().end(); j+=vret.getIndex().stride() ) {
    vret[j] = a*va[j] + b*vb[j];
  }
  #else
  GSIZET nn = vret.getIndex().end() - vret.getIndex().beg() + 1;
  dzaxpby(vret.data(), const_cast<GDOUBLE*>(va.data()), &a, 
                       const_cast<GDOUBLE*>(vb.data()), &b, &nn, &szVecCache_);
  #endif

} // end, add 


//**********************************************************************************
//**********************************************************************************
// METHOD : add (quad)
// DESC   : point-by-point addition, returned in specified GTVector:
//            vret = a*va + b*vb, for GFLOAT types
// ARGS   : vret: GTVector<T> return vector
//          va  : first vector operatnd
//          vb  : second vector operatnd
//          a,b: factors multiplying va and vb, respectively
// RETURNS: none
//**********************************************************************************
#pragma acc routine vector
template<>
void add<GQUAD>(GTVector<GQUAD> &vret, const GTVector<GQUAD> &va, const GTVector<GQUAD> &vb, GQUAD a, GQUAD b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( va.size() < vret.size() || vb.size() < vret.size() ) {
    cout << "GTVector<T>::add: " << "incompatible size" << endl;
while(1){};
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
  for ( GLLONG j=vret.getIndex().beg(); j<=vret.getIndex().end(); j+=vret.getIndex().stride() ) {
    vret[j] = a*va[j] + b*vb[j];
  }
  #else
  GSIZET nn = vret.getIndex().end() - vret.getIndex().beg() + 1;
  qzaxpby(vret.data(), const_cast<GQUAD*>(va.data()), &a, 
                       const_cast<GQUAD*>(vb.data()), &b, &nn, &szVecCache_);
  #endif

} // end, add 



//**********************************************************************************
//**********************************************************************************
// METHOD : operator * mat-vec (float)
// DESC   : matrix-vector product, returns product,
//          without destroying *this data, for GFLOAT type
// ARGS   : GTVector &
// RETURNS: none
//**********************************************************************************
template<>
void matvec_prod<GFLOAT>(GTVector<GFLOAT> &vret, const GTMatrix<GFLOAT> &A, const GTVector<GFLOAT> &b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( b.size() < A.size(2) ) {
    cout << "GMTK::matvec_prod: " << "incompatible size" << endl;
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
   for ( GSIZET i=0; i<n1_; i++ ) {
     vret[i] = 0;
     for ( GSIZET j=0; j<n2_; j++ ) {
       vret[i] += A(i,j) * b(j);
     }
   }
  #else
  GSIZET n1 = A.size(1);
  GSIZET n2 = A.size(2);
  fmxv(vret.data(), const_cast<GFLOAT*>(A.data().data()), 
                    const_cast<GFLOAT*>(b.data())       , 
                    &n1, &n2, &szMatCache_);
  #endif


} // end of operator * mat-vec, GFLOAT


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * mat-vec (double)
// DESC   : matrix-vector product, returns product,
//          without destroying *this data, for GDOUBLE type
// ARGS   : GTVector &
// RETURNS: none
//**********************************************************************************
template<>
void matvec_prod<GDOUBLE>(GTVector<GDOUBLE> &vret, const GTMatrix<GDOUBLE> &A, const GTVector<GDOUBLE> &b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( b.size() < A.size(2) ) {
    cout << "GMTK::matvec_prod: " << "incompatible size" << endl;
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
   for ( GSIZET i=0; i<n1_; i++ ) {
     vret[i] = 0;
     for ( GSIZET j=0; j<n2_; j++ ) {
       vret[i] += A(i,j) * b(j);
     }
   }
  #else
  GSIZET n1 = A.size(1);
  GSIZET n2 = A.size(2);
  dmxv(vret.data(), const_cast<GDOUBLE*>(A.data().data()), 
                    const_cast<GDOUBLE*>(b.data())       , 
                    &n1, &n2, &szMatCache_);
  #endif


} // end of operator * mat-vec, GDOUBLE


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * mat-vec (quad)
// DESC   : matrix-vector product, returns product,
//          without destroying *this data, for GQUAD type
// ARGS   : GTVector &
// RETURNS: none
//**********************************************************************************
template<>
void matvec_prod<GQUAD>(GTVector<GQUAD> &vret, const GTMatrix<GQUAD> &A, const GTVector<GQUAD> &b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( b.size() < A.size(2) ) {
    cout << "GMTK::matvec_prod: " << "incompatible size" << endl;
    exit(1);
  }
  #endif

  #if !defined(_G_USE_GBLAS)
   for ( GSIZET i=0; i<n1_; i++ ) {
     vret[i] = 0;
     for ( GSIZET j=0; j<n2_; j++ ) {
       vret[i] += A(i,j) * b(j);
     }
   }
  #else
  GSIZET n1 = A.size(1);
  GSIZET n2 = A.size(2);
  qmxv(vret.data(), const_cast<GQUAD*>(A.data().data()), 
                    const_cast<GQUAD*>(b.data())       , 
                    &n1, &n2, &szMatCache_);
  #endif


} // end of operator * mat-vec, GQUAD



//**********************************************************************************
//**********************************************************************************
// METHOD : mat-mat prod (float)
// DESC   : Multiplies C = A * B and returns matrix prod C
//          result, for GFLOAT types
// ARGS   : GTMatrix m factor
// RETURNS: none
//**********************************************************************************
template<>
void matmat_prod<GFLOAT>(GTMatrix<GFLOAT> &C, const GTMatrix<GFLOAT> &A, const GTMatrix<GFLOAT> &B) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( A.size(2) != B.size(1) ) {
    cout << "GMTK::matmat_prod:incompatible matrix"<< endl;
    exit(1);
  }
  #endif


  #if !defined(_G_USE_GBLAS)
  for ( GSIZET i=0; i<C.size(1); i++ ) {
    for ( GSIZET j=0; j<C.size(2); j++ ) {
      C(i,j) = 0.0;
      for ( GSIZET k=0; k<A.size(2); k++ ) {
        C(i,j) += A(i,k) * B(k,j);
      }
    }
  }
  #else
  GSIZET a1=A.size(1), a2 = A.size(2);
  GSIZET b1=B.size(1), b2 = B.size(2);
  fmxm(C.data().data(),
       const_cast<GFLOAT*>(A.data().data()),
       &a1,&a2,
       const_cast<GFLOAT*>(B.data().data()),
       &b1, &b2, &szMatCache_);
  #endif
  

} // end of operator * mat-mat, GFLOAT


//**********************************************************************************
//**********************************************************************************
// METHOD : mat-mat prod (double)
// DESC   : Multiplies C = A * B and returns matrix prod C
//          result, for GDOUBLE types
// ARGS   : GTMatrix m factor
// RETURNS: none
//**********************************************************************************
template<>
void matmat_prod<GDOUBLE>(GTMatrix<GDOUBLE> &C, const GTMatrix<GDOUBLE> &A, const GTMatrix<GDOUBLE> &B) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( A.size(2) != B.size(1) ) {
    cout << "GMTK::matmat_prod:incompatible matrix"<< endl;
    exit(1);
  }
  #endif


  #if !defined(_G_USE_GBLAS)
  for ( GSIZET i=0; i<C.size(1); i++ ) {
    for ( GSIZET j=0; j<C.size(2); j++ ) {
      C(i,j) = 0.0;
      for ( GSIZET k=0; k<A.size(2); k++ ) {
        C(i,j) += A(i,k) * B(k,j);
      }
    }
  }
  #else
  GSIZET a1=A.size(1), a2 = A.size(2);
  GSIZET b1=B.size(1), b2 = B.size(2);
  dmxm(C.data().data(),
       const_cast<GDOUBLE*>(A.data().data()),
       &a1,&a2,
       const_cast<GDOUBLE*>(B.data().data()),
       &b1, &b2, &szMatCache_);
  #endif
  

} // end of operator * mat-mat, GDOUBLE


//**********************************************************************************
//**********************************************************************************
// METHOD : mat-mat prod (quad)
// DESC   : Multiplies C = A * B and returns matrix prod C
//          result, for GQUAD types
// ARGS   : GTMatrix m factor
// RETURNS: none
//**********************************************************************************
template<>
void matmat_prod<GQUAD>(GTMatrix<GQUAD> &C, const GTMatrix<GQUAD> &A, const GTMatrix<GQUAD> &B) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( A.size(2) != B.size(1) ) {
    cout << "GMTK::matmat_prod:incompatible matrix"<< endl;
    exit(1);
  }
  #endif


  #if !defined(_G_USE_GBLAS)
  for ( GSIZET i=0; i<C.size(1); i++ ) {
    for ( GSIZET j=0; j<C.size(2); j++ ) {
      C(i,j) = 0.0;
      for ( GSIZET k=0; k<A.size(2); k++ ) {
        C(i,j) += A(i,k) * B(k,j);
      }
    }
  }
  #else
  GSIZET a1=A.size(1), a2 = A.size(2);
  GSIZET b1=B.size(1), b2 = B.size(2);
  qmxm(C.data().data(),
       const_cast<GQUAD*>(A.data().data()),
       &a1,&a2,
       const_cast<GQUAD*>(B.data().data()),
       &b1, &b2, &szMatCache_);
  #endif
  

} // end of operator * mat-mat, GQUAD



//**********************************************************************************
//**********************************************************************************
// METHOD : constrain2sphere(1)
// DESC   : Project/constrain input 3-vector to sphere:
//          -  -     -                           -   -  -
//          |vx|     |(r^2-x^2)   -xy     -xz    |   |vx|
//        P |vy| =   |   -yx    (r^2-y^2) -yz    |   |vy|
//          |vz|     |   -zx      -zy   (r^2-z^2)|   |vz|
//          -  -     -                           -   -  -
//
//        This is derived from a Lagrange multiplier constraint
//        that requires all vectors, v, to be normal to radial
//        vector, s.t. x.v = 0. 
//          
// ARGS   : grid : Grid. If not of the correct type, nothing is done
//          v    : Array of vector components to be constrained
//          Pv   : Projected vector. Must be different from v
// RETURNS: none
//**********************************************************************************
template<>
void constrain2sphere(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &v, GTVector<GTVector<GFTYPE>*> &Pv)
{

  if ( grid.gtype() != GE_2DEMBEDDED || v.size() != 3 ) return;

  assert( v.size() >= 3 && "Incompatible dimensionality");

  GSIZET nxy = grid.ndof();
  GFTYPE r2, x, y, z;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  for ( GSIZET j=0; j<nxy; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    r2 = x*x + y*y + z*z;
    (*Pv[0])[j] =  (*v[0])[j]*(r2-x*x) - (*v[1])[j]*x*y      - (*v[2])[j]*x*z;
    (*Pv[1])[j] = -(*v[0])[j]*y*x      + (*v[1])[j]*(r2-y*y) - (*v[2])[j]*y*z;
    (*Pv[2])[j] = -(*v[0])[j]*z*x      - (*v[1])[j]*z*y      + (*v[2])[j]*(r2-z*z);
   }

} // end of method constrain2sphere (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : constrain2sphere(2)
// DESC   : Project/constrain input 3-vector to sphere:
//          -  -     -                           -   -  -
//          |vx|     |(r^2-x^2)   -xy     -xz    |   |vx|
//        P |vy| =   |   -yx    (r^2-x^2) -yz    |   |vy|
//          |vz|     |   -zx      -zy   (r^2-z^2)|   |vz|
//          -  -     -                           -   -  -
//
//        This is derived from a Lagrange multiplier constraint
//        that requires all vectors, v, to be normal to radial
//        vector, s.t. x.v = 0. 
//          
// ARGS   : grid : Grid. If not of the correct type, nothing is done
//          v    : Array of vector components to be constrained, 
//                 modified on output to be projected components
// RETURNS: none
//**********************************************************************************
template<>
void constrain2sphere(GGrid &grid, GTVector<GTVector<GFTYPE>*> &v)
{

  if ( grid.gtype() != GE_2DEMBEDDED || v.size() < 3 ) return;

  GSIZET nxy = grid.ndof();
  GFTYPE ri, r2, x, y, z;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  x = (*xnodes)[0][0]; y = (*xnodes)[1][0]; z = (*xnodes)[2][0];
  r2 = x*x + y*y + z*z;
  ri = 1.0/r2;
  for ( GSIZET j=0; j<nxy; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    (*v[0])[j] = ( (*v[0])[j]*(r2-x*x) - (*v[1])[j]*x*y      - (*v[2])[j]*x*z )*ri;
    (*v[1])[j] = (-(*v[0])[j]*y*x      + (*v[1])[j]*(r2-y*y) - (*v[2])[j]*y*z )*ri;
    (*v[2])[j] = (-(*v[0])[j]*z*x      - (*v[1])[j]*z*y      + (*v[2])[j]*(r2-z*z) )*ri;
   }

} // end of method constrain2sphere (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderiv
// DESC   : Compute tensor product derivative in specified direction
//          of specified field, u, in ref space, using grid object.
//          Compute
//            du = [ I_X_I_X_Dx, or
//                   I_X_Dy_X_I, or
//                   Dz_X_I_X_I].
//     
//          depending on whether idir = 1, 2, or 3, respectively,
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          grid   : GGrid object 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          idir   : coordinate direction (1, 2, or 3)
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          du     : vector of length of u containing the derivative.
//             
// RETURNS:  none
//**********************************************************************************
template<>
void compute_grefderiv(GGrid &grid, GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                       GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du)
{
  GSIZET               ibeg, iend; // beg, end indices for global array
  GBOOL                bembedded;
  GTVector<GSIZET>     N(GDIM);
  GTMatrix<GFTYPE>    *Di;         // element-based 1d derivative operators
  GElemList           *gelems = &grid.elems();

  bembedded = grid.gtype() == GE_2DEMBEDDED;

#if defined(_G_IS2D)
  switch (idir) {
  case 1:
    for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
      GMTK::I2_X_D1(*Di, u, N[0], N[1], du); 
    }
    break;
  case 2:
    for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
      GMTK::D2_X_I1(*Di, u, N[0], N[1], du); 
    }
    break;
  case 3:
    assert( GDIM == 3
         && "Only GDIM reference derivatives");
    du = 0.0; //u;
    break;
  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  u.range_reset(); // reset to global range
  du.range_reset();

#elif defined(_G_IS3D)

  switch (idir) {
  case 1:
    for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( GSIZET k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
      GMTK::I3_X_I2_X_D1(*Di, u, N[0], N[1], N[2], du); 
    }
    break;

  case 2:
    for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( GSIZET k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
      GMTK::I3_X_D2_X_I1(*Di, u, N[0], N[1], N[2], du); 
    }
    break;

  case 3:
    for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( GSIZET k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(2)->getDerivMatrix(!dotrans); 
      GMTK::D3_X_I2_X_I1(*Di, u, N[0], N[1], N[2], du); 
    }
    break;

  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  u.range_reset(); // reset global vec to globalrange
  du.range_reset();

#endif

} // end of method compute_grefderiv


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivs
// DESC   : Compute tensor product derivs of specified field, u, in ref space
//          for grid, using grid object to determine which to compute. Compute:
//            du = [ I_X_I_X_Dx
//                   I_X_Dy_X_I
//                   Dz_X_I_X_I].
//     
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          grid   : GGrid object 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          du     : vector of length 2 or 3 containing the derivatives.
//                   If using GE_REGULAR in 2D, we only need to vector
//                   elements; else we need 3. These should be allocated globally.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
template<>
void compute_grefderivs(GGrid &grid, GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                        GBOOL dotrans, GTVector<GTVector<GFTYPE>*> &du)
{
  assert(du.size() >= GDIM
       && "Insufficient number of derivatives specified");
  


  GBOOL                        bembedded;
  GINT                         nxy;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<GFTYPE>*>  Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &grid.elems();

  bembedded = grid.gtype() == GE_2DEMBEDDED;
  assert(( (bembedded && du.size()>=3) 
        || (!bembedded&& du.size()>=GDIM) )
       && "Insufficient number of derviatves provided");

  nxy = bembedded ? GDIM+1 : GDIM;

#if defined(_G_IS2D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( GSIZET k=0; k<nxy ; k++ ) du[k]->range(ibeg, iend);
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::I2_X_D1(*Di[0], u, N[0], N[1], *du[0]); 
    GMTK::D2_X_I1(*Di[1], u, N[0], N[1], *du[1]); 
#if 0
    if ( bembedded ) { // ref 3-deriv is just W u:
      *du[2] = u;  
    }
#endif
  }
  u.range_reset(); // reset to global range
  for ( GSIZET k=0; k<GDIM; k++ ) du[k]->range_reset();

#elif defined(_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( GSIZET k=0; k<GDIM; k++ ) du[k]->range(ibeg, iend);
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(2)->getDerivMatrix(!dotrans); 
    GMTK::I3_X_I2_X_D1(*Di[0], u, N[0], N[1], N[2], *du[0]); 
    GMTK::I3_X_D2_X_I1(*Di[1], u, N[0], N[1], N[2], *du[1]); 
    GMTK::D3_X_I2_X_I1(*Di[2], u, N[0], N[1], N[2], *du[2]); 
  }
  u.range_reset(); // reset global vec to globalrange
  for ( GSIZET k=0; k<nxy; k++ ) du[k]->range_reset();

#endif

} // end of method compute_grefderivs


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivsW
// DESC   : Compute tensor product derivs of specified field, u, in ref space 
//          for grid using grid object to determine which to compute, and 
//          include weights.
//          Compute:
//            du = (Mz_X_My_Mx) [ I_X_I_X_Dx
//                                I_X_Dy_X_I
//                                Dz_X_I_X_I].
//     
//          where Dx, Dy, Dz are 1d derivative objects from basis functions, and
//          Mi are the (diagonal) 1d-weights (or 'mass functions'). This can be 
//          re-written as
//            du = [ Mz_X_My_X_(Mx Dx)
//                   Mz_X_(My Dy)_X_Mx
//                  (Mz Dz)_X_My_X_Mx],
//           with comparable expressions for 2d.
// ARGS   : 
//          grid   : GGrid object 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          du     : vector of length 2 or 3 containing the derivatives.
//                   If using GE_REGULAR in 2D, we only need two vector
//                   elements; else we need 3. These should be allocated globally.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
template<>
void compute_grefderivsW(GGrid &grid, GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                         GBOOL dotrans, GTVector<GTVector<GFTYPE>*> &du)
{
  assert(du.size() >= GDIM
       && "Insufficient number of derivatives specified");
  

  GBOOL                         bembedded;
  GINT                          nxy;
  GSIZET                        k;
  GSIZET                        ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>              N(GDIM);
  GTVector<GTVector<GFTYPE>*>   W(GDIM);    // element weights
  GTVector<GTMatrix<GFTYPE>*>   Di(GDIM);   // element-based 1d derivative operators
  GElemList                    *gelems = &grid.elems();

  bembedded = grid.gtype() == GE_2DEMBEDDED;
  assert(( (bembedded && du.size()>=3)
        || (!bembedded&& du.size()>=GDIM) )
       && "Insufficient number of derviatves provided");

  nxy = bembedded ? GDIM+1 : GDIM;

#if defined(_G_IS2D)
  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( k=0; k<nxy ; k++ ) du[k]->range(ibeg, iend);
    for ( k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    for ( k=0; k<GDIM; k++ ) W[k]= (*gelems)[e]->gbasis(k)->getWeights();
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrixW (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrixW(!dotrans);
    GMTK::Dg2_X_D1(*Di[0], *W [1], u, etmp, *du[0]); 
    GMTK::D2_X_Dg1(*W [0], *Di[1], u, etmp, *du[1]); 
#if 0
    if ( bembedded ) { // ref 3-deriv is just W u:
      k = 0;
      for ( GSIZET j=0; j<N[1]; j++ ) {
        for ( GSIZET i=0; i<N[0]; i++ ) {
          (*du[2])[k] = u[k]*(*W[0])[i]*(*W[1])[j];  
          k++;
        }
      }
    }
#endif
  }
  u.range_reset(); // reset global vec to global range
  for ( k=0; k<nxy; k++ ) du[k]->range_reset();

#elif defined(_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( k=0; k<GDIM; k++ ) du[k]->range(ibeg, iend);
    for ( k=0; k<GDIM; k++ ) {
      N[k]= (*gelems)[e]->size(k);
      W[k]= (*gelems)[e]->gbasis(k)->getWeights();
    }
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrixW (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrixW(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(2)->getDerivMatrixW(!dotrans); 
    GMTK::Dg3_X_Dg2_X_D1(*Di[0], *W [1], *W [2], u, etmp, *du[0]); 
    GMTK::Dg3_X_D2_X_Dg1(*W [0], *Di[1], *W [2], u, etmp, *du[1]); 
    GMTK::D3_X_Dg2_X_Dg1(*W [0], *W [1], *Di[2], u, etmp, *du[2]); 
  }
  u.range_reset(); // reset global vec to globalrange
  for ( k=0; k<GDIM; k++ ) du[k]->range_reset();

#endif

} // end of method compute_grefderivsW


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefdiv
// DESC   : Compute tensor product 'divergence' of input fields in ref space
//          for grid:
//             Div u = [I_X_I_X_Dx     |u1|
//                      I_X_Dy_X_I   . |u2|
//                      Dz_X_I_X_I]    |u3|
//          or
//     
//             Div u = [I_X_I_X_DxT)   |u1|
//                      I_X_DyT_X_I .  |u2|
//                      DzT_X_I_X_I]   |u3|
//          where Dx(T), Dy(T), Dz(T) are 1d derivative objects from basis functions
// ARGS   : 
//          grid   : GGrid object 
//          u      : input vector field whose divergence we want, allocated globally 
//                   (e.g., for all elements). Must have GDIM components, unless
//                   we're using an embedded grid, when GDIM=2, when it should have
//                   3 components. Will not be altered on exit.
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          divu   : scalar result
//             
// RETURNS:  none
//**********************************************************************************
template<>
void compute_grefdiv(GGrid &grid, GTVector<GTVector<GFTYPE>*> &u, GTVector<GFTYPE> &etmp,
                     GBOOL dotrans, GTVector<GFTYPE> &divu)
{

  GBOOL                        bembedded;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GTVector<GFTYPE>*>  W(GDIM);    // element 1/weights
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<GFTYPE>*>  Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &grid.elems();

  bembedded = grid.gtype() == GE_2DEMBEDDED;
#if 0
  assert(( (bembedded && u.size()==GDIM) 
        || (!bembedded&& u.size()==GDIM) )
       && "Insufficient number of vector field components provided");
#else
  assert(  u.size()==GDIM  
       && "Insufficient number of vector field components provided");
#endif

  divu = 0.0;

#if defined(_G_IS2D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    // restrict global vecs to local range
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    etmp.resizem((*gelems)[e]->nnodes());
    GMTK::I2_X_D1(*Di[0], *u[0], N[0], N[1], etmp); // D1 u1
    divu += etmp;
    GMTK::D2_X_I1(*Di[1], *u[1], N[0], N[1], etmp); // D2 u2
    divu += etmp;
#if 0
    if ( bembedded ) divu += *u[2]; // D3 acts as I
#endif
  }
  divu.range_reset();  // reset range to global scope
  for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range_reset(); 

#elif defined(_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    etmp.resizem((*gelems)[e]->nnodes());
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 

    GMTK::I3_X_I2_X_D1(*Di[0], *u[0], N[0], N[1], N[2], etmp); // D1 u1
    divu += etmp;
    GMTK::I3_X_D2_X_I1(*Di[1], *u[1], N[0], N[1], N[2], etmp); // D2 u2
    divu += etmp;
    GMTK::D3_X_I2_X_I1(*Di[2], *u[2], N[0], N[1], N[2], etmp); // D3 u3
    divu += etmp;
  }
  divu.range_reset();  // reset range to global scope
  for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range_reset(); 

#endif


} // end, method compute_grefdiv


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefdiviW
// DESC   : Compute tensor product 'weighted divergence' of input fields in ref space
//          for grid:
//             Div u = [iM_X_iM_X_ iWDx    |u1|
//                      iM_X_iMDy_X_iM   . |u2|
//                      iMDz_X_iM_X_iM]    |u3|
//          or
//     
//             Div u = [iM_X_iM_X_iM DxT)    |u1|
//                      iM_X_iM DyT)_X_iM .  |u2|
//                      iM DzT)_X_iM_X_iM]   |u3|
//          where Dx(T), Dy(T), Dz(T) are 1d derivative objects from basis functions,
//          and iM are the 1d inverse weights     
// ARGS   : 
//          grid   : GGrid object 
//          u      : input vector field whose divergence we want, allocated globally 
//                   (e.g., for all elements). Must have GDIM components, unless
//                   we're using an embedded grid, when GDIM=2, when it should have
//                   3 components. Will be altered on exit.
//          etmp   : tmp array (possibly resized here) for element-based ops.
//                   Is not global.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//          divu   : scalar result
//             
// RETURNS:  none
//**********************************************************************************
template<>
void compute_grefdiviW(GGrid &grid, GTVector<GTVector<GFTYPE>*> &u, GTVector<GFTYPE> &etmp,
                       GBOOL dotrans, GTVector<GFTYPE> &divu)
{

  GBOOL                        bembedded;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GTVector<GFTYPE>*>  iW(GDIM);   // element 1/weights
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<GFTYPE>*>  Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &grid.elems();

  bembedded = grid.gtype() == GE_2DEMBEDDED;
  assert(( (bembedded && u.size()==3) 
        || (!bembedded&& u.size()==GDIM) )
       && "Insufficient number of vector field components provided");

  divu = 0.0;

#if defined(_G_IS2D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    // restrict global vecs to local range
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    for ( GSIZET k=0; k<GDIM; k++ ) iW[k]= (*gelems)[e]->gbasis(k)->getiWeights();
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrixiW (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrixiW(!dotrans);
    etmp.resizem((*gelems)[e]->nnodes());
    GMTK::Dg2_X_D1(*Di[0], *iW[1], *u[0], etmp, divu); // W X W^-1 D1 u1
    GMTK::D2_X_Dg1(*iW[0], *Di[1], *u[1], etmp, *u[0]); // W^^-1 D2 X W u2
    divu += *u[0];
#if 0
    if ( bembedded ) divu += *u[2]; // D3 acts as I
#endif
#if 0
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::I2_X_D1(*Di[0], u, N[0], N[1], *du[0]); 
    GMTK::D2_X_I1(*Di[1], u, N[0], N[1], *du[1]); 
#endif
  }
  divu.range_reset();  // reset range to global scope
  for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range_reset(); 

#elif defined(_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    for ( GSIZET k=0; k<GDIM; k++ ) iW[k]= (*gelems)[e]->gbasis(k)->getiWeights();
    etmp.resizem((*gelems)[e]->nnodes());
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrixiW (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrixiW(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(1)->getDerivMatrixiW(!dotrans); 

    GMTK::Dg3_X_Dg2_X_D1(*Di[0], *iW [1], *iW [2], *u[0], etmp, divu); 
    GMTK::Dg3_X_D2_X_Dg1(*iW [0], *Di[1], *iW [2], *u[1], etmp, *u[0]); 
    divu += *u[0];
    GMTK::D3_X_Dg2_X_Dg1(*iW [0], *iW [1], *Di[2], *u[2], etmp, *u[0]); 
    divu += *u[0];
  }
  divu.range_reset();  // reset range to global scope
  for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range_reset(); 

#endif


} // end, method compute_grefdiviW



//**********************************************************************************
//**********************************************************************************
// METHOD : vsphere2cart 
// DESC   : Convert vector field from spherical coords to Cartesian.
//
// ARGS   : grid : Grid. If not of the correct type, nothing is done
//          vsph : Array of vector components. If we have GE_2DEMBEDDED grid,
//                 there must be at least 2 components, and only the firs 2
//                 are used, and assumed to be latitudual, and longitudinal
//                 respectively. If grid is a 3D spherical grid, then
//                 vector components are assumed to be (r, lat, long).
//                 If grid is REGULAR, then this transformation cannot be done.
//          vtype: Vector type of spherical coords GVECTYPE_(PHYS, CONTRAVAR, COVAR)
//          vcart: Cartesian vector field. Must have at least 3 components, and
//                 are returned as (x, y, z).
// RETURNS: none
//**********************************************************************************
template<>
void vsphere2cart(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &vsph, GVectorType vtype, GTVector<GTVector<GFTYPE>*> &vcart)
{

  if      ( GDIM == 2 && grid.gtype() == GE_2DEMBEDDED ) {
    assert( vsph.size() >= 2 && "GE_2DEMBEDDED grid requires 2 spherical components");
  }
  else if ( grid.gtype() == GE_DEFORMED ) {
    assert( vsph.size() >= 3 && "GE_DEFORMED grid requires 3 spherical components");
  }
  else if ( grid.gtype() != GE_REGULAR ) {
    assert( FALSE && "GE_REGULAR grid will not allow this transformation");
  }


  GSIZET           nxy = grid.ndof();
  GFTYPE           x, y, z, tiny;
  GFTYPE           phi, r, theta;
  GFTYPE           gpp, gtt;
  GFTYPE           vthcontra, vphicontra;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  tiny = numeric_limits<GFTYPE>::epsilon();

  //   v_i_cart = vtheta dx_i/dtheta + vphi dx_i/dphi
  // where
  //   vtheta, vhi are _contravariant_ (upper-index) components
  // Note: Metric is orthogonal:
  //   g_ij = (1, h_theta^2, h_phi^2) = (1, r^2, (r cos(theta))^2 )
  if ( grid.gtype() == GE_2DEMBEDDED ) {
    if ( vtype == GVECTYPE_PHYS ) { // vsph are physical components
      for ( GSIZET j=0; j<nxy; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
        r     = x*x + y*y + z*z;
        theta = asin(z/r);
        phi   = atan2(y,x);
        (*vcart[0])[j] = -(*vsph[0])[j]*  sin(theta)*cos(phi) 
                       -  (*vsph[1])[j]*             sin(phi);
        (*vcart[1])[j] = -(*vsph[0])[j]*  sin(theta)*sin(phi) 
                       +  (*vsph[1])[j]*             cos(phi);
        (*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta);
      }
    }
    else if ( vtype == GVECTYPE_COVAR ) { // vsph are covar. components
      for ( GSIZET j=0; j<nxy; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
        r     = x*x + y*y + z*z;
        theta = asin(z/r);
        phi   = atan2(y,x);
        vthcontra  = (*vsph[0])[j];
        vphicontra = (*vsph[1])[j];
        gtt        = r*r; gpp = pow(r*cos(theta),2);
        vthcontra  = (*vsph[0])[j]/(gtt+tiny);
        vphicontra = (*vsph[1])[j]/(gpp+tiny);
        (*vcart[0])[j] = -vthcontra *r*sin(theta)*cos(phi) 
                       -  vphicontra*r*cos(theta)*sin(phi);
        (*vcart[1])[j] = -vthcontra *r*sin(theta)*sin(phi) 
                       +  vphicontra*r*cos(theta)*cos(phi);
        (*vcart[2])[j] =  vthcontra *r*cos(theta);
      }
    }
    else if ( vtype == GVECTYPE_CONTRAVAR ) { // vsph are contravar. components
      for ( GSIZET j=0; j<nxy; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
        r     = x*x + y*y + z*z;
        theta = asin(z/r);
        phi   = atan2(y,x);
        vthcontra  = (*vsph[0])[j];
        vphicontra = (*vsph[1])[j];
        (*vcart[0])[j] = -vthcontra *r*sin(theta)*cos(phi) 
                       -  vphicontra*r*cos(theta)*sin(phi);
        (*vcart[1])[j] = -vthcontra *r*sin(theta)*sin(phi) 
                       +  vphicontra*r*cos(theta)*cos(phi);
        (*vcart[2])[j] =  vthcontra *r*cos(theta);
      }
    }
  }

  if ( grid.gtype() == GE_DEFORMED ) {
    if      ( vtype == GVECTYPE_PHYS ) {
      for ( GSIZET j=0; j<nxy; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
        r     = x*x + y*y + z*z;
        theta = asin(z/r);
        phi   = atan2(y,x);
        gtt        = r; gpp = r*cos(theta);
        vthcontra  = (*vsph[1])[j]/(gtt+tiny);
        vphicontra = (*vsph[2])[j]/(gpp+tiny);
        (*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
                       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
                       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
        (*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
                       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
                       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
        (*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
                       +  (*vsph[1])[j]*r*cos(theta);
      }
    }
    else if ( vtype == GVECTYPE_COVAR ) {
      for ( GSIZET j=0; j<nxy; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
        r     = x*x + y*y + z*z;
        theta = asin(z/r);
        phi   = atan2(y,x);
        gtt        = r*r; gpp = pow(r*cos(theta),2);
        vthcontra  = (*vsph[1])[j]/(gtt+tiny);
        vphicontra = (*vsph[2])[j]/(gpp+tiny);
        (*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
                       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
                       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
        (*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
                       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
                       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
        (*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
                       +  (*vsph[1])[j]*r*cos(theta);
      }
    }
    else if ( vtype == GVECTYPE_CONTRAVAR ) {
      for ( GSIZET j=0; j<nxy; j++ ) {
        x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
        r     = x*x + y*y + z*z;
        theta = asin(z/r);
        phi   = atan2(y,x);
        vthcontra  = (*vsph[0])[j];
        vphicontra = (*vsph[1])[j];
        (*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
                       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
                       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
        (*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
                       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
                       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
        (*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
                       +  (*vsph[1])[j]*r*cos(theta);
      }
    }
  }

} // end of method vsphere2cart


//**********************************************************************************
//**********************************************************************************
// METHOD : energy
// DESC   : 
//             Compute volume-integrated mean of energy from input 
//             vector field:
//               energy = 0.5 * Int _u_^2 dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 2, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : do max norm, instead of L2
// RETURNS: GFTYPE energy
//**********************************************************************************
template<>
GFTYPE energy(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
  GDOUBLE                     ener, local;
  GC_COMM                     comm = grid.get_comm();

  // Find _u_^2 = Sum_l u_l ^2
 *tmp[1] = *u[0]; tmp[1]->rpow(2);
  for ( GINT l=1; l<u.size(); l++ ) {
    *tmp[0] = *u[l]; tmp[0]->rpow(2);
    *tmp[1] += *tmp[0];
  }

  if ( ismax ) {
    ener =  0.5*tmp[1]->amax();
    if ( isglobal ) {
      local = ener;
      GComm::Allreduce(&local, &ener, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
    }
  }
  else {
    ener  = static_cast<GDOUBLE>(grid.integrate(*tmp[1], *tmp[0], isglobal));
    ener *=  0.5*static_cast<GDOUBLE>(grid.ivolume());
  }

  return static_cast<GFTYPE>(ener);

} // end of method energy


//**********************************************************************************
//**********************************************************************************
// METHOD : enstrophy
// DESC   : 
//             Compute volume-integrated mean 
//                 1/2 Int |curl u |^2 dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 4, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE enstrophy
//**********************************************************************************
template<>
GFTYPE enstrophy(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
  assert(tmp.size() >= 4 && "Insufficient temp space");

  
  GINT                        ibeg, iend;
  GDOUBLE                     enst, local;
  GC_COMM                     comm = grid.get_comm();
  GTVector<GFTYPE>           *cc;
  GTVector<GTVector<GFTYPE>*> utmp(3);

  utmp[0] = tmp[0];
  utmp[1] = tmp[1];
  cc      = tmp[2];

 *tmp[3] = 0.0;
  if ( u.size() == 3 ) {
    for ( GINT l=0; l<u.size(); l++ ) {
      GMTK::curl(grid, u, l+1, utmp, *cc);
      cc->rpow(2);
     *tmp[3] += *cc;
    }
  }
  else if ( u.size() == 2 ) {
    GMTK::curl(grid, u, 3, utmp, *cc);
    cc->rpow(2);
   *tmp[3] += *cc;
  }

  if ( ismax ) {
    enst =  0.5*static_cast<GDOUBLE>(tmp[3]->amax());
    if ( isglobal ) {
      local = enst;
      GComm::Allreduce(&local, &enst, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
    }
  }
  else {
    enst  = static_cast<GDOUBLE>(grid.integrate(*tmp[3], *tmp[0], isglobal));
    enst *=  0.5*grid.ivolume();
  }

  return static_cast<GFTYPE>(enst);

} // end of method enstrophy



//**********************************************************************************
//**********************************************************************************
// METHOD : helicity
// DESC   : 
//             Compute volume-integrated mean 
//                 Int |curl u \dot u| dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 4, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute abs max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE helicity
//**********************************************************************************
template<>
GFTYPE helicity(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
  assert(tmp.size() >= 4 && "Insufficient temp space");

  
  GDOUBLE                     hel, local;
  GC_COMM                     comm = grid.get_comm();
  GTVector<GFTYPE>           *cc;
  GTVector<GTVector<GFTYPE>*> utmp(3);

  utmp[0] = tmp[0];
  utmp[1] = tmp[1];
  cc      = tmp[2];

 *tmp[3] = 0.0;
  if ( u.size() == 3 ) {
    for ( GINT l=0; l<3; l++ ) {
      GMTK::curl(grid, u, l+1, utmp, *cc);
      if ( u.size() > l ) cc->pointProd(*u[l]);
     *tmp[3] += *cc;
    }
  }

  if ( ismax ) {
    hel =  static_cast<GDOUBLE>(tmp[3]->amax());
    if ( isglobal ) {
      local = hel;
      GComm::Allreduce(&local, &hel, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
    }
  }
  else {
    hel  = static_cast<GDOUBLE>(grid.integrate(*tmp[3], *tmp[0], isglobal));
    hel *=  static_cast<GDOUBLE>(grid.ivolume());
  }

  return static_cast<GFTYPE>(hel);

} // end of method helicity 



//**********************************************************************************
//**********************************************************************************
// METHOD : relhelicity
// DESC   : 
//             Compute volume-integrated mean 
//                 Int |curl u \dot u| /(|u| {curl u|) dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 5, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute abs max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE helicity
//**********************************************************************************
template<>
GFTYPE relhelicity(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
  assert(tmp.size() >= 5 && "Insufficient temp space");

  
  GDOUBLE                     local, rhel;
  GC_COMM                     comm = grid.get_comm();
  GTVector<GFTYPE>           *cc;
  GTVector<GTVector<GFTYPE>*> utmp(3);

  utmp[0] = tmp[0];
  utmp[1] = tmp[1];
  cc      = tmp[2];

 *tmp[3] = 0.0;
 *tmp[4] = 0.0;

  // Compute u. curl u:
  if ( u.size() == 3 ) {
    for ( GINT l=0; l<3; l++ ) {
      GMTK::curl(grid, u, l+1, utmp, *cc);
      cc->pointProd(*u[l]);
     *tmp[3] += *cc;
    }
  }
  
  // Compute |curl u|:
  if ( u.size() == 3 ) {
    for ( GINT l=0; l<3; l++ ) {
      GMTK::curl(grid, u, l+1, utmp, *cc);
      cc->rpow(2);
     *tmp[4] += *cc;
    }
  }
  else if ( u.size() == 2 ) {
    GMTK::curl(grid, u, 3, utmp, *cc);
    cc->rpow(2);
   *tmp[4] += *cc;
  } 
  tmp[4]->rpow(0.5);


  // Compute |u|:
  *tmp[0] = 0.0;
  for ( GINT l=0; l<u.size(); l++ ) {
   *tmp[1] = *u[l];
    tmp[1]->rpow(2);
    *tmp[0] += *tmp[1];
  }
  tmp[0]->rpow(0.5);


  // Compute u. (curl u) /|u| |curl u| integrand:
  GFTYPE tiny = 100.0*numeric_limits<GFTYPE>::epsilon();

  tmp[0]->pointProd(*tmp[4]); // compute |u| |curl u|
  for ( GSIZET k=0; k<utmp[0]->size(); k++ ) {
    // (*tmp[1])[k] = fabs((*tmp[0])[k]) <= tiny ? 0.0 : (*tmp[3])[k]/(*tmp[0])[k];  
    if ( fabs((*tmp[0])[k]) <= tiny ) {
      (*tmp[1])[k] = 0.0;
    }
    else {
      (*tmp[1])[k] = (*tmp[3])[k]/(*tmp[0])[k];
    }
  }


  if ( ismax ) {
    rhel =  static_cast<GDOUBLE>(tmp[1]->amax());
    if ( isglobal ) {
      local = rhel;
      GComm::Allreduce(&local, &rhel, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
    }
  }
  else {
    rhel  = static_cast<GDOUBLE>(grid.integrate(*tmp[1], *tmp[0], isglobal));
    rhel *= static_cast<GDOUBLE>(grid.ivolume());
  }

  return static_cast<GFTYPE>(rhel);

} // end of method relhelicity 



//**********************************************************************************
//**********************************************************************************
// METHOD : energyinj
// DESC   : 
//             Compute volume-integrated mean energy injection
//                 Int |u \dot f| dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : velocity field; entire field used to compute energy
//          uf      : forcing field; each must be non-NULL on entry; if
//                 any are NULL, or uf.size == 0, then return value is 0.
//          tmp     : tmp vector of length at least 2, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute abs max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE energy injection rate
//**********************************************************************************
template<>
GFTYPE energyinj(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &u,  const GTVector<GTVector<GFTYPE>*> &uf, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{

  if ( uf.size() == 0 ) return 0.0;

  GBOOL bnull = FALSE;
  for ( GINT l=0; l<u.size(); l++ ) bnull = bnull || uf[l] == NULLPTR;

  if ( bnull ) return 0.0;

  assert(tmp.size() >= 2 && "Insufficient temp space");

  
  GDOUBLE                     einj, local;
  GC_COMM                     comm = grid.get_comm();


  u[0]->pointProd(*uf[0], *tmp[0]);
  for ( GINT l=1; l<u.size(); l++ ) {
	  assert(uf[l]!= NULLPTR && "NULL force not allowed");
    u[l]->pointProd(*uf[l], *tmp[1]);
   *tmp[0] += *tmp[1];
  }

  if ( ismax ) {
    einj =  static_cast<GDOUBLE>(tmp[0]->amax());

    if ( isglobal ) {
      local = einj;
      GComm::Allreduce(&local, &einj, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
    }

  }
  else {
    einj  = static_cast<GDOUBLE>(grid.integrate(*tmp[0], *tmp[1], isglobal));
    einj *=  static_cast<GDOUBLE>(grid.ivolume());
  }

  return static_cast<GFTYPE>(einj);

} // end of method energyinj



//**********************************************************************************
//**********************************************************************************
// METHOD : domathop
// DESC   : 
//             Carry out math (mainly differential) operation, based on 
//             specified string description.
//
//             Return in appropriate component of output array
//          
// ARGS   : 
//          grid    : grid object
//          uin     : input state. May be vector or scalar, in an
//                    order and number appropriate for specified operation
//          sop     : string operation
//          tmp     : tmp vector of length at least 1, each
//                    of same length as uin
//          uout    : output array
//          iuout   : which indices of uout contain valid output data; 
//                    sized to reflect valid number of components
// RETURNS: none
//**********************************************************************************
template<>
void domathop(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &uin,  const GString sop, GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GTVector<GFTYPE>*> &uout, GTVector<GINT> &iuout)
{

  GINT                        nxy;
  GTVector<GTVector<GFTYPE>*> tmp(3);

  if      ( "div"  == sop ) { // operates on a vector field...
    // produces a scalar field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 1   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    GMTK::grad(grid, *uin[0], 1, utmp, *uout[0]);
    iuout.resize(1); iuout[0] = 0;
    for ( auto j=1; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[j], j+1, utmp, *utmp[utmp.size()-1]);
      *uout[0] += *utmp[utmp.size()-1];
    }
  }
  else if ( "grad" == sop ) {  // operates on a scalar field...
    // produces a vector field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 1   && "Insufficient temp space");
    assert(uin  .size() >= 1   && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    iuout.resize(nxy); 
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[0], j+1, utmp, *uout[j]);
      iuout[j] = j; 
    }
  }
  else if ( "gradmag" == sop ) {  // operates on a scalar field...
    // produces a scalar field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 2   && "Insufficient temp space");
    assert(uin  .size() >= 1   && "Insufficient no. input components");
    assert(uout .size() >= 1   && "Insufficient no. output components");
    iuout.resize(1); iuout[0] = 0; 
    tmp[0] = utmp[0]; tmp[1] = utmp[1];
    GMTK::grad(grid, *uin[0], 1, tmp, *uout[0]);
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[0], j+1, tmp, *tmp[1]);
      tmp[1]->rpow(2);
      *uout[0] += *tmp[1];
    }
    uout[0]->rpow(0.5);
  }
  else if ( "curl" == sop ) { // operates on a vector field...
    // produces a vector field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 2   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    iuout.resize(nxy); 
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::curl(grid, uin, j+1, utmp, *uout[j]);
      iuout[j] = j; 
    }
  }
  else if ( "curlmag" == sop ) { // operates on a vector field...
    // produces a vector field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 3   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    iuout.resize(1); iuout[0] = 0; 
    tmp[0] = utmp[0]; tmp[1] = utmp[1]; tmp[2] = utmp[2];
    GMTK::curl(grid, uin, 1, utmp, *uout[0]);
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::curl(grid, uin, j+1, tmp, *tmp[2]);
      tmp[2]->rpow(2);
      *uout[0] += *tmp[2];
    }
    uout[0]->rpow(0.5);
  }
  else if ( "lapderivs" == sop ) { // operates on a scalar field
    // produces a vector field of each termin Laplacnin:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 3   && "Insufficient temp space");
    assert(uin  .size() >= 1   && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    for ( auto j=0; j<3; j++ ) tmp[j] = utmp[j];
    iuout.resize(nxy); 
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[0], j+1, tmp , *tmp[2]);
      GMTK::grad(grid, *tmp[2], j+1, utmp, *uout[j]);
      iuout[j] = j; 
    }
  }
  else {
    assert(FALSE && "Invalid math operation");
  }

} // end of method domathop 


} // end, namespace GMTK
