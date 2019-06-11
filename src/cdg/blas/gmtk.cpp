//==================================================================================
// Module       : gmtk.cpp
// Date         : 1/31/18 (DLR)
// Description  : Math TooKit: namespace of C routines for various
//                math functions
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gtypes.h"
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
    std::cout << "GMTK::D2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::I2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::I2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::I2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::Dg_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::Dg_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::Dg_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_I1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_I1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_I1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_Dg1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_Dg1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_Dg1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_D1" << "incompatible size" << std::endl;
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
    std::cout << "GMTK::D2_X_D1" << "incompatible size" << std::endl;
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
void add<GFLOAT>(GTVector<GFLOAT> &vret, GTVector<GFLOAT> &va, GTVector<GFLOAT> &vb, GFLOAT a, GFLOAT b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( va.size() < vret.size() || vb.size() < vret.size() ) {
    std::cout << "GTVector<T>::add: " << "incompatible size" << std::endl;
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
  fzaxpby(vret.data(), va.data(), &a, vb.data(), &b, &nn, &szVecCache_);
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
void add<GDOUBLE>(GTVector<GDOUBLE> &vret, GTVector<GDOUBLE> &va, GTVector<GDOUBLE> &vb, GDOUBLE a, GDOUBLE b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( va.size() < vret.size() || vb.size() < vret.size() ) {
    std::cout << "GTVector<T>::add: " << "incompatible size" << std::endl;
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
  dzaxpby(vret.data(), va.data(), &a, vb.data(), &b, &nn, &szVecCache_);
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
void add<GQUAD>(GTVector<GQUAD> &vret, GTVector<GQUAD> &va, GTVector<GQUAD> &vb, GQUAD a, GQUAD b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( va.size() < vret.size() || vb.size() < vret.size() ) {
    std::cout << "GTVector<T>::add: " << "incompatible size" << std::endl;
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
  qzaxpby(vret.data(), va.data(), &a, vb.data(), &b, &nn, &szVecCache_);
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
void matvec_prod<GFLOAT>(GTVector<GFLOAT> &vret, GTMatrix<GFLOAT> &A, GTVector<GFLOAT> &b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( b.size() < A.size(2) ) {
    std::cout << "GMTK::matvec_prod: " << "incompatible size" << std::endl;
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
  fmxv(vret.data(), A.data().data(), b.data(), &n1, &n2, &szMatCache_);
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
void matvec_prod<GDOUBLE>(GTVector<GDOUBLE> &vret, GTMatrix<GDOUBLE> &A, GTVector<GDOUBLE> &b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( b.size() < A.size(2) ) {
    std::cout << "GMTK::matvec_prod: " << "incompatible size" << std::endl;
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
  dmxv(vret.data(), A.data().data(), b.data(), &n1, &n2, &szMatCache_);
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
void matvec_prod<GQUAD>(GTVector<GQUAD> &vret, GTMatrix<GQUAD> &A, GTVector<GQUAD> &b) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( b.size() < A.size(2) ) {
    std::cout << "GMTK::matvec_prod: " << "incompatible size" << std::endl;
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
  qmxv(vret.data(), A.data().data(), b.data(), &n1, &n2, &szMatCache_);
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
void matmat_prod<GFLOAT>(GTMatrix<GFLOAT> &C, GTMatrix<GFLOAT> &A, GTMatrix<GFLOAT> &B) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( A.size(2) != B.size(1) ) {
    std::cout << "GMTK::matmat_prod:incompatible matrix"<< std::endl;
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
  fmxm(C.data().data(),A.data().data(),&a1,&a2,
       B.data().data(),&b1, &b2, &szMatCache_);
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
void matmat_prod<GDOUBLE>(GTMatrix<GDOUBLE> &C, GTMatrix<GDOUBLE> &A, GTMatrix<GDOUBLE> &B) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( A.size(2) != B.size(1) ) {
    std::cout << "GMTK::matmat_prod:incompatible matrix"<< std::endl;
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
  dmxm(C.data().data(),A.data().data(),&a1,&a2,
       B.data().data(),&b1, &b2, &szMatCache_);
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
void matmat_prod<GQUAD>(GTMatrix<GQUAD> &C, GTMatrix<GQUAD> &A, GTMatrix<GQUAD> &B) 
{
  #if defined(_G_BOUNDS_CHK)
  if ( A.size(2) != B.size(1) ) {
    std::cout << "GMTK::matmat_prod:incompatible matrix"<< std::endl;
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
  qmxm(C.data().data(),A.data().data(),&a1,&a2,
       B.data().data(),&b1, &b2, &szMatCache_);
  #endif
  

} // end of operator * mat-mat, GQUAD


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
//          curl : result
// RETURNS: none.
//**********************************************************************************
template<>
void curl(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &u, const GINT idir, 
          GTVector<GTVector<GFTYPE>*> &tmp, GTVector<GFTYPE> &curl)
{

  if ( GDIM == 2 && u.size() > GDIM ) {
    switch (idir) {
      case 1:
        grid.deriv(*u[2], 2, *tmp[0], curl);
        curl *= -1.0;
        break;
      case 2:
        grid.deriv(*u[2], 1, *tmp[0], curl);
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curl);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curl -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
    return;
  }

  if ( GDIM == 2 ) {
    switch (idir) {
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curl);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curl -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
    return;
  }

  if ( GDIM == 3 ) {
    switch (idir) {
      case 1:
        grid.deriv(*u[1], 3, *tmp[0], curl);
        grid.deriv(*u[2], 2, *tmp[0], *tmp[1]);
        curl -= *tmp[1];
        break;
      case 2:
        grid.deriv(*u[2], 1, *tmp[0], curl);
        grid.deriv(*u[0], 3, *tmp[0], *tmp[1]);
        curl -= *tmp[1];
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curl);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curl -= *tmp[1];
        break;
    }
  }
  return;

} // end of method curl


//**********************************************************************************
//**********************************************************************************
// METHOD : constrain2sphere(1)
// DESC   : Project input 3-vector to sphere:
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
//          v    : Array of vector components
//          Pv   : Projected vector. May be the same as V, in which case
//                 V is overwritten.
// RETURNS: none
//**********************************************************************************
template<typename T>
void constrain2sphere(GGrid &grid, GTVector<GTVector<GFTYPE>*> &v, GTVector<GTVector<GFTYPE>*> &Pv)
{

  if ( grid.gtype() != GE_2DEMBEDDED ) return;

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
  GTMatrix<GFTYPE>     *Di;         // element-based 1d derivative operators
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
#if 1
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
//                 vector components are assumed to be (r, theta, phi).
//                 If grid is REGULAR, then this transformation cannot be done.
//          vtype: Vector type of spherical coords (PHYS, CONTRAVAR, COVAR)
//          vcart: Cartesian vector field. Must have at least 3 components, and
//                 are returned as (x, y, z).
// RETURNS: none
//**********************************************************************************
template<>
void vsphere2cart(GGrid &grid, GTVector<GTVector<T>*> &vsph, GVectorType vtype, GTVector<GTVector<T>*> &vcart);
{

  if ( grid.gtype() != GE_2DEMBEDDED ) {
    assert( vsph.size() >= 2 && "GE_2DEMBEDDED grid requires 2 spherical components");
  }
  else if ( grid.gtype() != GE_DEFORMED ) {
    assert( vsph.size() >= 3 && "GE_DEFORMED grid requires 3 spherical components");
  else if ( grid.gtype() != GE_REGULAR ) {
    assert( FALSE && "GE_REGULAR grid will not allow this transformation");
  }


  GSIZET           nxy = grid.ndof();
  GFTYPE           x, y, z, tiny;
  GFTYPE           phi, r, theta;
  GFTYPE           hpp, hrr, htt;
  GFTYPE           vthcontra, vphicontra;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  tiny = std::numeric_limits<GFTY?PE>::epsilon();

  if ( grid.gtype() == GE_2DEMBEDDED ) {
    for ( GSIZET j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
      r     = x*x + y*y + z*z;
      theta = asin(z/rad);
      phi   = atan2(y,x);
      vthcontra  = (*vsph[0])[j];
      vphicontra = (*vsph[1])[j];
      if ( vtype == GVECTYPE_PHYS ) {
        htt        = r; hpp = r*cos(theta);
        vthcontra  = (*vsph[0])[j]/(htt+tiny);
        vphicontra = (*vsph[1])[j]/(hpp+tiny);
      }
      else if ( vtype == GVECTYPE_COVAR ) {
        htt        = r*r; hpp = pow(r*cos(theta),2);
        vthcontra  = (*vsph[0])[j]*(htt+tiny);
        vphicontra = (*vsph[1])[j]*(hpp+tiny);
      }

      (*vcart[0])[j] = -(*v[0])[j]*r*sin(theta)*cos(phi) 
                     -  (*v[1])[j]*r*cos(theta)*sin(phi);
      (*vcart[1])[j] = -(*v[0])[j]*r*sin(theta)*sin(phi) 
                     -  (*v[1])[j]*r*cos(theta)*cos(phi);
      (*vcart[2])[j] =  (*v[0])[j]*r*cos(theta);
     }
     return;
   }

  if ( grid.gtype() == GE_DEFORMED ) {
    for ( GSIZET j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
      r     = x*x + y*y + z*z;
      theta = asin(z/rad);
      phi   = atan2(y,x);
      vthcontra  = (*vsph[0])[j];
      vphicontra = (*vsph[1])[j];
      if ( vtype == GVECTYPE_PHYS ) {
        htt        = r; hpp = r*cos(theta);
        vthcontra  = (*vsph[1])[j]/(htt+tiny);
        vphicontra = (*vsph[2])[j]/(hpp+tiny);
      }
      else if ( vtype == GVECTYPE_COVAR ) {
        htt        = r*r; hpp = pow(r*cos(theta),2);
        vthcontra  = (*vsph[1])[j]/(htt+tiny);
        vphicontra = (*vsph[2])[j]/(hpp+tiny);
      }

      (*vcart[0])[j] =  (*v[0])[j]*  cos(theta)*cos(phi)
                     -  (*v[1])[j]*r*sin(theta)*cos(phi) 
                     -  (*v[2])[j]*r*cos(theta)*sin(phi);
      (*vcart[1])[j] =  (*v[0])[j]*  cos(theta)*sin(phi)
                     -  (*v[1])[j]*r*sin(theta)*sin(phi) 
                     -  (*v[2])[j]*r*cos(theta)*cos(phi);
      (*vcart[2])[j] =  (*v[0])[j]*  sin(theta)
                     +  (*v[1])[j]*r*cos(theta);
     }
     return;
   }

} // end of method vsphere2cart


} // end, namespace GMTK
