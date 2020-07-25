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

#include <assert.h>
#include <limits>
#include "gtypes.h"
#include "cff_blas.h"
//#include "ggrid.hpp"
//#include "ggrid_box.hpp"
//#include "ggrid_icos.hpp"

template<typename T> class GTVector;
template<typename T> class GTMatrix;
                     class GGrid;

extern GINT szMatCache_;
extern GINT szVecCache_;


namespace GMTK
{

  template<typename T>     
  void Pml(T l, T m, GTVector<T> &y, GTVector<T> &z);

  template<typename T>     
  void Rx3(T alpha, GTVector<T> &y, GTVector<T> &z);

  template<typename T>     
  void Ry3(T alpha, GTVector<T> &x, GTVector<T> &z);

  template<typename T>     
  void Rz3(T alpha, GTVector<T> &x, GTVector<T> &y);

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
  void matvec_prod(GTVector<T> &vret, const GTMatrix<T> &A, const GTVector<T> &b);

  template<typename T>     
  void matmat_prod(GTMatrix<T> &C, const GTMatrix<T> &A, const GTMatrix<T> &B);

  template<typename T>
  void cross_prod_k(GTVector<GTVector<T>*> &A, GINT *iind, GINT nind, 
                    GINT isgn, GTVector<GTVector<T>*> &C);

  template<typename T>
  void cross_prod_k(GTVector<T> &Ax, GTVector<T> &Ay,
                    GINT *iind, GINT nind, GINT isgn,
                    GTVector<T> &Cx, GTVector<T> &Cy);

  template<typename T>
  void cross_prod(GTVector<GTVector<T>*> &A, GTVector<GTVector<T>*> &B,
                  GINT *iind, GINT nind, GTVector<GTVector<T>*> &C);

  template<typename T>
  void cross_prod(GTVector<T> &Ax, GTVector<T> &Ay, GTVector<T> &Az,
                  GTVector<T> &Bx, GTVector<T> &By, GTVector<T> &Bz,
                  GINT *iind, GINT nind,
                  GTVector<T> &Cx, GTVector<T> &Cy, GTVector<T> &Cz);

  template<typename T>  
  void    normalize_euclidean(GTVector<GTVector<T>*> &x, GINT *iind, GINT nind, T x0=1);

  template<typename T>  
  void    normalizeL2(GGrid &grid, GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*    > &tmp, T u0);

  template<typename T>
  void saxpby(GTVector<T> &x, T a, GTVector<T> &y, T b); 

  template<typename T>  
  void     add(GTVector<T> &vret, const GTVector<T> &va, const GTVector<T> &vb, T a, T b);

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
  template<typename T>  
  void    grad(GGrid &grid, GTVector<T> &u, const GINT idir, 
               GTVector<GTVector<T>*> &tmp, GTVector<T> &grad);
  template<typename T>  
  void    curl(GGrid &grid, const GTVector<GTVector<T>*> &u, const GINT idir, 
               GTVector<GTVector<T>*> &tmp, GTVector<T> &curl);

  template<typename T>  
  void    constrain2sphere(GGrid &grid, const GTVector<GTVector<T>*> &v, GTVector<GTVector<T>*    > &Pv);
  template<typename T>  
  void    constrain2sphere(GGrid &grid, GTVector<GTVector<T>*> &v);
  template<typename T>  
  void    vsphere2cart(GGrid &grid, const GTVector<GTVector<T>*> &vsph, GVectorType vtype, GTVector<GTVector<T>*> &vcart);
  template<typename T>  
  void    vcart2sphere(GGrid &grid, const GTVector<GTVector<T>*> &vcart, GVectorType vtype, GTVector<GTVector<T>*> &vsph);
  template<typename T>  
  void    cart2latlon(const GTVector<GTVector<T>*> &vcart, GTVector<GTVector<T>*> &latlon);
  template<typename T>  
  void    cart2spherical(const GTVector<GTVector<T>*> &vcart, GTVector<GTVector<T>*> &rlatlon);
  template<typename T>  
  void    zero(GTVector<T> &v);
  template<typename T>  
  void    zero(GTVector<GTVector<T>*> &v);
  template<typename T>  
  void    zero(GTVector<GTVector<T>> &v);
  template<typename T>  
  GDOUBLE energy(GGrid &grid, const GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax=FALSE);
  template<typename T>  
  GDOUBLE enstrophy(GGrid &grid, const GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax=FALSE);
  template<typename T>  
  GDOUBLE helicity(GGrid &grid, const GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax=FALSE);
  template<typename T>  
  GDOUBLE relhelicity(GGrid &grid, const GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax=FALSE);
  template<typename T>  
  GDOUBLE energyinj(GGrid &grid, const GTVector<GTVector<T>*> &u, const GTVector<GTVector<T>*> &uf, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax=FALSE);
  template<typename T>  
  void domathop(GGrid &grid, const GTVector<GTVector<T>*> &uin, const GString sop, GTVector<GTVector<T>*> &utmp, GTVector<GTVector<T>*> &uout, GTVector<GINT> &iuout);

};

#include "gmtk.ipp"

#endif
