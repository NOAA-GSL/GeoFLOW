//==================================================================================
// Module       : gllbasis.hpp
// Date         : 1/19/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a spectral element nodal basis
//                Note: This basis is the Gauss-Lobatto-Legendre
//                basis, defined as 
//                  h(xi)_j = -1/(N(N+1)*L_N(xi_j))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_j)
//                where N is the order of the expansion, and L_N is the Legendre
//                polynomial of order N, and the xi_j are the nodal points. L_N
//                is obtained from the generalized Jacobi polynomial computed in
//                computeJacobi.
//
//                Template argument T refers to the typename in which the basis
//                is computed (e.g., FLOAT, DOUBLE, QUAD), and TE refers to the
//                argument type in which the computed quantities are _evaluated_,
//                usually at the same or lower size.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GNBasis
//==================================================================================

#if !defined(_GLLBASIS_HPP)
#define _GLLBASIS_HPP
#include "gtypes.h"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"


template<typename T, typename TE> 
class GLLBasis: public GNBasis<T,TE>
{
        static_assert(std::is_floating_point<T>::value && std::is_floating_point<TE>::value,"Illegal template type");


public:

                          GLLBasis();
                          GLLBasis(GINT  inorder);
                          GLLBasis(GINT  inorder, GINT  maxorder);
                          GLLBasis(const GLLBasis &);
virtual                  ~GLLBasis();

virtual  GBOOL            resize(GINT  order);

virtual  void             operator=(const GLLBasis &);
                   
virtual  T                getXimin ();
virtual  T                getXimax (); 
virtual  GINT             getOrder();
virtual  void             setOrder(GINT );

virtual GTVector<T>      *getXiNodesComp();
virtual GTVector<TE>     *getXiNodes();
virtual TE               *getXiNodes(TE *ret, GINT  num);
virtual void              getXiNodes(GTVector<TE> &ret);

virtual  GTVector<T>     *getWeightsComp();
virtual  GTVector<TE>    *getWeights();
virtual  TE              *getWeights(TE *ret, GINT  num);
virtual  void             getWeights(GTVector<TE> &ret);

virtual  GTVector<TE>    *getiWeights();
virtual  void             getiWeights(GTVector<TE> &ret);

virtual  GTMatrix<T>     *getStiffMatrixComp();
virtual  GTMatrix<TE>    *getStiffMatrix();
virtual  void             getStiffMatrix(GTMatrix<TE> &ret);

virtual  GTMatrix<T>     *getDerivMatrixComp(GBOOL btranspose=FALSE);
virtual  GTMatrix<TE>    *getDerivMatrix(GBOOL btranspose=FALSE);
virtual  GTMatrix<TE>    *getDerivMatrixW(GBOOL btranspose=FALSE);
virtual  GTMatrix<TE>    *getDerivMatrixiW(GBOOL btranspose=FALSE);
virtual  void             getDerivMatrix(GTMatrix<TE> &ret, GBOOL btranspose=FALSE);
virtual  void             getDerivMatrixW(GTMatrix<TE> &ret, GBOOL btranspose=FALSE);
virtual  void             getDerivMatrixiW(GTMatrix<TE> &ret, GBOOL btranspose=FALSE);

virtual  void             getLegMatrix(GTMatrix<TE> &ret);

// Evaluation methods:
virtual  TE               evalBasis(GINT  i, TE eta);
virtual  GTVector<TE>    *evalBasis(GINT  i, GTVector<TE> &eta, GTVector<TE> &vret);
virtual  GTMatrix<TE>    *evalBasis(GTVector<TE> &eta, GTMatrix<TE> &mret);
virtual  GTVector<TE>    *evalBasis(GTVector<TE> &eta, GTVector<TE> &vret);
virtual  GTMatrix<TE>    *evalBasis(TE eta[], GINT neta, GTMatrix<TE> &mret);
virtual  TE               evalDBasis(GINT i, TE eta);
virtual  GTMatrix<TE>    *evalDBasis(GTVector<TE> &eta, GTMatrix<TE> &mret);
virtual  GTMatrix<TE>    *evalDBasis(TE eta[], GINT n, GTMatrix<TE> &mret);
virtual  GTVector<TE>    *evalDBasis (GINT i, GTVector<TE> &eta, GTVector<TE> &vret);



// Protected methods:
protected:

virtual  GBOOL             computeNodes       ();
virtual  GBOOL             computeWeights     ();
virtual  GBOOL             computeDerivMatrix ();
virtual  GBOOL             computeLegendreMatrix ();
         void              computeJacobi(GINT  & , T  alpha, T  beta , T &Pn, 
                                         T &dPn  , T &Pnm1 , T &dPnm1, T &Pnm2,
                                         T &dPnm2, T &xi);
virtual  GBOOL             computeStiffMatrix();

virtual  GBOOL             init();


// Protected data:
GINT             Np_;
GINT             kstop_;
GBOOL            bInit_;
GBOOL            bNeedDerivMatrix_;
GBOOL            bNeedLegMatrix_;
GBOOL            bNeedNodes_;
GBOOL            bNeedWeights_;
T                alpha_;
T                beta_;
T                ximin_;
T                ximax_;
T                eps_;
T                ttiny_;
TE               tetiny_;
GTVector<T>      xiNodes_;
GTVector<T>      weights_;
GTVector<T>      Pn_;
GTVector<T>      dPn_;
GTMatrix<T>      Phi_;
GTMatrix<T>      dPhi_;
GTMatrix<T>      dPhiT_;
GTMatrix<T>      stiffMatrix_;
GTMatrix<T>      LegMatrix_;

// Data evaluated at the TE type:
GTVector<TE>     xiNodesEv_;
GTVector<TE>     weightsEv_;
GTVector<TE>     iweightsEv_;
GTMatrix<TE>     dPhiEv_;
GTMatrix<TE>     dPhiTEv_;
GTMatrix<TE>     dPhiWEv_;
GTMatrix<TE>     dPhiWTEv_;
GTMatrix<TE>     dPhiiWEv_;
GTMatrix<TE>     dPhiiWTEv_;
GTMatrix<TE>     stiffMatrixEv_;

};

#include "gllbasis.ipp"

#endif
