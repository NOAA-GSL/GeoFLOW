//==================================================================================
// Module       : gnbasis.hpp
// Date         : 1/19/18 (DLR)
// Description  : Forms pure virtual abstract base class for all 
//                allowed nodal basis objects
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GNBasis
//==================================================================================

#if !defined(_GNBASIS_HPP)
#define _GNBASIS_HPP
#include "gtypes.h"
#include "gtvector.hpp"
#include "gtmatrix.hpp"

// 2 template args: T is the type used for computations; 
// TE, for evaluations (deep copies, etc.). Must
// have sizeof(T) >= sizeof(TE).
template<typename T, typename TE> 
class GNBasis
{
        static_assert(std::is_floating_point<T>::value,"Illegal template type");
        static_assert(sizeof(T)>=sizeof(TE),"Invalid computational template type");

public:

                          GNBasis(){};
                          GNBasis(GINT  inorder){};
                          GNBasis(GINT  inorder, GINT  maxorder){};
                          GNBasis(const GNBasis &){};
virtual                  ~GNBasis(){};

virtual  void             operator=(const GNBasis &){}
virtual  GBOOL            resize(GINT  order)=0;

                   
// Data retrieval methods:
virtual  T                getXimin ()=0;
virtual  T                getXimax ()=0; 
virtual  GINT             getOrder()=0;
virtual  void             setOrder(GINT)=0;

virtual  GTVector<T>     *getXiNodesComp()=0;
virtual  GTVector<TE>    *getXiNodes()=0;
virtual  TE              *getXiNodes(TE *ret, GINT  num)=0;
virtual  void             getXiNodes(GTVector<TE> &ret)=0;

virtual  GTVector<T>     *getWeightsComp()=0;
virtual  GTVector<TE>    *getWeights()=0;
virtual  TE              *getWeights(TE *ret, GINT num)=0;
virtual  void             getWeights(GTVector<TE> &ret)=0;

virtual  GTVector<TE>    *getiWeights()=0;
virtual  void             getiWeights(GTVector<TE> &ret)=0;

virtual  GTMatrix<T>     *getStiffMatrixComp()=0;
virtual  GTMatrix<TE>    *getStiffMatrix()=0;
virtual  void             getStiffMatrix(GTMatrix<TE> &ret)=0;

virtual  GTMatrix<T>     *getDerivMatrixComp(GBOOL btrans=FALSE)=0;
virtual  GTMatrix<TE>    *getDerivMatrix(GBOOL btrans=FALSE)=0;
virtual  GTMatrix<TE>    *getDerivMatrixW(GBOOL btrans=FALSE)=0;
virtual  GTMatrix<TE>    *getDerivMatrixiW(GBOOL btrans=FALSE)=0;
virtual  void             getDerivMatrix(GTMatrix<TE> &ret, GBOOL btrans=FALSE)=0;
virtual  void             getDerivMatrixW(GTMatrix<TE> &ret, GBOOL btrans=FALSE)=0;
virtual  void             getDerivMatrixiW(GTMatrix<TE> &ret, GBOOL btrans=FALSE)=0;

virtual  void             getLegMatrix(GTMatrix<TE> &ret)=0;

// Evaluation methods:
virtual  TE               evalBasis(GINT  i, TE eta)=0;
virtual  GTVector<TE>    *evalBasis(GINT  i, GTVector<TE> &eta, GTVector<TE> &vret)=0;
virtual  GTMatrix<TE>    *evalBasis(GTVector<TE> &eta, GTMatrix<TE> &mret)=0;
virtual  GTMatrix<TE>    *evalBasis(TE eta[], GINT neta, GTMatrix<TE> &mret)=0;
virtual  TE               evalDBasis(GINT i, TE eta)=0;
virtual  GTMatrix<TE>    *evalDBasis(GTVector<TE> &eta, GTMatrix<TE> &mret)=0;
virtual  GTMatrix<TE>    *evalDBasis(TE eta[], GINT n, GTMatrix<TE> &mret)=0;
virtual  GTVector<TE>    *evalDBasis (GINT i, GTVector<TE> &eta, GTVector<TE> &vret)=0;


};
#endif
