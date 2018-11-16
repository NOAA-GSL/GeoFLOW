//==================================================================================
// Module       : glbasis.hpp
// Date         : 1/19/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a spectral element nodal Gauss-Legendre basis
//                Note: This basis is the Gauss-Legendre
//                basis, defined as
//                  h(xi)_j = L_N(xi) / ( dL_N(xi_j)/dxi * (xi-xi_j) )
//                where N is the order of the expansion, and L_N is the Legendre
//                polynomial of order N, and the xi_j are the nodal points. L_N
//                is obtained from the generalized Jacobi polynomial computed in
//                ComputeJacobi.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLLBasis
//==================================================================================
#if !defined(_GLBASIS_HPP)
#define _GLBASIS_HPP
#include "gtypes.h"
#include "gllbasis.hpp"


template<typename T, typename TE> // T is evaluation type; TE is compute type
class GLBasis: public GLLBasis<T,TE>
{
        static_assert(std::is_floating_point<T>::value && std::is_floating_point<TE>::value, "Illegal template type");
         using GLLBasis<T,TE>::Np_;
         using GLLBasis<T,TE>::kstop_;
         using GLLBasis<T,TE>::bInit_;
         using GLLBasis<T,TE>::bNeedDerivMatrix_;
         using GLLBasis<T,TE>::bNeedLegMatrix_;
         using GLLBasis<T,TE>::alpha_;
         using GLLBasis<T,TE>::beta_;
         using GLLBasis<T,TE>::ximin_;
         using GLLBasis<T,TE>::ximax_;
         using GLLBasis<T,TE>::eps_;
         using GLLBasis<T,TE>::ttiny_;
         using GLLBasis<T,TE>::teiny_;
         using GLLBasis<T,TE>::xiNodes_;
         using GLLBasis<T,TE>::weights_;
         using GLLBasis<T,TE>::Pn_;
         using GLLBasis<T,TE>::dPn_;
         using GLLBasis<T,TE>::Phi_;
         using GLLBasis<T,TE>::dPhi_;
         using GLLBasis<T,TE>::dPhiT_;
         using GLLBasis<T,TE>::stiffMatrix_;
         using GLLBasis<T,TE>::LegMatrix_;

         using GLLBasis<T,TE>::xiNodesEv_;
         using GLLBasis<T,TE>::weightsEv_;
         using GLLBasis<T,TE>::dPhi_Ev;
         using GLLBasis<T,TE>::dPhiTEv_;
         using GLLBasis<T,TE>::stiffMatrixEv_;

         using GLLBasis<T,TE>::init;
public:
                           GLBasis();
                           GLBasis(GINT );
                           GLBasis(GINT  , GINT );
                           GLBasis(const GLBasis &);
virtual                   ~GLBasis();

         GBOOL             resize(GINT); 
         void              operator=(const GLBasis &);
                   

// Private methods:
private:

         GBOOL             computeNodes       ();
         GBOOL             computeWeights     ();
         GBOOL             computeDerivMatrix ();
         GBOOL             computeBasisAtNodes();
         GBOOL             computeLegendreMatrix();

};
#include "glbasis.ipp"

#endif
