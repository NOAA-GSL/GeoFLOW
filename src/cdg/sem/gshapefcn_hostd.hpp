//==================================================================================
// Module       : gshapefcn_hostd.hpp
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

#if !defined(_GSHAPEFCN_HOSTD_HPP)
#define _GSHAPEFCN_HOSTD_HPP
#include "gtvector.hpp"
#include "gshapefcn_base.hpp"

class GShapeFcn_hostd: public GShapeFcn_base
{

public:

                          GShapeFcn_hostd();
                          GShapeFcn_hostd(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b);
                          GShapeFcn_hostd(const GShapeFcn_hostd &);
                         ~GShapeFcn_hostd();

         void             operator=(const GShapeFcn_hostd &obj) 
                          { gbasis_ = obj.gbasis_; }
                   
// Methods:
        void              Ni(GTVector<GINT> &I,
                             GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N); // N^i
        void              dNdXi(GTVector<GINT> &I, GINT j, 
                                GTVector<GTVector<GFTYPE>*> &xi,
                                GTVector<GFTYPE> &out);   // dN^i/Xi^j for i, j

private:
        void              Ni_1d(GTVector<GINT> &I,
                             GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N); 
        void              Ni_2d(GTVector<GINT> &I,
                             GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N); 
        void              Ni_3d(GTVector<GINT> &I,
                             GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N); 
        void              dNdXi_1d(GTVector<GINT> &I, GINT j, 
                                                   GTVector<GTVector<GFTYPE>*> &xi,
                                                   GTVector<GFTYPE> &out); 
        void              dNdXi_2d(GTVector<GINT> &I, GINT j,
                                                   GTVector<GTVector<GFTYPE>*> &xi,
                                                   GTVector<GFTYPE> &out); 
        void              dNdXi_3d(GTVector<GINT> &I, GINT j,
                                                   GTVector<GTVector<GFTYPE>*> &xi,
                                                   GTVector<GFTYPE> &out); 

        GTVector<GTVector<GFTYPE>>   d_; // tmp vectors for multiplications


};
#endif
