//==================================================================================
// Module       : gshapefcn_linear.hpp
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

#if !defined(_GSHAPEFCN_LINEAR_HPP)
#define _GSHAPEFCN_LINEAR_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "gshapefcn_base.hpp"

class GShapeFcn_linear: public GShapeFcn_base
{

public:

                          GShapeFcn_linear();
                          GShapeFcn_linear(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b) {gbasis_ = b;}
                          GShapeFcn_linear(const GShapeFcn_linear &);
                         ~GShapeFcn_linear();

// Methods:
        void              Ni(GTVector<GINT> &ishape,
                          GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &N); // N_i
        void              dNdXi(GTVector<GINT> &I, GINT jder,
                                GTVector<GTVector<GFTYPE>*> &xi,
                                GTVector<GFTYPE> &dNdxi);                         // dN^i/Xi^j for i, j

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

};
#endif
