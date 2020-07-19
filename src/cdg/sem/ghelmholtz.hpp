//==================================================================================
// Module       : ghelmholtz.hpp
// Date         : 11/1/18 (DLR)
// Description  : Represents the SEM generalized Helmholtz operator:
//                  H = qM + pL,
//                where M = mass operator, L is Laplacian operator, and
//                q, p are scalars that may or may not be constant. 
//                The mass term is added only if calls are made to 
//                set_mass_scalar and p is applied only if set_Lap_scalar, 
//                respectively. In fact, only if the mass operator is set
//                is this operator a real Helmholtz operator; otherwise, it's
//                really just a weak Laplacian operator. 
//                Note: this operator will fail if the grid contains more 
//                      than a single element type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLinOp
//==================================================================================

#if !defined(_GHELMHOLTZOP_HPP)
#define _GHELMHOLTZOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "gmass.hpp"
#include "ggrid.hpp"
#include "glinop.hpp"

class GHelmholtz: public GLinOp
{

public:

                          GHelmholtz(GGrid &grid);
                          GHelmholtz(const GHelmholtz &);
                         ~GHelmholtz();

        void              opVec_prod(GTVector<GFTYPE> &in, 
                                     GTVector<GTVector<GFTYPE>*> &utmp,
                                     GTVector<GFTYPE> &out);                  // Operator-vector product 
        void              set_Lap_scalar(GTVector<GFTYPE> &p);                // Scalar multipliying Laplacian
        void              set_mass_scalar(GTVector<GFTYPE> &q);               // Scalar multiplying Mass
        void              init();                                             // must call after all 'sets'
        void              use_metric(GBOOL flag) {buse_metric_ = flag;}       // set flag to use metric & Jacobian
//      void              set_tmp(GTVector<GTVector<GFTYPE>*> &utmp) 
//                        { utmp_.resize(utmp.size()); utmp_ = utmp; }       // Set temp space 

private:
        void              def_init();
        void              reg_init();
        void              def_prod(GTVector<GFTYPE> &in, 
                                   GTVector<GTVector<GFTYPE>*> &utmp,
                                   GTVector<GFTYPE> &out);
        void              reg_prod(GTVector<GFTYPE> &in, 
                                   GTVector<GTVector<GFTYPE>*> &utmp,
                                   GTVector<GFTYPE> &out);
        void              embed_prod(GTVector<GFTYPE> &in, 
                                     GTVector<GTVector<GFTYPE>*> &utmp,
                                     GTVector<GFTYPE> &out);
        void              compute_refderivs(GTVector<GFTYPE> &, 
                                            GTVector<GTVector<GFTYPE>*> &, GBOOL btrans=FALSE);
        void              compute_refderivsW(GTVector<GFTYPE> &, 
                                             GTVector<GTVector<GFTYPE>*> &, GBOOL btrans=FALSE);

        void              compute_div(GTVector<GTVector<GFTYPE>*> &, 
                                      GTVector<GFTYPE> &, GBOOL btrans=TRUE); 

        GBOOL                         buse_metric_;   // use metric terms?
        GBOOL                         bown_q_;        // does object own q?
        GBOOL                         bown_p_;        // does object own p?
        GBOOL                         bown_mass_;     // does object own massop?
        GBOOL                         bcompute_helm_; // compute full Helm, not just Lap?
        GTVector<GFTYPE>             *p_;    // scalar multiplying L
        GTVector<GFTYPE>             *q_;    // scalar multiplying M
        GTVector<GFTYPE>              etmp1_;// elem-based (not global) tmp vector
        GTVector<GTVector<GFTYPE>*>   utmp_; // global array of temp vectors
        GTMatrix<GTVector<GFTYPE>*>   G_;    // metric components
        GGrid                        *grid_; // grid set on construction


};
#endif
