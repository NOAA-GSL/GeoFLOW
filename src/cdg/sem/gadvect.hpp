//==================================================================================
// Module       : gadvect.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the advection operator:
//                u.Grad p  This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none
//==================================================================================

#if !defined(_GADVECTOP_HPP)
#define _GADVECTOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"

class GAdvect
{

public:

                          GAdvect(GGrid &grid);
                          GAdvect(const GAdvect &);
                         ~GAdvect();

        void              apply(GTVector<GFTYPE> &p, const GTVector<GTVector<GFTYPE>*> &u, 
                                GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GFTYPE> &po);                       // Operator-field evaluation
        void              init();                                            // must call after all 'sets'

private:
        void              def_init();
        void              reg_init();
        void              def_prod(GTVector<GFTYPE> &p, const GTVector<GTVector<GFTYPE>*> &u, 
                                   GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GFTYPE> &po);
        void              reg_prod(GTVector<GFTYPE> &p, const GTVector<GTVector<GFTYPE>*> &u, 
                                   GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GFTYPE> &po);

        GBOOL                         bInitialized_;
        GTVector<GFTYPE>              etmp1_;  // elem-based (non-global) tmp vector
        GTVector<GTVector<GFTYPE>*>   G_;      // metric components
        GGrid                        *grid_;   // grid set on construction


};
#endif
