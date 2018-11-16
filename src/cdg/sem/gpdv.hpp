//==================================================================================
// Module       : gpdv.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the 'pdV' operator:
//                p Div u. This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none
//==================================================================================

#if !defined(_GLAPLACIANOP_HPP)
#define _GLAPLACIANOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"

class GpdV
{

public:

                          GpdV(GGrid &grid, GMass &massop);
                          GpdV(const GpdV &);
                         ~GpdV();

        void              apply(GTVector<GFTYPE> &p, GTVector<GTVector<GFTYPE>*> &u, 
                                GTVector<GFTYPE> &po);                       // Operator-field evaluation
        void              init();                                            // must call after all 'sets'
        void              set_tmp(GTVector<GTVector<GFTYPE>*> &utmp) 
                          { utmp_.resize(utmp.size()); utmp_ = utmp; }       // Set temp space 

private:
        void              def_init();
        void              reg_init();
        void              def_prod(GTVector<GFTYPE> &p, GTVector<GTVector<GFTYPE>*> &u, 
                                   GTVector<GFTYPE> &po);
        void              reg_prod(GTVector<GFTYPE> &p, GTVector<GTVector<GFTYPE>*> &u, 
                                   GTVector<GFTYPE> &po);

        GBOOL                         bInitialized_;
        GTVector<GFTYPE>              etmp1_;  // elem-based (non-global) tmp vector
        GTVector<GTVector<GFTYPE>*>   utmp_;   // global array of temp vectors
        GTVector<GTVector<GFTYPE>*>   G_;      // metric components
        GMass                        *massop_; // mass matrix, required
        GGrid                        *grid_;   // grid set on construction


};
#endif
