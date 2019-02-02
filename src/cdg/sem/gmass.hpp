//==================================================================================
// Module       : gmass.hpp
// Date         : 10/19/18 (DLR)
// Description  : Represents the SEM mass operator.  Mass-lumping is assumed.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLinOp
//==================================================================================

#if !defined(_GMASSOP_HPP)
#define _GMASSOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "glinop.hpp"

class GMass: public GLinOp
{

public:

                          GMass(GGrid &grid, GBOOL doinverse=FALSE);
                          GMass(const GMass &);
                         ~GMass();

        void              opVec_prod(GTVector<GFTYPE> &in, GTVector<GFTYPE> &out); // Operator-vector product
//      void              do_mass_lumping(GBOOL bml);                              // Set mass lumping flag
        void              init();

private:
        void              init1d();
        void              init2d();
        void              init3d();


        GBOOL             bdoinverse_;
        GBOOL             bmasslumped_;
        GTVector<GFTYPE>  mass_;
        GGrid            *grid_;


};
#endif
