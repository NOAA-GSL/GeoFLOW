//==================================================================================
// Module       : glinop.hpp
// Date         : 10/19/18 (DLR)
// Description  : Represents pure abstact base class for all SEM operators
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : 
//==================================================================================

#if !defined(_GLINOP_HPP)
#define _GLINOP_HPP
#include "gtvector.hpp"
#include "ggrid.hpp"

class GLinOp
{

public:

                          GLinOp(GGrid &grid) { grid_ = &grid; bInitialized_=FALSE;}
                          GLinOp(const GLinOp &op) { grid_=op.grid_; } ;
                         ~GLinOp(){};

virtual void              opVec_prod(GTVector<GFTYPE> &in, 
                                     GTVector<GTVector<GFTYPE>*> &utmp, 
                                     GTVector<GFTYPE> &out)=0; // Operator-vector product

virtual void              init()=0; // Init after all sets, before use

protected:

        GBOOL             bInitialized_;
        GGrid            *grid_;


};
#endif
