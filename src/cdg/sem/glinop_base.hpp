//==================================================================================
// Module       : glinop_base.hpp
// Date         : 10/19/18 (DLR)
// Description  : Represents pure abstact base class for all SEM operators
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : 
//==================================================================================

#if !defined(_GLINOP_BASE_HPP)
#define _GLINOP_BASE_HPP

#include "gtvector.hpp"


template<typename TypePack>
class GLinOpBase
{
public:
                          using Types      = TypePack;
                          using State      = typename Types::State;
                          using StateComp  = typename Types::StateComp;
                          using Grid       = typename Types::Grid;
                          using Mass       = typename Types::Mass;
                          using Ftype      = typename Types::Ftype;
                          using Size       = typename Types::Size;


                          GLinOpBase(Grid &grid) { grid_ = &grid; bInitialized_=FALSE;}
                          GLinOpBase(const GLinOpBase &op) { grid_=op.grid_; } ;
                         ~GLinOpBase(){};

virtual void              opVec_prod(StateComp  &in, 
                                     State      &utmp, 
                                     StateComp  &out)=0; // Operator-vector product

virtual void              init()=0; // Init after all sets, before use

protected:

        GBOOL             bInitialized_;
        Grid             *grid_;

};
#endif
