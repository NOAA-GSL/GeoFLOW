//==================================================================================
// Module       : ginits.hpp
// Date         : 7/16/19 (DLR)
// Description  : Scalar inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITS_HPP)
#define _GINITS_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gstateinfo.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gutils.hpp"


using namespace geoflow;
using namespace geoflow::tbox;




template<typename TypePack>
struct ginits 
{
        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using EqnBase    = EquationBase<Equation>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using Grid       = typename Types::Grid;
        using Time       = typename Types::Time;
        using Ftype      = typename Types::Ftype;

static GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, GGrid &grid, Time &time, State &utmp, State &u);


};


#include "inits.ipp"





#endif
