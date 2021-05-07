//==================================================================================
// Module       : ginitforce_direct_user.hpp
// Date         : 7/10/19 (DLR)
// Description  : Direct user force initialization function implementations. These
//                methods are called directly during configuration, and can set 
//                forcing for entire state (v+b+s) etc. The 'component' types
//                in which component groups (v, b, s, etc) are set individually
//                are contained in separate namespaces.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITFORCE_USER_HPP)
#define _GINITFORCE_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gstateinfo.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"
#include "gutils.hpp"

using namespace geoflow;
using namespace geoflow::tbox;



template<typename TypePack>
struct ginitforce
{
        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using EqnBase    = EquationBase<Types>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using Grid       = typename Types::Grid;
        using Time       = typename Types::Time;
        using Ftype      = typename Types::Ftype;

static GBOOL impl_null    (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &t, State &utmp, State &u, State &uf);
static GBOOL impl_rand    (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &t, State &utmp, State &u, State &uf);

};

#include "ginitforce_direct_user.ipp"

#endif
