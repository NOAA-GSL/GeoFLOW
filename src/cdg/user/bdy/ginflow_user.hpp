//==================================================================================
// Module       : gns_inflow_user.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary initialization function implementations provided
//                by user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITBDY_USER_HPP)
#define _GINITBDY_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gstateinfo.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::tbox;



template<typename TypePack>
struct GInflowUser
{
public:
        using Types      = TypePack;
        using Base       = UpdateBdyBase<Types>;
        using EqnBase    = EquationBase<TypePack>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using State      = typename Types::State;
        using Grid       = typename Types::Grid;
        using GridBox    = typename Types::GridBox;
        using GridIcos   = typename Types::GridIcos;
        using Ftype      = typename Types::Ftype;
        using Time       = typename Types::Time;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");

static GBOOL myinflow  (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, const GINT id, State &u, State &utmp, State &ub);

};


#include "ginflow_user.ipp"



#endif
