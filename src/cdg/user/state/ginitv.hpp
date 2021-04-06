//==================================================================================
// Module       : ginitv.hpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITV_HPP)
#define _GINITV_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gstateinfo.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gutils.hpp"


using namespace geoflow;
using namespace geoflow::tbox;



template<typename TypePack>
struct ginitv
{
        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using EqnBase    = EquationBase<Types>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using Grid       = typename Types::Grid;
        using Time       = typename Types::Time;
        using Ftype      = typename Types::Ftype;



static GBOOL impl_rand         (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

static GBOOL impl_abc_box      (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

static GBOOL impl_abc_icos     (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

static GBOOL impl_simpsum1d_box(const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

static GBOOL impl_simpsum_box  (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

static GBOOL impl_simpsum_icos (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);


};


#include "ginitv.ipp"


#endif
