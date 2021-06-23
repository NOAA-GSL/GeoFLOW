//==================================================================================
// Module       : ginitstate_user.hpp
// Date         : 7/10/19 (DLR)
// Description  : Direct user state initialization function implementations. These
//                methods are called directly during configuration, and can set 
//                forcing for entire state (v+b+s) etc. The 'component' types
//                in which component groups (v, b, s, etc) are set individually
//                are contained in separate namespaces.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITSTATE_USER_HPP)
#define _GINITSTATE_USER_HPP

#include "gtypes.h"
#include "gphysics.h"
#include "gtvector.hpp"
#include "gcomm.hpp"
#include "gutils.hpp"
#include "gmconv.hpp"
#include "gburgers.hpp"
#include "gmtk.hpp"
#include "ginitstate_direct_user.hpp"
#include <random>

#include "tbox/property_tree.hpp"

using namespace geoflow;
using namespace geoflow::tbox;




template<typename TypePack>
struct ginitstate
{
        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using EqnBase    = EquationBase<Types>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using Grid       = typename Types::Grid;
        using GridBox    = typename Types::GridBox;
        using GridIcos   = typename Types::GridIcos;
        using Time       = typename Types::Time;
        using Ftype      = typename Types::Ftype;

static GBOOL impl_boxnwaveburgers      (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_icosnwaveburgers     (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxdirgauss          (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxpergauss          (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_icosgauss            (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxdrywarmbubble     (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxdrybubble         (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxsod               (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxdryscharadv       (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_icosabcconv          (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
static GBOOL impl_boxmtnwave           (const PropertyTree& ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);


};


#include "ginitstate_direct_user.ipp"


#endif
