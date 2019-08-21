//==================================================================================
// Module       : gupdatebdy.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary update function implementations provided
//                by system.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GUPDATEBDY_HPP)
#define _GUPDATEBDY_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ginitstate_factory.hpp"
#include "gupdatebdy_user.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;

namespace gupdatebdy
{

GBOOL impl_bystateinit          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
GBOOL impl_simple_outflow       (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
GBOOL impl_noslip               (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
};

#include "gupdatebdy.ipp"

#endif
