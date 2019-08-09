//==================================================================================
// Module       : gspecb.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary specification function implementations provided
//                by system.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECB_HPP)
#define _GSPECB_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"
#include "gspecb_user.hpp"


using namespace geoflow;
using namespace pdeint;

namespace gspecb
{

GBOOL impl_bystateinit          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
GBOOL impl_simple_outflow       (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
GBOOL impl_noslip               (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
};

#include "gspecb.ipp"

#endif
