//==================================================================================
// Module       : ginitb.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary initialization function implementations
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITB_HPP)
#define _GINITB_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"
#include "ginits_factory.hpp"


using namespace geoflow;
using namespace pdeint;

namespace ginitb
{

GBOOL impl_bystateinit          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
};

#include "ginitb.ipp"

#endif
