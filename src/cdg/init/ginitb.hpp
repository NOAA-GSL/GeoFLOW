//==================================================================================
// Module       : ginitb.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary initialization function implementations provided
//                by system.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITB_HPP)
#define _GINITB_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ginitb_user.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

namespace ginitb
{

GBOOL impl_mystateinit          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
};

#include "ginitb.ipp"

#endif
