//==================================================================================
// Module       : ginitstate.hpp
// Date         : 7/11/19 (DLR)
// Description  : State initialization function implementations provided
//                by system.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITSTATE_HPP)
#define _GINITSTATE_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gupdatebdy.hpp"
#include "ginitstate_user.hpp"
#include "gupdatebdy_user.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;

namespace ginitstate
{

GBOOL impl_mystateinit          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
};

#include "ginitstate.ipp"

#endif
