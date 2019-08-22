//==================================================================================
// Module       : ginitforce.hpp
// Date         : 7/11/19 (DLR)
// Description  : Force initialization function implementations provided
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
#include "ginitforce_user.hpp"
#include "gupdatebdy_user.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;

namespace ginitforce
{

GBOOL impl_myforceinit          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &uf);
};

#include "ginitforce.ipp"

#endif
