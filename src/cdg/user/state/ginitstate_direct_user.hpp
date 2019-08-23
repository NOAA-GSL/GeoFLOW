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

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"

using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;

namespace ginitstate
{

GBOOL impl_boxnwaveburgers      (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
GBOOL impl_boxdirgauss          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
GBOOL impl_boxpergauss          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
GBOOL impl_icosgauss            (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
};

#include "ginitstate_direct_user.ipp"

#endif
