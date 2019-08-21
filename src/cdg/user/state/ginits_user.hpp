//==================================================================================
// Module       : ginits_user.hpp
// Date         : 7/10/19 (DLR)
// Description  : User state initialization function implementations
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITS_USER_HPP)
#define _GINITS_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"

using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;

namespace ginits
{

GBOOL impl_boxnwaveburgers      (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
GBOOL impl_boxdirgauss          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
GBOOL impl_boxpergauss          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
GBOOL impl_icosgauss            (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);
};

#include "ginits_burgers.ipp"

#endif
