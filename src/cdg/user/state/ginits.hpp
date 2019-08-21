//==================================================================================
// Module       : ginits.hpp
// Date         : 7/10/19 (DLR)
// Description  : State initialization function implementations
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITS_HPP)
#define _GINITS_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"

using namespace geoflow;
using namespace pdeint;

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
