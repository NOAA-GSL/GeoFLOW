//==================================================================================
// Module       : ginitb.hpp
// Date         : 7/16/19 (DLR)
// Description  : Mag. field inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITB_HPP)
#define _GINITB_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gstateinfo.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gutils.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;



namespace ginitb
{

GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u);

};



#endif
