//==================================================================================
// Module       : ginitv.hpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITV_HPP)
#define _GINITV_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "geoflow.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;


namespace ginitv
{

GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
GBOOL impl_abc_box   (const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
GBOOL impl_abc_icos  (const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u);

};



#endif
