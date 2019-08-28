//==================================================================================
// Module       : ginitfps.hpp
// Date         : 7/16/19 (DLR)
// Description  : Passive scalar forcing inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITFPS_HPP)
#define _GINITFPS_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "geoflow.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;


namespace ginitfps
{

GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &uf);

};



#endif
