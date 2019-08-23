//==================================================================================
// Module       : ginitfv.hpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity force inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITFV_HPP)
#define _GINITFV_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "geoflow.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;


namespace ginitfv
{

GBOOL impl_rand      (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &uf);

};

#include "ginitfv.ipp"

#endif
